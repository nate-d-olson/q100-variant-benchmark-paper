"""
VCF processing rules for variant characterization
"""


rule download_stratification:
    """
    Download GIAB genome stratification BED files.
    
    Downloads stratification BED files from GIAB FTP for genomic contexts
    (tandem repeats, homopolymers, segmental duplications, low mappability).
    Files are downloaded once per reference genome (GRCh37, GRCh38, CHM13)
    and shared across all benchmarks using that reference.
    
    Stratifications: v3.6 from GIAB
    """
    output:
        bed="resources/stratifications/{ref}_{context}.bed.gz",
    params:
        url=lambda wildcards: (
            f"https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/"
            f"genome-stratifications/v3.6/{wildcards.ref}@all/"
            f"{CONTEXT_PATHS[wildcards.context].format(ref=wildcards.ref)}"
        ),
    log:
        "logs/downloads/stratifications/{ref}_{context}.log",
    message:
        "Downloading {wildcards.ref} {wildcards.context} stratification BED"
    retries: 3
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        # Log download start
        echo "Downloading stratification BED: {wildcards.ref} {wildcards.context}" > {log}
        echo "URL: {params.url}" >> {log}
        echo "Started at $(date)" >> {log}
        
        # Download BED file
        wget -O {output.bed} {params.url} 2>> {log}
        
        # Verify download
        if [ ! -s {output.bed} ]; then
            echo "ERROR: Downloaded file is empty" >> {log}
            exit 1
        fi
        
        # Log completion
        echo "Completed at $(date)" >> {log}
        echo "File size: $(du -h {output.bed} | cut -f1)" >> {log}
        """


rule generate_svlen:
    """
    Generate SV length table for structural variant benchmark sets.
    
    Extracts variant information including chromosome, position, end, SV type,
    VKX annotation, SV length, tandem repeat features, and repeat masker classification
    from VCF files filtered to benchmark regions.
    """
    input:
        vcf="results/subset_vcfs/{benchmark}.vcf.gz",
        tbi="results/subset_vcfs/{benchmark}.vcf.gz.tbi",
    output:
        tsv=ensure(
            "results/sv_len/{benchmark}_svlen.tsv",
            non_empty=True,
        ),
    log:
        "logs/sv_len/{benchmark}_svlen.log",
    message:
        "Extracting SV lengths for {wildcards.benchmark}"
    threads: 1
    resources:
        mem_mb=4096,
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        # Log execution start
        echo "Generating SV length table for {wildcards.benchmark}" > {log}
        echo "Started at $(date)" >> {log}
        
        # Extract SV information from VCF
        bcftools query -HH \
            -f '%CHROM\\t%POS0\\t%END\\t%SVTYPE\\t%VKX\\t%SVLEN\\t%TRF\\t%TRFstart\\t%TRFend\\t%RM_clsfam\\n' \
            {input.vcf} \
            > {output.tsv} 2>> {log}
        
        # Verify output is not empty
        if [ ! -s {output.tsv} ]; then
            echo "ERROR: Output file is empty" >> {log}
            exit 1
        fi
        
        # Log completion
        echo "Completed at $(date)" >> {log}
        echo "Output size: $(wc -l < {output.tsv}) lines" >> {log}
        """


rule index_vcf:
    """
    Index VCF files using bcftools index.
    
    Generic rule to create tabix index (.tbi) for any compressed VCF file.
    Uses the Snakemake wrapper for standardized, reproducible indexing.
    """
    input:
        "{prefix}.vcf.gz"
    output:
        "{prefix}.vcf.gz.tbi"
    log:
        "logs/vcf_processing/{prefix}_index.log"
    threads: 1
    resources:
        mem_mb=2048
    conda:
        "../envs/bcftools.yaml"
    params:
        extra="-t"  # Optional bcftools index parameters
    shell:
        "bcftools index --threads {threads} {params.extra} {input} > {log} 2>&1"


rule subset_vcf_to_benchmark_regions:
    """
    Subset VCF to benchmark regions for historical benchmarks.
    
    For benchmark sets that have separate VCF and BED files,
    filter the VCF to only include variants within benchmark regions.
    Draft benchmarks that are already filtered are simply symlinked.
    Index is created separately by the index_vcf rule.
    """
    input:
        vcf=lambda wildcards: config["benchmarksets"][wildcards.benchmark]["vcf"],
        tbi=lambda wildcards: config["benchmarksets"][wildcards.benchmark]["vcf"] + ".tbi",
    output:
        vcf=ensure(
            "results/subset_vcfs/{benchmark}.vcf.gz",
            non_empty=True,
        ),
    params:
        bed=lambda wildcards: config["benchmarksets"][wildcards.benchmark]["bed"],
    log:
        "logs/subset_vcfs/{benchmark}_subset.log",
    message:
        "Filtering {wildcards.benchmark} to benchmark regions"
    threads: 2
    resources:
        mem_mb=4096,
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        # Log execution start
        echo "Subsetting {wildcards.benchmark} to benchmark regions" > {log}
        echo "Started at $(date)" >> {log}
        
        if [ -z "{params.bed}" ] || [ "{params.bed}" = "None" ]; then
            # No BED file - VCF is already filtered, create symlink
            echo "No BED file provided, creating symlink to original VCF" >> {log}
            ln -sf $(realpath {input.vcf}) {output.vcf}
        else
            # Filter VCF to benchmark regions
            echo "Filtering VCF to benchmark regions" >> {log}
            bcftools view -R {params.bed} {input.vcf} -Oz -o {output.vcf} 2>> {log}
        fi
        
        # Log completion
        echo "Completed at $(date)" >> {log}
        echo "Output: {output.vcf}" >> {log}
        """


rule rtg_vcfstats:
    """
    Generate per-chromosome RTG vcfstats for a benchmark set.
    
    Calculates comprehensive variant statistics using RTG Tools vcfstats
    on VCFs that have been filtered to benchmark regions, for a single chromosome.
    
    Chromosome naming:
    - CHM13/GRCh38: chr1-chr22, chrX, chrY
    - GRCh37: 1-22, X, Y
    """
    input:
        vcf="results/subset_vcfs/{benchmark}.vcf.gz",
        tbi="results/subset_vcfs/{benchmark}.vcf.gz.tbi",
    output:
        stats=ensure(
            "results/vcfstats/{benchmark}_{chrom}_vcfstats.txt",
            non_empty=True,
        ),
    log:
        "logs/vcfstats/{benchmark}_{chrom}_vcfstats.log",
    message:
        "Computing vcfstats for {wildcards.benchmark} chr {wildcards.chrom}"
    threads: 1
    resources:
        mem_mb=2048,
    conda:
        "../envs/rtg-tools.yaml"
    shell:
        """
        # Log execution start
        echo "Running RTG vcfstats for {wildcards.benchmark} chromosome {wildcards.chrom}" > {log}
        echo "Started at $(date)" >> {log}
        
        # Run RTG vcfstats for specific chromosome
        rtg vcfsubset --region={wildcards.chrom} -i {input.vcf} -o - |
            rtg vcfstats - \
            > {output.stats} 2>> {log}
        
        # Log completion
        echo "Completed at $(date)" >> {log}
        echo "Output size: $(wc -l < {output.stats}) lines" >> {log}
        """


rule intersect_vcf_with_context:
    """
    Intersect benchmark VCF with genomic stratification contexts.
    
    Filters benchmark VCF to only variants within specific genomic contexts
    (tandem repeats, homopolymers, segmental duplications, low mappability).
    Uses bcftools view with -R to restrict to stratification BED regions.
    
    Note: Output marked for potential temp() after validation to save disk space.
    """
    input:
        vcf="results/subset_vcfs/{benchmark}.vcf.gz",
        tbi="results/subset_vcfs/{benchmark}.vcf.gz.tbi",
        bed=lambda wildcards: (
            f"resources/stratifications/"
            f"{get_reference_for_benchmark(wildcards.benchmark)}_{wildcards.context}.bed.gz"
        ),
    output:
        vcf="results/context_vcfs/{benchmark}_{context}.vcf.gz",  # TODO: add temp() after validation
    log:
        "logs/context_vcfs/{benchmark}_{context}.log",
    threads: 1
    resources:
        mem_mb=4096,
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        # Log execution start
        echo "Intersecting {wildcards.benchmark} with {wildcards.context} context" > {log}
        echo "Started at $(date)" >> {log}
        
        # Filter VCF to stratification regions
        # Note: -R option allows empty output if no variants overlap with regions
        bcftools view -R {input.bed} {input.vcf} -Oz -o {output.vcf} 2>> {log}
        
        # Check if output contains variants
        VARIANT_COUNT=$(bcftools view -H {output.vcf} | wc -l)
        echo "Variants in context: $VARIANT_COUNT" >> {log}
        
        if [ "$VARIANT_COUNT" -eq 0 ]; then
            echo "WARNING: No variants found in {wildcards.context} context for {wildcards.benchmark}" >> {log}
        fi
        
        # Log completion
        echo "Completed at $(date)" >> {log}
        echo "Output: {output.vcf}" >> {log}
        """


rule count_context_variants:
    """
    Count variants by type within genomic stratification contexts.
    
    Extracts variant types using bcftools query and counts occurrences.
    Uses TYPE field for small variants (snp, ins, del, indel, mnp)
    and SVTYPE field for structural variants (DEL, INS, DUP, INV, BND).
    
    Outputs long-format TSV with columns:
    - benchmark: Benchmark set name
    - ref: Reference genome
    - context: Genomic context (TR, HP, SD, etc.)
    - variant_type: Type of variant
    - count: Number of variants of that type
    """
    input:
        vcf="results/context_vcfs/{benchmark}_{context}.vcf.gz",
    output:
        counts="results/context_counts/{benchmark}_{context}_counts.tsv",
    params:
        ref=lambda wildcards: get_reference_for_benchmark(wildcards.benchmark),
        script="workflow/scripts/count_variants_by_type.py",
    log:
        "logs/context_counts/{benchmark}_{context}.log",
    message:
        "Counting variants for {wildcards.benchmark} in {wildcards.context}"
    threads: 1
    resources:
        mem_mb=2048,
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        # Log execution start
        echo "Counting variants for {wildcards.benchmark} in {wildcards.context} context" > {log}
        echo "Started at $(date)" >> {log}
        
        # Determine query field based on benchmark type
        if [[ "{wildcards.benchmark}" == *"stvar"* ]]; then
            QUERY_FIELD="%SVTYPE"
            echo "Using SVTYPE field for structural variants" >> {log}
        else
            QUERY_FIELD="%TYPE"
            echo "Using TYPE field for small variants" >> {log}
        fi
        
        # Check if VCF has any variants
        VARIANT_COUNT=$(bcftools view -H {input.vcf} | wc -l)
        echo "Total variants to count: $VARIANT_COUNT" >> {log}
        
        if [ "$VARIANT_COUNT" -eq 0 ]; then
            echo "No variants in context - creating empty count file" >> {log}
            # Python script handles empty input and outputs appropriate zero counts
            python {params.script} {wildcards.benchmark} {params.ref} {wildcards.context} \
                < /dev/null > {output.counts} 2>> {log}
        else
            # Extract variant types and count
            bcftools query -f "$QUERY_FIELD\\n" {input.vcf} 2>> {log} | \
                python {params.script} {wildcards.benchmark} {params.ref} {wildcards.context} \
                > {output.counts} 2>> {log}
        fi
        
        # Log completion
        echo "Completed at $(date)" >> {log}
        echo "Output lines: $(wc -l < {output.counts})" >> {log}
        """

