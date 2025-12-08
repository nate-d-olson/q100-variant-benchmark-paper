"""
VCF processing rules for variant characterization
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
