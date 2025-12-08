"""
VCF processing rules for variant characterization
"""


rule generate_sv_len:
    input:
        vcf=config["vcf_files"]["chm13_stvar"],
        bed=config["bed_files"]["chm13_stvar"],
    output:
        tsv=config["outputs"]["sv_len"],
    log:
        "logs/vcf_processing/sv_len.log",
    threads: 1
    resources:
        mem_mb=4096,
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        bcftools query -HH \
            -f '%CHROM\\t%POS0\\t%END\\t%SVTYPE\\t%VKX\\t%SVLEN\\t%TRF\\t%TRFstart\\t%TRFend\\t%RM_clsfam\\n' \
            -R {input.bed} {input.vcf} \
            > {output.tsv} 2> {log}
        
        if [ ! -s {output.tsv} ]; then
            echo "ERROR: Output file is empty" >> {log}
            exit 1
        fi
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
            "{subset_dir}/{{benchmark}}.vcf.gz".format(
                subset_dir=config["outputs"]["subset_vcf_dir"]
            ),
            non_empty=True,
        ),
    params:
        bed=lambda wildcards: config["benchmarksets"][wildcards.benchmark]["bed"],
    log:
        "logs/vcf_processing/{benchmark}_subset.log",
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
    Generate RTG vcfstats for a benchmark set.
    
    Calculates comprehensive variant statistics using RTG Tools vcfstats
    on VCFs that have been filtered to benchmark regions.
    """
    input:
        vcf="{subset_dir}/{{benchmark}}.vcf.gz".format(
            subset_dir=config["outputs"]["subset_vcf_dir"]
        ),
        tbi="{subset_dir}/{{benchmark}}.vcf.gz.tbi".format(
            subset_dir=config["outputs"]["subset_vcf_dir"]
        ),
    output:
        stats=ensure(
            "{rtg_dir}/{{benchmark}}_vcfstats.txt".format(
                rtg_dir=config["outputs"]["rtg_stats_dir"]
            ),
            non_empty=True,
        ),
    log:
        "logs/vcf_processing/{benchmark}_rtg_vcfstats.log",
    threads: 1
    resources:
        mem_mb=4096,
    conda:
        "../envs/rtg-tools.yaml"
    shell:
        """
        # Log execution start
        echo "Running RTG vcfstats for {wildcards.benchmark}" > {log}
        echo "Started at $(date)" >> {log}
        
        # Run RTG vcfstats
        rtg vcfstats {input.vcf} > {output.stats} 2>> {log}
        
        # Log completion
        echo "Completed at $(date)" >> {log}
        echo "Output size: $(wc -l < {output.stats}) lines" >> {log}
        """
