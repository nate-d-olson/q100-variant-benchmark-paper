"""
VCF processing rules for variant characterization
"""


rule index_vcf:
    """
    Index VCF files using bcftools index.

    Generic rule to create tabix index (.tbi) for any compressed VCF file.
    Uses bcftools for standardized, reproducible indexing.
    """
    input:
        "{prefix}.vcf.gz",
    output:
        "{prefix}.vcf.gz.tbi",
    log:
        "logs/vcf_processing/{prefix}_index.log",
    threads: 1
    resources:
        mem_mb=2048,
    conda:
        "../envs/bcftools.yaml"
    params:
        extra="-t",
    shell:
        "bcftools index --threads {threads} {params.extra} {input} > {log} 2>&1"
