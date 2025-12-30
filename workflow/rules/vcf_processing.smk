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


rule split_multiallelics:
    """Split multiallelic variants BEFORE annotation."""
    input:
        vcf="resources/benchmarksets/{benchmark}_benchmark.vcf.gz",
        vcfidx="resources/benchmarksets/{benchmark}_benchmark.vcf.gz.tbi",
        ref=lambda w: f"resources/references/{config["benchmarksets"][w.benchmark].get("ref")}.fa.gz",
        refidx=lambda w: f"resources/references/{config["benchmarksets"][w.benchmark].get("ref")}.fa.gz.fai",
    output:
        vcf="results/split_multiallelics/{benchmark}/split.vcf.gz",
    params:
        old_rec_tag=lambda w: "--old-rec-tag SVLEN" if "stvar" in w.benchmark else "",
    log:
        "logs/split_multiallelics/{benchmark}.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        echo "Splitting multiallelics for {wildcards.benchmark}" > {log}

        bcftools norm --multiallelics -any \
            --check-ref w \
            {params.old_rec_tag} \
            --fasta-ref {input.ref} \
            {input.vcf} \
            | bcftools sort -Oz -o {output.vcf} 2>> {log}

        echo "Completed at $(date)" >> {log}
        """
