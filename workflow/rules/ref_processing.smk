"""
Reference fasta processing rules
"""


rule get_ref_assembled_bases:
    """
    Calculating the number of non-N bases in reference genomes
    """
    input:
        fa="resources/references/{ref}.fa.gz",
        fai="resources/references/{ref}.fa.gz.fai",
    output:
        "results/ref_genome_sizes/{ref}_size.tsv",
    log:
        "logs/ref_genome_sizes/{ref}_size.log",
    threads: 2
    resources:
        mem_mb=2048,
    conda:
        "../envs/samtools.yaml"
    params:
        chroms=get_chromosomes,
    shell:
        """
        samtools faidx {input.fa} {params.chroms} |
            seqkit fx2tab -inlH -j {threads} \
            -C ACGT -C N \
            --skip-file-check \
            -o {output} > {log} 2>&1
        """
