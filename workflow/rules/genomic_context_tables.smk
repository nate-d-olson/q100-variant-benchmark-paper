rule genomic_context_coverage:
    """Compute base-level coverage for genomic context regions against benchmark."""
    input:
        strat_bed="resources/stratifications/{ref}_{genomic_context}.bed.gz",
        bench_bed="resources/benchmarksets/{bench_version}_{ref}_{bench_type}_benchmark.bed",
    output:
        cov_bed="results/genomic-context_coverage/{bench_version}_{ref}_{bench_type}/{genomic_context}_cov.bed",
    log:
        "logs/genomic_context_coverage/{bench_version}_{ref}_{bench_type}/{ref}_{genomic_context}.log",
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        bedtools coverage -a {input.strat_bed} -b {input.bench_bed} 1> {output.cov_bed} 2> {log}
        """
