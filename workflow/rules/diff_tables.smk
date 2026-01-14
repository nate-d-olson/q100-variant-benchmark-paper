rule diff_region_coverage:
    """Combine stratification BEDs with IDs for bcftools annotation."""
    input:
        strat_bed="resources/stratifications/{ref}_{strat}.bed.gz",
        bench_bed="resources/benchmarksets/{bench_version}_{ref}_{var_type}_benchmark.bed",
    output:
        cov_bed="results/diff_region_coverage/{bench_version}_{ref}_{var_type}/{strat}_cov.bed",
    log:
        "logs/diff_region_coverage/{bench_version}_{ref}_{var_type}/{ref}_{strat}.log",
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        bedtools coverage -a {input.strat_bed} -b {input.bench_bed} 1> {output.cov_bed} 2> {log}
        """
