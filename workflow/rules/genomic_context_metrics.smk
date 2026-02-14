"""
Rules for computing genomic context (difficult region) coverage metrics.

Summarizes bedtools coverage output (_cov.bed files from genomic_context_tables.smk)
into a single CSV table with per-context overlap metrics:
- context_name: Name of the genomic context region
- context_bp: Total bases in the genomic context BED
- intersect_bp: Bases overlapping between genomic context and benchmark regions
- pct_of_context: Percent of genomic context covered by benchmark
- pct_of_bench: Percent of benchmark regions covered by genomic context
"""


rule compute_genomic_context_coverage_table:
    """
    Summarize bedtools coverage output into a genomic context coverage table.

    Reads _cov.bed files produced by the genomic_context_coverage rule and
    the benchmark BED to compute per-context overlap metrics.
    """
    input:
        cov_beds=get_genomic_context_cov_beds,
        bench_bed="resources/benchmarksets/{benchmark}_benchmark.bed",
    output:
        csv=ensure(
            "results/genomic_context/{benchmark}/genomic_context_coverage_table.csv",
            non_empty=True,
        ),
    params:
        context_names=lambda wc: get_genomic_context_ids(wc),
    log:
        "logs/compute_genomic_context_coverage_table/{benchmark}.log",
    message:
        "Computing genomic context coverage table for {wildcards.benchmark}"
    conda:
        "../envs/bedtools.yaml"
    script:
        "../scripts/compute_coverage_table.py"
