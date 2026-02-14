"""
Rules for computing genomic context (difficult region) coverage metrics.

This module computes overlap between genomic context BED files (TR, HP, SD, MAP, etc.)
and benchmark regions (benchmark.bed), generating a single CSV table with:
- context_name: Name of the genomic context region
- context_bp: Total bases in the genomic context BED
- intersect_bp: Bases overlapping between genomic context and benchmark regions
- pct_of_context: Percent of genomic context covered by benchmark
- pct_of_bench: Percent of benchmark regions covered by genomic context

Sources genomic context stratifications from references.stratifications in config.
"""


rule compute_genomic_context_coverage_table:
    """
    Compute coverage metrics for all genomic contexts against benchmark regions.

    Processes all context BED files in a single pass, writing one CSV output.
    Context BED files are resolved directly from resources/stratifications/
    via the get_genomic_context_beds input function.
    """
    input:
        context_beds=get_genomic_context_beds,
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
    resources:
        mem_mb=4096,
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/compute_coverage_table.py"
