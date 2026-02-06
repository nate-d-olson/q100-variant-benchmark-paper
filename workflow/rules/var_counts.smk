"""
Rules for counting variants by genomic context (difficult regions).

This module analyzes variant tables to count variants within each genomic context region.
Uses the CONTEXT_IDS annotation field added during the variant table generation pipeline.

Outputs:
- Per-genomic context variant counts by variant type
- Summary tables for cumulative plot generation
"""


rule count_variants_by_genomic_context:
    """
    Count variants in each genomic context region from variant table.

    Reads the variant table TSV and counts variants by:
    - Genomic context (from CONTEXT_IDS field)
    - Variant type (SNV, INDEL, etc.)

    Output format:
    - context_name: Genomic context region name (TR, HP, SD, MAP, etc.)
    - var_type: Variant type
    - count: Number of variants
    """
    input:
        tsv="results/variant_tables/{benchmark}/variants.tsv",
    output:
        csv=ensure(
            "results/var_counts/{benchmark}/variants_by_genomic_context.csv",
            non_empty=True,
        ),
    log:
        "logs/var_counts/{benchmark}/count_by_genomic_context.log",
    message:
        "Counting variants by genomic context for {wildcards.benchmark}"
    resources:
        mem_mb=8192,
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/count_variants_by_genomic_context.py"


rule summarize_variant_counts:
    """
    Create summary table with total variants per genomic context.

    Aggregates variant type counts into:
    - context_name
    - total_variants
    - snp_count (if applicable)
    - indel_count (if applicable)
    - other_count (if applicable)
    """
    input:
        csv="results/var_counts/{benchmark}/variants_by_genomic_context.csv",
    output:
        summary=ensure(
            "results/var_counts/{benchmark}/genomic_context_summary.csv",
            non_empty=True,
        ),
    log:
        "logs/var_counts/{benchmark}/summarize.log",
    message:
        "Summarizing variant counts for {wildcards.benchmark}"
    resources:
        mem_mb=2048,
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/summarize_var_counts.py"


rule combine_metrics_and_counts:
    """
    Combine genomic context metrics (from strat_metrics.smk) with variant counts.

    Creates a unified table with:
    - context_name
    - context_bp: Total genomic context size
    - intersect_bp: Overlap with benchmark
    - pct_of_context: % of genomic context in benchmark
    - pct_of_bench: % of benchmark in genomic context
    - total_variants: Total variant count
    - variant_density: Variants per Mb of benchmark-covered genomic context

    This table is ready for cumulative plot generation.
    """
    input:
        metrics="results/genomic_context_metrics/{benchmark}/genomic_context_coverage_table.csv",
        counts="results/var_counts/{benchmark}/genomic_context_summary.csv",
    output:
        combined=ensure(
            "results/var_counts/{benchmark}/genomic_context_combined_metrics.csv",
            non_empty=True,
        ),
    log:
        "logs/var_counts/{benchmark}/combine_metrics.log",
    message:
        "Combining metrics and counts for {wildcards.benchmark}"
    resources:
        mem_mb=2048,
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/combine_metrics_counts.py"
