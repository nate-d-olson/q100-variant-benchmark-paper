"""
Rules for counting variants by stratification (difficult regions).

This module analyzes variant tables to count variants within each stratification region.
Uses the STRAT_IDS annotation field added during the variant table generation pipeline.

Outputs:
- Per-stratification variant counts by variant type
- Summary tables for cumulative plot generation
"""


rule count_variants_by_stratification:
    """
    Count variants in each stratification region from variant table.

    Reads the variant table TSV and counts variants by:
    - Stratification (from STRAT_IDS field)
    - Variant type (SNP, INDEL, etc.)

    Output format:
    - strat_name: Stratification region name (TR, HP, SD, MAP, etc.)
    - var_type: Variant type
    - count: Number of variants
    """
    input:
        tsv="results/variant_tables/{benchmark}/variants.tsv",
    output:
        csv=ensure(
            "results/var_counts/{benchmark}/variants_by_stratification.csv",
            non_empty=True,
        ),
    log:
        "logs/var_counts/{benchmark}/count_by_strat.log",
    message:
        "Counting variants by stratification for {wildcards.benchmark}"
    resources:
        mem_mb=8192,
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/count_variants_by_strat.py"


rule summarize_variant_counts:
    """
    Create summary table with total variants per stratification.

    Aggregates variant type counts into:
    - strat_name
    - total_variants
    - snp_count (if applicable)
    - indel_count (if applicable)
    - other_count (if applicable)
    """
    input:
        csv="results/var_counts/{benchmark}/variants_by_stratification.csv",
    output:
        summary=ensure(
            "results/var_counts/{benchmark}/stratification_summary.csv",
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
    Combine stratification metrics (from strat_metrics.smk) with variant counts.

    Creates a unified table with:
    - strat_name
    - strat_bp: Total stratification size
    - intersect_bp: Overlap with benchmark
    - pct_of_strat: % of stratification in benchmark
    - pct_of_bench: % of benchmark in stratification
    - total_variants: Total variant count
    - variant_density: Variants per Mb of benchmark-covered stratification

    This table is ready for cumulative plot generation.
    """
    input:
        metrics="results/strat_metrics/{benchmark}/stratification_coverage_table.csv",
        counts="results/var_counts/{benchmark}/stratification_summary.csv",
    output:
        combined=ensure(
            "results/var_counts/{benchmark}/stratification_combined_metrics.csv",
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
