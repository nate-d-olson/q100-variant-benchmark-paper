"""
Rules for genomic context (difficult region) analysis.

This module computes coverage and variant metrics for genomic context regions
(TR, HP, SD, MAP, etc.) across benchmark sets. The analysis flow is:

1. genomic_context_coverage: bedtools coverage of each context against benchmark
2. compute_genomic_context_coverage_table: summarize coverage into a CSV table
3. generate_variant_parquet: convert annotated VCF to Parquet variant table
4. count_variants_by_genomic_context: count variants per context from Parquet

Outputs:
- results/genomic_context/{benchmark}/coverage/{context}_cov.bed
- results/genomic_context/{benchmark}/genomic_context_coverage_table.csv
- results/variant_tables/{benchmark}/variants.parquet
- results/genomic_context/{benchmark}/variants_by_genomic_context.parquet
"""


rule genomic_context_coverage:
    """Compute base-level coverage for genomic context regions against benchmark."""
    input:
        strat_bed="resources/stratifications/{ref}_{genomic_context}.bed.gz",
        bench_bed="resources/benchmarksets/{bench_version}_{ref}_{bench_type}_benchmark.bed",
    output:
        cov_bed="results/genomic_context/{bench_version}_{ref}_{bench_type}/coverage/{genomic_context}_cov.bed",
    log:
        "logs/genomic_context_coverage/{bench_version}_{ref}_{bench_type}_{genomic_context}.log",
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        bedtools coverage -a {input.strat_bed} -b {input.bench_bed} 1> {output.cov_bed} 2> {log}
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


rule generate_variant_parquet:
    """
    Generate variant Parquet table from annotated VCF using Truvari.

    Reads fully annotated VCF (with CONTEXT_IDS, REGION_IDS) and produces
    a Parquet table with:
    - Correct variant type classification (SNV/INDEL for smvar, INS/DEL for stvar)
    - Size filtering (smvar <50bp, stvar >=50bp)
    - Truvari size bins (SZBINS)
    - All INFO fields included automatically
    """
    input:
        vcf="results/annotate_vcf_regions/{benchmark}/fully_annotated.vcf.gz",
        vcfidx="results/annotate_vcf_regions/{benchmark}/fully_annotated.vcf.gz.tbi",
    output:
        parquet=ensure(
            "results/variant_tables/{benchmark}/variants.parquet",
            non_empty=True,
        ),
    params:
        bench_type=lambda w: "stvar" if "stvar" in w.benchmark else "smvar",
    log:
        "logs/generate_variant_parquet/{benchmark}.log",
    message:
        "Generating variant Parquet for {wildcards.benchmark}"
    resources:
        mem_mb=16384,
    conda:
        "../envs/truvari.yaml"
    script:
        "../scripts/generate_variant_parquet.py"


rule count_variants_by_genomic_context:
    """
    Count variants in each genomic context region from Parquet variant table.

    Variant types are already correctly classified:
    - smvar: SNV, INDEL (all non-SNV <50bp)
    - stvar: INS, DEL (>=50bp)

    Size filtering is already applied during Parquet generation.

    Output columns:
    - context_name: Genomic context region name (TR, HP, SD, MAP, etc.)
    - var_type: Classified variant type
    - count: Number of variants
    - szbin: Truvari size bin (when counting by size)
    """
    input:
        parquet="results/variant_tables/{benchmark}/variants.parquet",
    output:
        parquet=ensure(
            "results/genomic_context/{benchmark}/variants_by_genomic_context.parquet",
            non_empty=True,
        ),
    log:
        "logs/count_variants_by_genomic_context/{benchmark}.log",
    message:
        "Counting variants by genomic context for {wildcards.benchmark}"
    resources:
        mem_mb=8192,
    conda:
        "../envs/truvari.yaml"
    script:
        "../scripts/count_variants_by_genomic_context.py"
