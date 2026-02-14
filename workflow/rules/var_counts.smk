"""
Rules for counting variants by genomic context (difficult regions).

Reads the Parquet variant table (with correct variant type classification
and size filtering already applied) and counts variants per genomic context.

Outputs:
- Per-genomic context variant counts by variant type (Parquet)
"""


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
        "logs/genomic_context/{benchmark}/count_by_genomic_context.log",
    message:
        "Counting variants by genomic context for {wildcards.benchmark}"
    resources:
        mem_mb=8192,
    conda:
        "../envs/truvari.yaml"
    script:
        "../scripts/count_variants_by_genomic_context.py"
