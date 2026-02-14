"""
Rules for generating variant tables using Truvari with column-wise BED annotations.

This module generates Parquet variant tables by reading VCFs with Truvari's
VariantRecord API and annotating variants against genomic context and region
BED files using Python interval overlap queries. Each BED region produces a
boolean column in the output Parquet.

Replaces the previous bcftools annotate approach that used comma-separated
INFO fields (CONTEXT_IDS, REGION_IDS).
"""


rule generate_variant_parquet:
    """
    Generate variant Parquet table from VCF with column-wise BED annotations.

    Reads the pre-processed VCF (split multiallelics or svinfo-annotated),
    loads genomic context and region BED files into interval indices, and
    annotates each variant with boolean columns indicating overlap.

    Output columns:
    - Metadata: bench_version, ref, bench_type
    - Variant: chrom, pos, end, gt, var_type, var_size, szbin
    - Quality: ref_len, alt_len, qual, filter, is_pass
    - Genomic contexts: HP, TR, SD, MAP, etc. (boolean)
    - Regions: in_benchmark (boolean), excl_* (boolean, v5.0q only)
    """
    input:
        vcf=lambda w: (
            "results/run_truvari_anno_svinfo/{benchmark}/svinfo.vcf.gz"
            if w.benchmark.startswith("v06_")
            else "results/split_multiallelics/{benchmark}/split.vcf.gz"
        ),
        vcfidx=lambda w: (
            "results/run_truvari_anno_svinfo/{benchmark}/svinfo.vcf.gz.tbi"
            if w.benchmark.startswith("v06_")
            else "results/split_multiallelics/{benchmark}/split.vcf.gz.tbi"
        ),
        context_beds=lambda w: [
            b.split(":")[0] for b in get_genomic_context_beds(w)
        ],
        region_beds=lambda w: [
            b.split(":")[0] for b in get_region_beds(w)
        ],
    output:
        parquet=ensure(
            "results/variant_tables/{benchmark}/variants.parquet",
            non_empty=True,
        ),
    params:
        bench_type=lambda w: "stvar" if "stvar" in w.benchmark else "smvar",
        context_bed_specs=lambda w: get_genomic_context_beds(w),
        region_bed_specs=lambda w: get_region_beds(w),
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
