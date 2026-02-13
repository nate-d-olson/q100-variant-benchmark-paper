"""
Rules for generating variant tables via bcftools annotate and Truvari.

This module annotates VCFs with genomic context and region BED overlaps via
bcftools, then generates comprehensive Parquet variant tables using Truvari's
vcf_to_df() with correct variant type classification and size filtering.
"""

# Note: Stratification BEDs are downloaded using storage.http()
# The combine_stratification_beds rule accepts URLs directly


rule combine_genomic_context_beds:
    """Combine genomic context BEDs with IDs for bcftools annotation."""
    input:
        beds=lambda w: [b.split(":")[0] for b in get_genomic_context_beds(w)],
    output:
        bed=temp(
            "results/combine_genomic_context_beds/{benchmark}/context_combined.bed"
        ),
        bed_gz=temp(
            ensure(
                "results/combine_genomic_context_beds/{benchmark}/context_combined.bed.gz",
                non_empty=True,
            )
        ),
        tbi=temp(
            "results/combine_genomic_context_beds/{benchmark}/context_combined.bed.gz.tbi"
        ),
    params:
        bed_specs=lambda w: get_genomic_context_beds(w),
    log:
        "logs/combine_genomic_context_beds/{benchmark}.log",
    conda:
        "../envs/python.yaml"
    shell:
        """
        echo "Combining genomic context BEDs for {wildcards.benchmark}" > {log}

        python workflow/scripts/combine_beds_with_id.py \
            --beds {params.bed_specs} \
            --output {output.bed} 2>> {log}

        bgzip -c {output.bed} > {output.bed_gz}
        tabix -p bed {output.bed_gz}

        echo "Completed at $(date)" >> {log}
        """


rule combine_region_beds:
    """Combine benchmark region and exclusion BEDs for annotation."""
    input:
        beds=lambda w: [b.split(":")[0] for b in get_region_beds(w)],
    output:
        bed=temp("results/combine_region_beds/{benchmark}/region_combined.bed"),
        bed_gz=temp(
            ensure(
                "results/combine_region_beds/{benchmark}/region_combined.bed.gz",
                non_empty=True,
            )
        ),
        tbi=temp("results/combine_region_beds/{benchmark}/region_combined.bed.gz.tbi"),
    params:
        bed_specs=lambda w: get_region_beds(w),
    log:
        "logs/combine_region_beds/{benchmark}.log",
    conda:
        "../envs/python.yaml"
    shell:
        """
        echo "Combining region BEDs for {wildcards.benchmark}" > {log}

        python workflow/scripts/combine_beds_with_id.py \
            --beds {params.bed_specs} \
            --output {output.bed} 2>> {log}

        bgzip -c {output.bed} > {output.bed_gz}
        tabix -p bed {output.bed_gz}

        echo "Completed at $(date)" >> {log}
        """


rule generate_annotation_headers:
    """Generate VCF header lines for annotation fields."""
    output:
        headers=temp(
            "results/generate_annotation_headers/{benchmark}/annotation_headers.txt"
        ),
    log:
        "logs/generate_annotation_headers/{benchmark}.log",
    conda:
        "../envs/python.yaml"
    shell:
        """
        python workflow/scripts/generate_header_lines.py \
            --output {output.headers} 2> {log}
        """


rule annotate_vcf_genomic_contexts:
    """Annotate VCF with genomic context region IDs."""
    input:
        vcf=lambda w: (
            (
                "results/run_truvari_anno_svinfo/{benchmark}/svinfo.vcf.gz"
                if w.benchmark.startswith("v06_")
                else "results/split_multiallelics/{benchmark}/split.vcf.gz"
            ),
        ),
        tbi=lambda w: (
            "results/run_truvari_anno_svinfo/{benchmark}/svinfo.vcf.gz.tbi"
            if w.benchmark.startswith("v06_")
            else "results/split_multiallelics/{benchmark}/split.vcf.gz.tbi"
        ),
        context_bed="results/combine_genomic_context_beds/{benchmark}/context_combined.bed.gz",
        context_tbi="results/combine_genomic_context_beds/{benchmark}/context_combined.bed.gz.tbi",
        headers="results/generate_annotation_headers/{benchmark}/annotation_headers.txt",
    output:
        vcf=temp(
            "results/annotate_vcf_genomic_contexts/{benchmark}/context_annotated.vcf.gz"
        ),
    log:
        "logs/annotate_vcf_genomic_contexts/{benchmark}.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        echo "Annotating {wildcards.benchmark} with genomic contexts" > {log}

        bcftools annotate \
            --annotations {input.context_bed} \
            --columns CHROM,FROM,TO,INFO/CONTEXT_IDS \
            --header-lines {input.headers} \
            --output-type z --output {output.vcf} {input.vcf} 2>> {log}

        echo "Completed at $(date)" >> {log}
        """


rule annotate_vcf_regions:
    """Annotate VCF with benchmark region and exclusion IDs."""
    input:
        vcf="results/annotate_vcf_genomic_contexts/{benchmark}/context_annotated.vcf.gz",
        tbi="results/annotate_vcf_genomic_contexts/{benchmark}/context_annotated.vcf.gz.tbi",
        region_bed="results/combine_region_beds/{benchmark}/region_combined.bed.gz",
        region_tbi="results/combine_region_beds/{benchmark}/region_combined.bed.gz.tbi",
    output:
        vcf=temp(
            ensure(
                "results/annotate_vcf_regions/{benchmark}/fully_annotated.vcf.gz",
                non_empty=True,
            )
        ),
    log:
        "logs/annotate_vcf_regions/{benchmark}.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        echo "Annotating {wildcards.benchmark} with regions" > {log}

        bcftools annotate \
            --annotations {input.region_bed} \
            --columns CHROM,FROM,TO,INFO/REGION_IDS \
            --output-type z --output {output.vcf} {input.vcf} 2>> {log}

        echo "Completed at $(date)" >> {log}
        """


rule generate_variant_parquet:
    """
    Generate variant Parquet table from annotated VCF using Truvari.

    Reads fully annotated VCF (with CONTEXT_IDS, REGION_IDS) and produces
    a Parquet table with:
    - Correct variant type classification (SNV/INDEL for smvar, INS/DEL for stvar)
    - Size filtering (smvar <50bp, stvar >=50bp)
    - Truvari size bins (SZBINS)
    - All INFO fields included automatically

    Replaces: extract_info_fields + generate_var_table rules.
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
