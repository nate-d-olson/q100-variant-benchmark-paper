"""
Rules for generating variant tables via bcftools annotate and query.

This module replaces the old rtg-tools vcfstats, context counts, and coverage
calculations with a unified bcftools-based approach that generates comprehensive
TSV tables with variant-level annotations.
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


rule extract_info_fields:
    """Extract INFO field names from VCF header."""
    input:
        vcf="results/annotate_vcf_regions/{benchmark}/fully_annotated.vcf.gz",
    output:
        fields=temp("results/extract_info_fields/{benchmark}/info_fields.txt"),
    params:
        exclude=["CONTEXT_IDS", "REGION_IDS"],
    log:
        "logs/extract_info_fields/{benchmark}.log",
    conda:
        "../envs/python.yaml"
    shell:
        """
        python workflow/scripts/extract_info_fields.py \
            --vcf {input.vcf} \
            --exclude {params.exclude} \
            --output {output.fields} 2> {log}
        """


rule generate_var_table:
    """Generate variant table with all annotations."""
    input:
        vcf="results/annotate_vcf_regions/{benchmark}/fully_annotated.vcf.gz",
        vcfidx="results/annotate_vcf_regions/{benchmark}/fully_annotated.vcf.gz.tbi",
        fields="results/extract_info_fields/{benchmark}/info_fields.txt",
        region_bed="results/combine_region_beds/{benchmark}/region_combined.bed.gz",
        region_tbi="results/combine_region_beds/{benchmark}/region_combined.bed.gz.tbi",
    output:
        tsv=protected(
            ensure("results/variant_tables/{benchmark}/variants.tsv", non_empty=True)
        ),
    params:
        strat_ids=lambda w: get_strat_ids(w),
        region_ids=lambda w: get_region_ids(w),
    log:
        "logs/generate_var_table/{benchmark}.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        echo "Generating variant table for {wildcards.benchmark}" > {log}

        # Build base format string
        BASE_FMT="%CHROM\\t%POS\\t%END\\t[%GT]\\t%VKX\\t%TYPE\\t%STRLEN(REF)\\t%STRLEN(ALT)\\t%INFO/CONTEXT_IDS\\t%INFO/REGION_IDS"

        # Add dynamic INFO fields
        DYN_FIELDS=$(cat {input.fields} | sed 's@^@%INFO/@' | paste -sd'\\t' -)

        if [ -n "$DYN_FIELDS" ]; then
            FULL_FMT="$BASE_FMT\\t$DYN_FIELDS\\n"
        else
            FULL_FMT="$BASE_FMT\\n"
        fi

        # Generate table and expand annotations
        bcftools query -HH -f "$FULL_FMT" \
            --exclude 'GT="0/0" | GT="./."' \
            --regions-file {input.region_bed} \
            -o {output.tsv} {input.vcf} 1>> {log} 2>&1

        echo "Completed at $(date)" >> {log}
        echo "Output lines: $(wc -l < {output.tsv})" >> {log}
        """
