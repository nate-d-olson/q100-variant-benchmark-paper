"""
Rules for VCF annotation with genomic context and region BED overlaps.

This module prepares BED files and annotates VCFs using bcftools annotate:
1. combine_genomic_context_beds: merge context BEDs into one indexed file
2. combine_region_beds: merge benchmark region + exclusion BEDs
3. generate_annotation_headers: create VCF INFO header lines
4. annotate_vcf_genomic_contexts: add CONTEXT_IDS to VCF
5. annotate_vcf_regions: add REGION_IDS to VCF
"""


rule combine_genomic_context_beds:
    """Combine genomic context BEDs with IDs for bcftools annotation."""
    input:
        beds=lambda w: [b.split(":")[0] for b in get_genomic_context_bed_specs(w)],
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
        bed_specs=lambda w: get_genomic_context_bed_specs(w),
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
        cat > {output.headers} <<'EOF'
##INFO=<ID=CONTEXT_IDS,Number=.,Type=String,Description="Comma-separated list of genomic context region IDs overlapping variant">
##INFO=<ID=REGION_IDS,Number=.,Type=String,Description="Comma-separated list of region IDs (Benchmark, Exclusions) overlapping variant">
EOF
        echo "Generated annotation headers" > {log}
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
