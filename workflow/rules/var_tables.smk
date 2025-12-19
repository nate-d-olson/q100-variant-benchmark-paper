"""
Rules for generating variant tables via bcftools annotate and query.

This module replaces the old rtg-tools vcfstats, context counts, and coverage
calculations with a unified bcftools-based approach that generates comprehensive
TSV tables with variant-level annotations.
"""

rule download_stratification_bed:
    """
    Download GIAB genome stratification BED files.

    Downloads stratification BED files from GIAB FTP for genomic contexts
    (tandem repeats, homopolymers, segmental duplications, low mappability).
    """
    output:
        bed="resources/stratifications/{ref}_{strat_id}.bed.gz",
    params:
        url=lambda wildcards: config["stratifications"][wildcards.strat_id]["url"].format(ref=wildcards.ref),
    log:
        "logs/downloads/stratifications/{ref}_{strat_id}.log",
    message:
        "Downloading {wildcards.ref} {wildcards.strat_id} stratification BED"
    retries: 3
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        echo "Downloading stratification BED: {wildcards.ref} {wildcards.strat_id}" > {log}
        echo "URL: {params.url}" >> {log}
        echo "Started at $(date)" >> {log}

        wget -O {output.bed} {params.url} 2>> {log}

        if [ ! -s {output.bed} ]; then
            echo "ERROR: Downloaded file is empty" >> {log}
            exit 1
        fi

        echo "Completed at $(date)" >> {log}
        echo "File size: $(du -h {output.bed} | cut -f1)" >> {log}
        """


rule combine_stratification_beds:
    """Combine stratification BEDs with IDs for bcftools annotation."""
    input:
        beds=lambda w: [b.split(":")[0] for b in get_stratification_beds(w)],
    output:
        bed=temp("results/var_tables/{benchmark}/strat_combined.bed"),
        bed_gz=ensure("results/var_tables/{benchmark}/strat_combined.bed.gz", non_empty=True),
        tbi="results/var_tables/{benchmark}/strat_combined.bed.gz.tbi",
    params:
        bed_specs=lambda w: get_stratification_beds(w),
    log:
        "logs/var_tables/{benchmark}/combine_strat_beds.log",
    conda:
        "../envs/python.yaml"
    shell:
        """
        echo "Combining stratification BEDs for {wildcards.benchmark}" > {log}

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
        bed=temp("results/var_tables/{benchmark}/region_combined.bed"),
        bed_gz=ensure("results/var_tables/{benchmark}/region_combined.bed.gz", non_empty=True),
        tbi="results/var_tables/{benchmark}/region_combined.bed.gz.tbi",
    params:
        bed_specs=lambda w: get_region_beds(w),
    log:
        "logs/var_tables/{benchmark}/combine_region_beds.log",
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


rule split_multiallelics:
    """Split multiallelic variants BEFORE annotation."""
    input:
        vcf=lambda w: get_benchmark_vcf(w.benchmark),
        tbi=lambda w: get_benchmark_vcf(w.benchmark) + ".tbi",
        ref=lambda w: f"resources/references/{get_reference_for_benchmark(w.benchmark).lower()}.fa.gz",
        fai=lambda w: f"resources/references/{get_reference_for_benchmark(w.benchmark).lower()}.fa.gz.fai",
    output:
        vcf="results/annotated_vcfs/{benchmark}/split.vcf.gz",
    params:
        old_rec_tag=lambda w: "--old-rec-tag SVLEN" if "stvar" in w.benchmark else "",
    log:
        "logs/var_tables/{benchmark}/split_multiallelic.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        echo "Splitting multiallelics for {wildcards.benchmark}" > {log}

        bcftools norm --multiallelics -any \
            --check-ref w \
            {params.old_rec_tag} \
            --fasta-ref {input.ref} \
            --output-type z --output {output.vcf} {input.vcf} 2>> {log}

        echo "Completed at $(date)" >> {log}
        """

rule run_truvari_anno_svinfo:
    input:
        vcf=lambda w: get_benchmark_vcf(w.benchmark),
        tbi=lambda w: get_benchmark_vcf(w.benchmark) + ".tbi",
    output:
        vcf="results/annotated_vcfs/{benchmark}_svinfo.vcf.gz",
    log:
        "logs/truvari_anno_svinfo/{benchmark}.log",
    conda:
        "../envs/truvari.yaml"
    params:
        minsize=20,
        vcf=lambda w, output: output.vcf[:-3],  # Remove .gz extension
    shell:
        """
        truvari anno svinfo \
            -o {params.vcf} \
            --minsize {params.minsize} \
            {input.vcf} >> {log} 2>&1
        bgzip {params.vcf} >> {log}
        """

rule generate_annotation_headers:
    """Generate VCF header lines for annotation fields."""
    output:
        headers="results/var_tables/{benchmark}/annotation_headers.txt"
    log:
        "logs/var_tables/{benchmark}/generate_headers.log",
    conda:
        "../envs/python.yaml"
    shell:
        """
        python workflow/scripts/generate_header_lines.py \
            --output {output.headers} 2> {log}
        """


rule annotate_vcf_stratifications:
    """Annotate VCF with stratification region IDs."""
    input:
        vcf=lambda w: "results/annotated_vcfs/{benchmark}_svinfo.vcf.gz" if w.benchmark.startswith("v06_") else "results/annotated_vcfs/{benchmark}_svinfo.vcf.gz",
        tbi=lambda w: "results/annotated_vcfs/{benchmark}_svinfo.vcf.gz" + ".tbi" if w.benchmark.startswith("v06_") else "results/annotated_vcfs/{benchmark}_svinfo.vcf.gz.tbi",
        strat_bed="results/var_tables/{benchmark}/strat_combined.bed.gz",
        strat_tbi="results/var_tables/{benchmark}/strat_combined.bed.gz.tbi",
        headers="results/var_tables/{benchmark}/annotation_headers.txt",
    output:
        vcf="results/annotated_vcfs/{benchmark}/strat_annotated.vcf.gz",
    log:
        "logs/var_tables/{benchmark}/annotate_strat.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        echo "Annotating {wildcards.benchmark} with stratifications" > {log}

        bcftools annotate \
            --annotations {input.strat_bed} \
            --columns CHROM,FROM,TO,INFO/STRAT_IDS \
            --header-lines {input.headers} \
            --output-type z --output {output.vcf} {input.vcf} 2>> {log}

        echo "Completed at $(date)" >> {log}
        """


rule annotate_vcf_regions:
    """Annotate VCF with benchmark region and exclusion IDs."""
    input:
        vcf="results/annotated_vcfs/{benchmark}/strat_annotated.vcf.gz",
        tbi="results/annotated_vcfs/{benchmark}/strat_annotated.vcf.gz.tbi",
        region_bed="results/var_tables/{benchmark}/region_combined.bed.gz",
        region_tbi="results/var_tables/{benchmark}/region_combined.bed.gz.tbi",
    output:
        vcf=ensure("results/annotated_vcfs/{benchmark}/fully_annotated.vcf.gz", non_empty=True),
    log:
        "logs/var_tables/{benchmark}/annotate_regions.log",
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
        vcf="results/annotated_vcfs/{benchmark}/fully_annotated.vcf.gz"
    output:
        fields="results/var_tables/{benchmark}/info_fields.txt"
    params:
        exclude=["STRAT_IDS", "REGION_IDS"]
    log:
        "logs/var_tables/{benchmark}/extract_fields.log",
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
        vcf="results/annotated_vcfs/{benchmark}/fully_annotated.vcf.gz",
        fields="results/var_tables/{benchmark}/info_fields.txt",
    output:
        tsv=ensure("results/var_tables/{benchmark}/variants.tsv", non_empty=True),
    params:
        strat_ids=lambda w: get_strat_ids(w),
        region_ids=lambda w: get_region_ids(w),
    log:
        "logs/var_tables/{benchmark}/generate_table.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        echo "Generating variant table for {wildcards.benchmark}" > {log}

        # Build base format string
        BASE_FMT="%CHROM\\t%POS\\t%END\\t[%GT]\\t%VKX\\t%TYPE\\t%STRLEN(REF)\\t%STRLEN(ALT)\\t%INFO/STRAT_IDS\\t%INFO/REGION_IDS"

        # Add dynamic INFO fields
        DYN_FIELDS=$(cat {input.fields} | sed 's@^@%INFO/@' | paste -sd'\\t' -)

        if [ -n "$DYN_FIELDS" ]; then
            FULL_FMT="$BASE_FMT\\t$DYN_FIELDS\\n"
        else
            FULL_FMT="$BASE_FMT\\n"
        fi

        # # Build header
        # BASE_HDR="CHROM\\tPOS\\tEND\\tGT\\tVKX\\tN_ALT\\tTYPE\\tSTRLEN_REF\\tSTRLEN_ALT\\tILEN\\tSTRAT_IDS\\tREGION_IDS"
        # DYN_HDR=$(cat {input.fields} | paste -sd'\\t' -)

        # if [ -n "$DYN_HDR" ]; then
        #     HEADER="$BASE_HDR\\t$DYN_HDR"
        # else
        #     HEADER="$BASE_HDR"
        # fi

        # Generate table and expand annotations
        bcftools query -HH -f "$FULL_FMT" -o {output.tsv} {input.vcf} 1>> {log} 2>&1
        # (echo -e "$HEADER"; bcftools query -f "$FULL_FMT" {input.vcf}) | \
        # python workflow/scripts/expand_annotations.py \
        #     --strat-ids {params.strat_ids} \
        #     --region-ids {params.region_ids} > {output.tsv} 2>> {log}

        echo "Completed at $(date)" >> {log}
        echo "Output lines: $(wc -l < {output.tsv})" >> {log}
        """
