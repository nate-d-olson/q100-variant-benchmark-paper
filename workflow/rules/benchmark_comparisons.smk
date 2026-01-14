# workflow/rules/benchmark_comparisons.smk

rule rtg_format:
    input:
        ref="resources/references/{ref_name}.fa.gz",
    output:
        directory("resources/references/{ref_name}.sdf"),
    log:
        "logs/rtg_format/{ref_name}.log",
    conda:
        "../envs/rtg-tools.yaml"
    shell:
        "rtg format -o {output} {input.ref} > {log} 2>&1"


rule compare_smvar:
    input:
        unpack(get_comparison_files),
        sdf=lambda w: f"resources/references/{config['comparisons'][w.comp_id]['ref']}.sdf",
    output:
        dir=directory("results/comparisons/smvar/{comp_id}"),
        tp="results/comparisons/smvar/{comp_id}/tp.vcf.gz",
        fp="results/comparisons/smvar/{comp_id}/fp.vcf.gz",
        fn="results/comparisons/smvar/{comp_id}/fn.vcf.gz",
    params:
        outdir="results/comparisons/smvar/{comp_id}",
    log:
        "logs/compare_smvar/{comp_id}.log",
    conda:
        "../envs/rtg-tools.yaml"
    shell:
        """
        rm -rf {params.outdir}
        rtg vcfeval \
            -b {input.old_vcf} \
            -c {input.new_vcf} \
            -t {input.sdf} \
            --evaluation-regions {input.new_bed} \
            --bed-regions {input.old_bed} \
            -o {params.outdir} \
            > {log} 2>&1
        """


rule run_truvari_bench:
    input:
        ## TODO use union of new and old benchmark regions as bed
        ## TODO add reference fasta to inputs
        unpack(get_comparison_files),
    output:
        dir=directory("results/comparisons/stvar/{comp_id}/bench"),
        bed="results/comparisons/stvar/{comp_id}/regions_union.bed",
        json="results/comparisons/stvar/{comp_id}/bench/summary.json",  # Marker file
    params:
        outdir="results/comparisons/stvar/{comp_id}/bench",
    log:
        "logs/compare_stvar_bench/{comp_id}.log",
    conda:
        "../envs/truvari.yaml"
    shell:
        """
        rm -rf {params.outdir}
        bedtools intersect -a {input.new_bed} -b {input.old_bed} | \
            bedtools sort -i - > {output.bed}
        truvari bench \
            -b {input.old_vcf} \
            -c {input.new_vcf} \
            --pick ac \
            --passonly \
            -r 2000 \
            -C 5000 \
            --reference {input.ref}
            --includebed {input.bed} \
            --output {params.outdir} \
            > {log} 2>&1
        """


rule run_truvari_refine:
    input:
        bench_dir="results/comparisons/stvar/{comp_id}/bench",
        ref=lambda w: f"resources/references/{config['comparisons'][w.comp_id]['ref']}.fa.gz",
    output:
        dir=directory("results/comparisons/stvar/{comp_id}/refine"),
        tp="results/comparisons/stvar/{comp_id}/refine/refine/tp-call.vcf.gz",
        fp="results/comparisons/stvar/{comp_id}/refine/refine/fp.vcf.gz",
        fn="results/comparisons/stvar/{comp_id}/refine/refine/fn.vcf.gz",
    params:
        outdir="results/comparisons/stvar/{comp_id}/refine",
    log:
        "logs/compare_stvar_refine/{comp_id}.log",
    conda:
        "../envs/truvari.yaml"
    shell:
        """
        rm -rf {params.outdir}
        cp -r {input.bench_dir} {params.outdir}
        
        truvari refine \
            --reference {input.ref} \
            --recount \
            --use-region-coords \
            --use-original-vcfs \
            --align mafft \
            {input.bench_dir} \
            --output {params.outdir} \
            > {log} 2>&1
        """




rule stratify_comparison:
    input:
        unpack(get_strat_inputs),
    output:
        variants_csv="results/stats/{comp_id}_variants.csv",
        regions_csv="results/stats/{comp_id}_regions.csv",
    log:
        "logs/stratify/{comp_id}.log",
    conda:
        "../envs/analysis.yaml"
    script:
        "../scripts/stratify_comparison.py"
