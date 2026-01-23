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
    resources:
        mem_mb=8192,
        runtime=120,  # 2 hours
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
        unpack(get_comparison_files),
    output:
        multiext(
            "results/comparisons/stvar/{comp_id}/",
            "fn.vcf.gz",
            "fn.vcf.gz.tbi",
            "fp.vcf.gz",
            "fp.vcf.gz.tbi",
            "tp-base.vcf.gz",
            "tp-base.vcf.gz.tbi",
            "tp-comp.vcf.gz",
            "tp-comp.vcf.gz.tbi",
            "candidate.refine.bed",
            "summary.json",
        ),
        union_bed="results/comparisons/stvar/{comp_id}_union.bed",
    params:
        outdir=directory("results/comparisons/stvar/{comp_id}/"),
    log:
        "logs/compare_stvar_bench/{comp_id}.log",
    resources:
        mem_mb=16384,
        runtime=240,  # 4 hours
    conda:
        "../envs/truvari.yaml"
    shell:
        """
        bedtools intersect -a {input.new_bed} -b {input.old_bed} | \
            bedtools sort -i - > {output.union_bed}
        rm -rf {params.outdir}
        truvari bench \
            -b {input.old_vcf} \
            -c {input.new_vcf} \
            --pick ac \
            --passonly \
            -r 2000 \
            -C 5000 \
            --reference {input.ref} \
            --includebed {output.union_bed} \
            --output {params.outdir} \
            > {log} 2>&1
        """


rule run_truvari_refine:
    input:
        unpack(get_comparison_files),
        union_bed="results/comparisons/stvar/{comp_id}_union.bed",
    output:
        multiext(
            "results/comparisons/stvar/{comp_id}/",
            "refine.comp.vcf.gz",
            "refine.base.vcf.gz",
            "refine.variant_summary.json",
        ),
    params:
        bench_dir="results/comparisons/stvar/{comp_id}",
    log:
        "logs/compare_stvar_refine/{comp_id}.log",
    threads: 16
    conda:
        "../envs/truvari.yaml"
    shell:
        """
        ## Note Truvari wiki suggests using --coord 0 for TR benchmark
        truvari refine \
            --reference {input.ref} \
            --use-original-vcfs \
            --threads {threads} \
            --debug \
            {params.bench_dir} \
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
