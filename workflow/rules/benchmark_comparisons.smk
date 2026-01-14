# workflow/rules/benchmark_comparisons.smk


wildcard_constraints:
    comp_id="[^/]+",


def get_comparison_files(wildcards):
    comp = config["comparisons"][wildcards.comp_id]
    return {
        "new_vcf": f"resources/benchmarksets/{comp['new_benchmark']}_benchmark.vcf.gz",
        "new_bed": f"resources/benchmarksets/{comp['new_benchmark']}_benchmark.bed",
        "old_vcf": f"resources/benchmarksets/{comp['old_benchmark']}_benchmark.vcf.gz",
        "old_bed": f"resources/benchmarksets/{comp['old_benchmark']}_benchmark.bed",
        "ref": f"resources/references/{comp['ref']}.fa.gz",
    }


def get_stratifications_for_comp(wildcards):
    comp = config["comparisons"][wildcards.comp_id]
    ref = comp["ref"]
    strats = config["references"][ref].get("stratifications", {})
    return [f"resources/stratifications/{ref}_{s}.bed.gz" for s in strats]


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
        unpack(get_comparison_files)
    output:
        dir=directory("results/comparisons/stvar/{comp_id}/bench"),
        json="results/comparisons/stvar/{comp_id}/bench/summary.json" # Marker file
    params:
        outdir="results/comparisons/stvar/{comp_id}/bench"
    log:
        "logs/compare_stvar_bench/{comp_id}.log",
    conda:
        "../envs/truvari.yaml"
    shell:
        """
        rm -rf {params.outdir}
        
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
            {params.outdir} \
            > {log} 2>&1
        """


def get_strat_inputs(wildcards):
    comp = config["comparisons"][wildcards.comp_id]
    ctype = comp["type"]
    if ctype == "smvar":
        base = f"results/comparisons/smvar/{wildcards.comp_id}"
        return {
            "tp": f"{base}/tp.vcf.gz",
            "fp": f"{base}/fp.vcf.gz",
            "fn": f"{base}/fn.vcf.gz",
            "new_bed": f"resources/benchmarksets/{comp['new_benchmark']}_benchmark.bed",
            "old_bed": f"resources/benchmarksets/{comp['old_benchmark']}_benchmark.bed",
            "strat_beds": get_stratifications_for_comp(wildcards),
        }
    else:
        # Fallback to bench results if refine is problematic or desired
        # The user asked for refine, but bench results are in results/comparisons/stvar/{comp_id}/bench
        base = f"results/comparisons/stvar/{wildcards.comp_id}/bench"
        return {
            "tp": f"{base}/tp-comp.vcf.gz",
            "fp": f"{base}/fp.vcf.gz",
            "fn": f"{base}/fn.vcf.gz",
            "new_bed": f"resources/benchmarksets/{comp['new_benchmark']}_benchmark.bed",
            "old_bed": f"resources/benchmarksets/{comp['old_benchmark']}_benchmark.bed",
            "strat_beds": get_stratifications_for_comp(wildcards),
        }


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
