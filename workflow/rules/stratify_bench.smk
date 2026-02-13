"""
Truvari Stratify Integration Rules

This module implements stratification analysis for structural variant benchmarks,
analyzing performance across difficult genomic regions (TR, HP, SD, MAP) and their
complements.

For each stratification, generates four analyses:
1. Overlap mode on regions - variants overlapping stratification regions
2. Within mode on regions - variants completely within stratification regions
3. Overlap mode on inverse - variants overlapping non-stratification regions
4. Within mode on inverse - variants completely within non-stratification regions

This provides comprehensive understanding of variant detection in difficult vs
non-difficult genomic contexts.
"""


rule create_inverse_strat_bed:
    """
    Generate complement BED files (genome regions NOT in stratification).

    Uses bedtools complement to create inverse of stratification regions,
    allowing analysis of variant performance in non-difficult regions.
    """
    input:
        strat_bed="resources/stratifications/{ref}_{strat_name}.bed.gz",
        refidx="resources/references/{ref}.fa.gz.fai",
    output:
        inverse_bed="results/stratifications/inverse/{ref}_{strat_name}_inverse.bed",
    log:
        "logs/create_inverse_strat_bed/{ref}_{strat_name}.log",
    conda:
        "../envs/truvari.yaml"
    shell:
        """
        bedtools sort \
            -i {input.strat_bed} \
            -faidx {input.refidx} \
            | bedtools complement -L \
                -g {input.refidx} \
                -i - \
                > {output.inverse_bed} 2> {log}
        """


rule run_truvari_stratify:
    """
    Execute truvari stratify with various configurations.

    Wildcards:
        comp_id: Comparison identifier (e.g., v5_vs_v0.6_stvar)
        strat_name: Stratification name (TR, HP, SD, etc.)
        region_type: "regions" or "inverse"
        mode: "within" or "complement"
    """
    input:
        unpack(get_stratify_inputs),
        multiext(
            "results/comparisons/stvar/{comp_id}/",
            "refine.comp.vcf.gz",
            "refine.base.vcf.gz",
            "refine.variant_summary.json",
        ),
    output:
        tsv="results/stratified_bench/{comp_id}/{strat_name}_{region_type}_{mode}.tsv",
    params:
        bench_dir=lambda w: f"results/comparisons/stvar/{w.comp_id}/",
        complement_flag=lambda w: "--complement" if w.mode != "within" else "",
    log:
        "logs/run_truvari_stratify/{comp_id}_{strat_name}_{region_type}_{mode}.log",
    conda:
        "../envs/truvari.yaml"
    shell:
        """
        truvari stratify \
            {params.complement_flag} \
            {input.bed} \
            {params.bench_dir} \
            > {output.tsv} 2> {log}
        """


rule aggregate_stratified_results:
    """
    Aggregate all TSVs into summary tables.

    Processes all 24 TSVs per comparison (6 strats × 4 modes) and generates:
    - summary_by_variant.tsv: TP, FP, FN, Precision, Recall, F1 per stratification
    - summary_by_region.tsv: Region counts and statistics per stratification
    """
    input:
        tsvs=get_all_stratified_tsvs,
    output:
        variant_summary=ensure(
            "results/stratified_bench/{comp_id}/summary_by_variant.tsv",
            non_empty=True,
        ),
        region_summary=ensure(
            "results/stratified_bench/{comp_id}/summary_by_region.tsv",
            non_empty=True,
        ),
    log:
        "logs/aggregate_stratified_results/{comp_id}.log",
    conda:
        "../envs/truvari.yaml"
    script:
        "../scripts/aggregate_stratified_bench.py"
