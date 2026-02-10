# AI Disclosure: This file was modified with assistance from Claude (Anthropic)
# for code generation and pipeline design.
"""
Rules for exclusion analysis of v5.0q benchmark sets.

Generates three types of output:

1. Per-exclusion impact (Q1): BED overlap metrics + variant counts
   Output: results/exclusions/{benchmark}/exclusion_impact.csv

2. Exclusion interactions (Q2): upset-style decomposition showing which
   combinations of exclusions overlap bases and variants
   Output: results/exclusions/{benchmark}/exclusion_interactions.csv

3. Cross-version comparison (Q3): why old benchmark regions/variants
   are absent in v5.0q (not in dip.bed vs excluded vs other)
   Output: results/exclusions/{comp_id}/old_only_*.csv

Helper functions are defined in rules/common.smk.
"""


rule materialize_exclusion:
    """
    Materialize exclusion BED file from source.
    
    For single files: copies the source BED.
    For start/end pairs: concatenates, sorts, and merges into a single BED.
    
    Input checksums are validated in the shell command.
    """
    input:
        files=get_exclusion_inputs,
    output:
        bed=ensure(
            "results/exclusions/{benchmark}/{exclusion}.bed",
            non_empty=True,
        ),
    params:
        exclusion_type=get_exclusion_type,
    log:
        "logs/exclusions/{benchmark}/{exclusion}_materialize.log",
    message:
        "Materializing exclusion {wildcards.exclusion} for {wildcards.benchmark}"
    resources:
        mem_mb=2048,
    threads: 1
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        echo "Materializing {wildcards.exclusion} ({params.exclusion_type})" > {log}
        echo "Input files: {input.files}" >> {log}
        echo "Started at $(date)" >> {log}
        
        # Materialize the BED file (always sort for bedtools compatibility)
        if [ "{params.exclusion_type}" == "single" ]; then
            bedtools sort -i {input.files} > {output.bed} 2>> {log}
        else
            # Pair: concatenate, sort, and merge
            cat {input.files} | bedtools sort -i - | bedtools merge -i - > {output.bed} 2>> {log}
        fi
        
        echo "Output lines: $(wc -l < {output.bed})" >> {log}
        echo "Completed at $(date)" >> {log}
        """


rule compute_dip_size:
    """
    Compute total size of benchmark regions (dip.bed) for a benchmark set.
    
    This is computed once per benchmark set and reused for all exclusion
    metrics calculations.
    """
    input:
        bed="resources/benchmarksets/{benchmark}_dip.bed",
    output:
        size="results/exclusions/{benchmark}/dip_size.txt",
    log:
        "logs/exclusions/{benchmark}/dip_size.log",
    message:
        "Computing dip.bed size for {wildcards.benchmark}"
    resources:
        mem_mb=2048,
    threads: 1
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        echo "Computing dip.bed size for {wildcards.benchmark}" > {log}
        echo "Input: {input.bed}" >> {log}
        echo "Started at $(date)" >> {log}
        
        # Calculate total size (sort then merge overlapping regions)
        bedtools sort -i {input.bed} | bedtools merge -i - | awk '{{sum+=$3-$2}} END {{print sum+0}}' > {output.size} 2>> {log}
        
        echo "Dip size: $(cat {output.size}) bp" >> {log}
        echo "Completed at $(date)" >> {log}
        """


rule compute_exclusion_metrics:
    """
    Compute intersection metrics for exclusion against benchmark regions.
    """
    input:
        exclusion="results/exclusions/{benchmark}/{exclusion}.bed",
        dip_bed="resources/benchmarksets/{benchmark}_dip.bed",
    output:
        tsv="results/exclusions/{benchmark}/coverage/{exclusion}.tsv",
    log:
        "logs/exclusions/{benchmark}/coverage_{exclusion}.log",
    message:
        "Computing metrics for {wildcards.exclusion} in {wildcards.benchmark}"
    resources:
        mem_mb=4096,
    threads: 1
    conda:
        "../envs/bedtools.yaml"
    script:
        "../scripts/compute_bed_metrics.py"


rule compute_exclusion_impact:
    """
    Compute per-exclusion impact: BED metrics + variant counts.

    Combines per-exclusion BED overlap metrics with variant counts from the
    variant table. Applies size filtering: smvar <50bp, stvar >=50bp.
    """
    input:
        unpack(get_exclusion_impact_inputs),
    output:
        csv=ensure(
            "results/exclusions/{benchmark}/exclusion_impact.csv",
            non_empty=True,
        ),
    params:
        excl_name_mapping=lambda wc: list(
            get_exclusion_name_mapping(wc.benchmark).items()
        ),
    log:
        "logs/exclusions/{benchmark}/impact.log",
    message:
        "Computing exclusion impact for {wildcards.benchmark}"
    resources:
        mem_mb=4096,
    threads: 1
    conda:
        "../envs/bedtools.yaml"
    script:
        "../scripts/count_exclusion_variants.py"


rule compute_exclusion_interactions:
    """
    Compute exclusion interaction (upset-style) decomposition.

    For each unique combination of overlapping exclusions, computes bases
    and variant counts. Uses bedtools multiinter for bases and variant table
    REGION_IDS for variants. Size filtering: smvar <50bp, stvar >=50bp.
    """
    input:
        unpack(get_exclusion_interaction_inputs),
    output:
        csv=ensure(
            "results/exclusions/{benchmark}/exclusion_interactions.csv",
            non_empty=True,
        ),
    params:
        excl_name_mapping=lambda wc: list(
            get_exclusion_name_mapping(wc.benchmark).items()
        ),
    log:
        "logs/exclusions/{benchmark}/interactions.log",
    message:
        "Computing exclusion interactions for {wildcards.benchmark}"
    resources:
        mem_mb=4096,
    threads: 1
    conda:
        "../envs/bedtools.yaml"
    script:
        "../scripts/compute_exclusion_interactions.py"


rule annotate_old_benchmark_status:
    """
    Annotate old benchmark variants/regions with v5.0q exclusion status.

    For variants in the old benchmark that are NOT in v5.0q benchmark regions,
    determines whether they are: not in dip.bed, excluded by an exclusion, or
    in dip.bed but not excluded.
    """
    input:
        unpack(get_old_benchmark_analysis_inputs),
    output:
        variants_tsv=ensure(
            "results/exclusions/{comp_id}/old_only_variants.tsv", non_empty=True
        ),
        summary_csv=ensure(
            "results/exclusions/{comp_id}/old_only_summary.csv", non_empty=True
        ),
        regions_csv=ensure(
            "results/exclusions/{comp_id}/old_only_regions.csv", non_empty=True
        ),
    params:
        comp_type=lambda wc: config["comparisons"][wc.comp_id]["type"],
        excl_names=lambda wc: [
            e["name"]
            for e in get_exclusion_config(
                config["comparisons"][wc.comp_id]["new_benchmark"]
            )
        ],
    log:
        "logs/exclusions/{comp_id}/old_benchmark_status.log",
    message:
        "Annotating old benchmark status for {wildcards.comp_id}"
    resources:
        mem_mb=8192,
    threads: 2
    conda:
        "../envs/bedtools.yaml"
    script:
        "../scripts/annotate_old_benchmark_status.py"
