"""
Rules for computing stratification (difficult region) coverage metrics.

This module computes overlap between stratification BED files (TR, HP, SD, MAP, etc.)
and benchmark regions (dip.bed), generating tables with:
- strat_name: Name of the stratification region
- strat_bp: Total bases in the stratification BED
- intersect_bp: Bases overlapping between stratification and benchmark regions
- pct_of_strat: Percent of stratification covered by benchmark
- pct_of_dip: Percent of benchmark regions covered by stratification

Mirrors the structure of exclusions.smk but for genomic context stratifications.
"""


rule materialize_stratification:
    """
    Materialize individual stratification BED file for metrics computation.

    Downloads are handled by download_stratification rule. This rule just
    creates a symlink or copy for consistency with the metrics pipeline.
    """
    input:
        bed=lambda w: f"resources/stratifications/{config['benchmarksets'][w.benchmark]['ref']}_{w.strat_name}.bed.gz",
    output:
        bed=ensure(
            "results/strat_metrics/{benchmark}/stratifications/{strat_name}.bed.gz",
            non_empty=True,
        ),
    params:
        benchmark_ref=lambda w: config["benchmarksets"][w.benchmark]["ref"],
    log:
        "logs/strat_metrics/{benchmark}/{strat_name}_materialize.log",
    message:
        "Materializing stratification {wildcards.strat_name} for {wildcards.benchmark}"
    resources:
        mem_mb=1024,
    shell:
        """
        echo "Materializing {wildcards.strat_name} for {wildcards.benchmark}" > {log}
        echo "Reference: {params.benchmark_ref}" >> {log}
        echo "Input: {input.bed}" >> {log}
        echo "Started at $(date)" >> {log}

        # Create symlink to save space
        ln -sf $(realpath {input.bed}) {output.bed}

        echo "Completed at $(date)" >> {log}
        """


rule compute_stratification_size:
    """
    Compute total size of a stratification region (after merging overlaps).

    This is used as the denominator for pct_of_strat calculation.
    """
    input:
        bed="results/strat_metrics/{benchmark}/stratifications/{strat_name}.bed.gz",
    output:
        size="results/strat_metrics/{benchmark}/sizes/{strat_name}_size.txt",
    log:
        "logs/strat_metrics/{benchmark}/{strat_name}_size.log",
    message:
        "Computing size for {wildcards.strat_name} in {wildcards.benchmark}"
    resources:
        mem_mb=2048,
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        echo "Computing stratification size for {wildcards.strat_name}" > {log}
        echo "Input: {input.bed}" >> {log}
        echo "Started at $(date)" >> {log}

        # Calculate total size (merge overlapping regions first)
        gzip -dc {input.bed} | bedtools sort -i - | bedtools merge -i - | \
            awk '{{sum+=$3-$2}} END {{print sum+0}}' > {output.size} 2>> {log}

        STRAT_SIZE=$(cat {output.size})
        echo "Stratification size: $STRAT_SIZE bp" >> {log}
        echo "Completed at $(date)" >> {log}
        """


rule compute_stratification_metrics:
    """
    Compute intersection metrics for stratification against benchmark regions.
    """
    input:
        strat_bed="results/strat_metrics/{benchmark}/stratifications/{strat_name}.bed.gz",
        dip_bed="resources/benchmarksets/{benchmark}_dip.bed",
    output:
        tsv="results/strat_metrics/{benchmark}/metrics/{strat_name}.tsv",
    log:
        "logs/strat_metrics/{benchmark}/{strat_name}_metrics.log",
    message:
        "Computing metrics for {wildcards.strat_name} in {wildcards.benchmark}"
    resources:
        mem_mb=4096,
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/compute_bed_metrics.py"


rule aggregate_stratification_metrics:
    """
    Aggregate all stratification metrics into a single CSV table.

    Output format:
    - strat_name: Name of the stratification
    - strat_bp: Total bases in stratification (merged)
    - intersect_bp: Overlap with benchmark regions
    - pct_of_strat: Percent of stratification covered by benchmark
    - pct_of_dip: Percent of benchmark covered by stratification
    """
    input:
        lambda wc: expand(
            "results/strat_metrics/{benchmark}/metrics/{strat_name}.tsv",
            benchmark=wc.benchmark,
            strat_name=get_stratification_ids(wc),
        ),
    output:
        csv=protected(ensure(
            "results/strat_metrics/{benchmark}/stratification_coverage_table.csv",
            non_empty=True,
        )),
    log:
        "logs/strat_metrics/{benchmark}/aggregate.log",
    message:
        "Aggregating stratification metrics for {wildcards.benchmark}"
    resources:
        mem_mb=1024,
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        echo "Aggregating stratification metrics for {wildcards.benchmark}" > {log}
        echo "Input files: {input}" >> {log}
        echo "Started at $(date)" >> {log}

        # Write header
        echo "strat_name,strat_bp,intersect_bp,pct_of_strat,pct_of_dip" > {output.csv}

        # Concatenate all metric files, converting tabs to commas
        cat {input} | tr '\\t' ',' >> {output.csv}

        echo "Output lines: $(wc -l < {output.csv})" >> {log}
        echo "Completed at $(date)" >> {log}
        """
