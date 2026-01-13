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
    threads: 1
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
    threads: 1
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
    Compute intersection metrics for a stratification against benchmark regions.

    Calculates:
    - strat_bp: Total bases in the stratification (after merge)
    - intersect_bp: Bases overlapping between stratification and benchmark regions
    - pct_of_strat: Percent of stratification covered by benchmark (useful for assessing representativeness)
    - pct_of_dip: Percent of benchmark covered by stratification (useful for characterizing difficulty)
    """
    input:
        strat_bed="results/strat_metrics/{benchmark}/stratifications/{strat_name}.bed.gz",
        strat_size="results/strat_metrics/{benchmark}/sizes/{strat_name}_size.txt",
        dip_bed="resources/benchmarksets/{benchmark}_dip.bed",
        dip_size="results/exclusions/{benchmark}/dip_size.txt",  # Reuse from exclusions
    output:
        tsv="results/strat_metrics/{benchmark}/metrics/{strat_name}.tsv",
    log:
        "logs/strat_metrics/{benchmark}/{strat_name}_metrics.log",
    message:
        "Computing metrics for {wildcards.strat_name} in {wildcards.benchmark}"
    threads: 1
    resources:
        mem_mb=4096,
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        echo "Computing metrics for {wildcards.strat_name}" > {log}
        echo "Started at $(date)" >> {log}

        STRAT_SIZE=$(cat {input.strat_size})
        DIP_SIZE=$(cat {input.dip_size})
        echo "Stratification size: $STRAT_SIZE bp" >> {log}
        echo "Dip size: $DIP_SIZE bp" >> {log}

        # Intersect size (overlap between stratification and dip)
        INTERSECT_SIZE=$(bedtools intersect \
            -a <(bedtools sort -i {input.dip_bed}) \
            -b <(gzip -dc {input.strat_bed} | bedtools sort -i -) | \
            awk '{{sum+=$3-$2}} END {{print sum+0}}')
        echo "Intersect size: $INTERSECT_SIZE bp" >> {log}

        # Calculate percentages
        if [ "$STRAT_SIZE" -gt 0 ]; then
            PCT_OF_STRAT=$(awk -v intersect="$INTERSECT_SIZE" -v strat="$STRAT_SIZE" \
                'BEGIN {{printf "%.6f", (intersect / strat) * 100}}')
        else
            PCT_OF_STRAT=0
        fi

        if [ "$DIP_SIZE" -gt 0 ]; then
            PCT_OF_DIP=$(awk -v intersect="$INTERSECT_SIZE" -v dip="$DIP_SIZE" \
                'BEGIN {{printf "%.6f", (intersect / dip) * 100}}')
        else
            PCT_OF_DIP=0
        fi

        echo "Percent of stratification in benchmark: $PCT_OF_STRAT" >> {log}
        echo "Percent of benchmark in stratification: $PCT_OF_DIP" >> {log}

        # Output tab-separated values (header added during aggregation)
        echo -e "{wildcards.strat_name}\\t$STRAT_SIZE\\t$INTERSECT_SIZE\\t$PCT_OF_STRAT\\t$PCT_OF_DIP" > {output.tsv}

        echo "Completed at $(date)" >> {log}
        """


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
        csv=ensure(
            "results/strat_metrics/{benchmark}/stratification_coverage_table.csv",
            non_empty=True,
        ),
    log:
        "logs/strat_metrics/{benchmark}/aggregate.log",
    message:
        "Aggregating stratification metrics for {wildcards.benchmark}"
    threads: 1
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
