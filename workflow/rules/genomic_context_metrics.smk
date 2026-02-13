"""
Rules for computing genomic context (difficult region) coverage metrics.

This module computes overlap between genomic context BED files (TR, HP, SD, MAP, etc.)
and benchmark regions (benchmark.bed), generating tables with:
- context_name: Name of the genomic context region
- context_bp: Total bases in the genomic context BED
- intersect_bp: Bases overlapping between genomic context and benchmark regions
- pct_of_context: Percent of genomic context covered by benchmark
- pct_of_bench: Percent of benchmark regions covered by genomic context

Sources genomic context stratifications from references.stratifications in config.
"""


## TODO - remove and replace with input function call for rule that use these
##        symlinked paths
rule materialize_genomic_context:
    """
    Materialize individual genomic context BED file for metrics computation.

    Downloads are handled by download_stratification rule. This rule just
    creates a symlink or copy for consistency with the metrics pipeline.
    """
    input:
        bed=lambda w: f"resources/stratifications/{config['benchmarksets'][w.benchmark]['ref']}_{w.genomic_context}.bed.gz",
    output:
        bed=ensure(
            "results/genomic_context/{benchmark}/coverage/{genomic_context}.bed.gz",
            non_empty=True,
        ),
    params:
        benchmark_ref=lambda w: config["benchmarksets"][w.benchmark]["ref"],
    log:
        "logs/genomic_context/{benchmark}/{genomic_context}_materialize.log",
    message:
        "Materializing genomic context {wildcards.genomic_context} for {wildcards.benchmark}"
    resources:
        mem_mb=1024,
    ## Place holder env to pass lint prior to removal of this rule and replacement with input function call
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        echo "Materializing {wildcards.genomic_context} for {wildcards.benchmark}" > {log}
        echo "Reference: {params.benchmark_ref}" >> {log}
        echo "Input: {input.bed}" >> {log}
        echo "Started at $(date)" >> {log}

        # Create symlink to save space
        ln -sf $(realpath {input.bed}) {output.bed}

        echo "Completed at $(date)" >> {log}
        """


rule compute_genomic_context_size:
    """
    Compute total size of a genomic context region (after merging overlaps).

    This is used as the denominator for pct_of_context calculation.
    Size files are marked as temporary since they're only needed during metrics computation.
    """
    input:
        bed="results/genomic_context/{benchmark}/coverage/{genomic_context}.bed.gz",
    output:
        size=temp(
            "results/genomic_context/{benchmark}/sizes/{genomic_context}_size.txt"
        ),
    log:
        "logs/genomic_context/{benchmark}/{genomic_context}_size.log",
    message:
        "Computing size for {wildcards.genomic_context} in {wildcards.benchmark}"
    resources:
        mem_mb=2048,
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        echo "Computing genomic context size for {wildcards.genomic_context}" > {log}
        echo "Input: {input.bed}" >> {log}
        echo "Started at $(date)" >> {log}

        # Calculate total size (merge overlapping regions first)
        gzip -dc {input.bed} | bedtools sort -i - | bedtools merge -i - | \
            awk '{{sum+=$3-$2}} END {{print sum+0}}' > {output.size} 2>> {log}

        CONTEXT_SIZE=$(cat {output.size})
        echo "Genomic context size: $CONTEXT_SIZE bp" >> {log}
        echo "Completed at $(date)" >> {log}
        """


rule compute_genomic_context_metrics:
    """
    Compute intersection metrics for genomic context against benchmark regions.

    Uses benchmark.bed file for all benchmark sets.
    """
    input:
        strat_bed="results/genomic_context/{benchmark}/coverage/{genomic_context}.bed.gz",
        bench_bed="resources/benchmarksets/{benchmark}_benchmark.bed",
    output:
        tsv="results/genomic_context/{benchmark}/metrics/{genomic_context}.tsv",
    log:
        "logs/genomic_context/{benchmark}/{genomic_context}_metrics.log",
    message:
        "Computing metrics for {wildcards.genomic_context} in {wildcards.benchmark}"
    resources:
        mem_mb=4096,
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/compute_bed_metrics.py"


rule aggregate_genomic_context_metrics:
    """
    Aggregate all genomic context metrics into a single CSV table.

    Output format:
    - context_name: Name of the genomic context
    - context_bp: Total bases in genomic context (merged)
    - intersect_bp: Overlap with benchmark regions
    - pct_of_context: Percent of genomic context covered by benchmark
    - pct_of_bench: Percent of benchmark covered by genomic context
    """
    input:
        lambda wc: expand(
            "results/genomic_context/{benchmark}/metrics/{genomic_context}.tsv",
            benchmark=wc.benchmark,
            genomic_context=get_genomic_context_ids(wc),
        ),
    output:
        csv=ensure(
            "results/genomic_context/{benchmark}/genomic_context_coverage_table.csv",
            non_empty=True,
        ),
    log:
        "logs/genomic_context/{benchmark}/aggregate.log",
    message:
        "Aggregating genomic context metrics for {wildcards.benchmark}"
    resources:
        mem_mb=1024,
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        echo "Aggregating genomic context metrics for {wildcards.benchmark}" > {log}
        echo "Input files: {input}" >> {log}
        echo "Started at $(date)" >> {log}

        # Write header
        echo "context_name,context_bp,intersect_bp,pct_of_context,pct_of_bench" > {output.csv}

        # Concatenate all metric files, converting tabs to commas
        cat {input} | tr '\\t' ',' >> {output.csv}

        echo "Output lines: $(wc -l < {output.csv})" >> {log}
        echo "Completed at $(date)" >> {log}
        """
