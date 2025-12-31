"""
Rules for calculating stratification overlap with benchmark regions.

Generates tables showing the percent overlap of benchmark regions with 
difficult stratification intervals. For each interval in the stratification 
BED files (TR, HP, SD, MAP, etc.), calculates what percentage of that interval 
is covered by the benchmark regions.

Output table columns:
- stratification: Name of the stratification (TR, HP, SD, etc.)
- chrom: Chromosome
- start: Interval start position
- end: Interval end position
- length: Length of the interval
- overlap_bp: Bases overlapping with benchmark regions
- percent_overlap: Percentage of interval covered by benchmark

This enables visualizations showing:
1. Overall coverage of difficult contexts in the benchmark
2. Distribution of coverage across intervals (cumulative plots)
3. Comparison between different benchmark versions
"""


rule calculate_stratification_overlap:
    """
    Calculate percent overlap of benchmark regions with stratification intervals.
    
    Uses bedtools coverage to find overlaps and calculates the percentage
    of each stratification interval that is covered by the benchmark regions.
    Outputs a table suitable for loading into R for visualization.
    """
    input:
        benchmark_bed="resources/benchmarksets/{benchmark}_benchmark.bed",
        strat_beds=lambda w: get_stratification_bed_paths(w),
    output:
        table=ensure(
            "results/stratification_overlap/{benchmark}/overlap_table.tsv",
            non_empty=True,
        ),
    params:
        strat_names=lambda w: get_stratification_names(w),
    log:
        "logs/stratification_overlap/{benchmark}.log",
    threads: 1
    resources:
        mem_mb=4096,
    conda:
        "../envs/bedtools.yaml"
    script:
        "../scripts/calculate_stratification_overlap.py"
