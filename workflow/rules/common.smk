"""
Common helper functions for the Q100 variant benchmark pipeline
"""


# Genomic stratification contexts and their file paths
CONTEXT_PATHS = {
    "TR": "LowComplexity/{ref}_AllTandemRepeats.bed.gz",
    "TR10kb": "LowComplexity/{ref}_AllTandemRepeats_ge101bp_slop5.bed.gz",
    "HP": "LowComplexity/{ref}_AllHomopolymers_ge7bp_imperfectge11bp_slop5.bed.gz",
    "SD": "SegmentalDuplications/{ref}_segdups.bed.gz",
    "SD10kb": "SegmentalDuplications/{ref}_segdups_gt10kb.bed.gz",
    "MAP": "Mappability/{ref}_lowmappabilityall.bed.gz",
}


def get_chromosomes_for_benchmark(benchmark):
    """
    Get chromosome list for a benchmark based on reference genome.
    
    Rules:
    - CHM13: chr1-chr22, chrX, chrY (all benchmarks)
    - GRCh38: chr1-chr22, chrX, chrY (all benchmarks)
    - GRCh37: 1-22, X, Y (all benchmarks)
    - v5q benchmarks include X and Y
    - Historical benchmarks (v421, cmrg, v06, tr) only chr1-22 or 1-22
    
    Args:
        benchmark: Name of the benchmark set
        
    Returns:
        List of chromosomes appropriate for the benchmark
    """
    # Determine reference genome from benchmark name

    # Get base chromosome list
    chroms = [str(i) for i in range(1, 23)]  # Chromosomes 1-22
    
    # For historical benchmarks, exclude X and Y
    if benchmark.startswith("v5q_"):
        chroms = chroms + ["X", "Y"]
    
    ## Add chr prefix for CHM13 and GRCh38
    if "grch38" in benchmark or "chm13" in benchmark:
        chroms = ["chr" + c for c in chroms]
    
    return chroms

def get_v5q_stvar_benchmarks():
    """
    Get list of structural variant (stvar) benchmark sets.
    
    Returns:
        List of benchmark names that are structural variant sets
    """
    return [b for b in config["benchmarksets"].keys() if "stvar" in b and b.startswith("v5q_")]


def get_subset_stats_inputs(wildcards):
    """
    Generate list of inputs for rtg_vcfstats rule for all benchmarks and chromosomes.
    
    Returns:
        List of file paths for rtg_vcfstats outputs
    """
    inputs = []
    for benchmark in config["benchmarksets"]:
        chroms = get_chromosomes_for_benchmark(benchmark)
        for chrom in chroms:
            stats_file = f"results/vcfstats/{benchmark}_{chrom}_vcfstats.txt"
            inputs.append(stats_file)
    return inputs


def get_sv_len_inputs(wildcards):
    """
    Generate list of sv_len.tsv files for all stvar benchmarks.
    
    Returns:
        List of file paths for sv_len outputs
    """
    return [f"results/sv_len/{benchmark}_svlen.tsv" for benchmark in get_v5q_stvar_benchmarks()]


def get_reference_for_benchmark(benchmark):
    """
    Map benchmark name to reference genome for stratifications.
    
    Args:
        benchmark: Name of the benchmark set
        
    Returns:
        Reference name compatible with GIAB stratifications URLs
        (GRCh37, GRCh38, or CHM13 - note no 'v2.0' suffix)
    """
    if "chm13" in benchmark.lower():
        return "CHM13"
    elif "grch38" in benchmark.lower():
        return "GRCh38"
    elif "grch37" in benchmark.lower():
        return "GRCh37"
    else:
        raise ValueError(f"Unknown reference for benchmark: {benchmark}")


def get_context_count_inputs(wildcards):
    """
    Generate list of stratification context count files for all benchmarks.
    
    Returns:
        List of file paths for context count outputs (72 files total)
    """
    contexts = ["TR", "TR10kb", "HP", "SD", "SD10kb", "MAP"]
    inputs = []
    
    for benchmark in config["benchmarksets"]:
        for context in contexts:
            inputs.append(f"results/context_counts/{benchmark}_{context}_counts.tsv")
    
    return inputs
