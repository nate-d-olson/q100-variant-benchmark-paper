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
        chroms = [f"chr{c}" for c in chroms]

    return chroms


def get_v5q_stvar_benchmarks():
    """
    Get list of structural variant (stvar) benchmark sets.

    Returns:
        List of benchmark names that are structural variant sets
    """
    return [
        b
        for b in config["benchmarksets"].keys()
        if "stvar" in b and b.startswith("v5q_")
    ]


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
    return [
        f"results/sv_len/{benchmark}_svlen.tsv"
        for benchmark in get_v5q_stvar_benchmarks()
    ]


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
        if "cmrg" not in benchmark.lower():
            for context in contexts:
                inputs.append(
                    f"results/context_counts/{benchmark}_{context}_counts.tsv"
                )

    return inputs


# ============================================================================
# Benchmark Region Coverage Helper Functions
# ============================================================================


def get_benchmark_bed(benchmark):
    """
    Get BED file path for a benchmark set.

    Checks config 'bed' key first, then falls back to draft benchmark BEDs
    in data/NIST_HG002_DraftBenchmark_defrabbV0.020-20250117/ for v5q benchmarks.

    Args:
        benchmark: Name of the benchmark set

    Returns:
        Path to benchmark BED file

    Raises:
        ValueError: If no BED file can be found for the benchmark
    """
    # Check config first
    bed_path = config["benchmarksets"].get(benchmark, {}).get("bed")
    if bed_path:
        return bed_path

    # For v5q benchmarks, construct path to draft benchmark BED
    if benchmark.startswith("v5q_"):
        ref = get_reference_for_benchmark(benchmark)
        # Map reference to filename prefix
        ref_prefix = {
            "CHM13": "CHM13v2.0",
            "GRCh38": "GRCh38",
            "GRCh37": "GRCh37",
        }[ref]

        # Determine variant type (smvar or stvar)
        if "smvar" in benchmark:
            vartype = "smvar"
        elif "stvar" in benchmark:
            vartype = "stvar"
        else:
            raise ValueError(f"Unknown variant type in benchmark: {benchmark}")

        return (
            f"data/NIST_HG002_DraftBenchmark_defrabbV0.020-20250117/"
            f"{ref_prefix}_HG2-T2TQ100-V1.1_{vartype}.benchmark.bed"
        )

    raise ValueError(f"No BED file found for benchmark: {benchmark}")


def get_chromosomes_for_coverage(ref, chrom_scope):
    """
    Get chromosome list for coverage calculation by reference and scope.

    Args:
        ref: Reference genome (CHM13, GRCh38, or GRCh37)
        chrom_scope: Either 'autosomes' or 'sexchrom'

    Returns:
        List of chromosome names with appropriate prefix convention
    """
    if chrom_scope == "autosomes":
        chroms = [str(i) for i in range(1, 23)]
    elif chrom_scope == "sexchrom":
        chroms = ["X", "Y"]
    else:
        raise ValueError(f"Unknown chrom_scope: {chrom_scope}")

    # Add chr prefix for CHM13 and GRCh38
    if ref in ["CHM13", "GRCh38"]:
        chroms = [f"chr{c}" for c in chroms]

    return chroms


def get_chromosome_grep_pattern(ref, chrom_scope):
    """
    Get grep pattern for filtering BED files by chromosome.

    Args:
        ref: Reference genome (CHM13, GRCh38, or GRCh37)
        chrom_scope: Either 'autosomes' or 'sexchrom'

    Returns:
        Grep -E compatible regex pattern matching target chromosomes
    """
    chroms = get_chromosomes_for_coverage(ref, chrom_scope)
    # Create pattern that matches chromosome at start of line with tab after
    # e.g., "^chr1\t|^chr2\t|..." or "^1\t|^2\t|..."
    return "|".join([f"^{c}\\t" for c in chroms])


def get_context_coverage_inputs(wildcards):
    """
    Generate list of context coverage files for all benchmarks.

    Coverage is calculated for:
    - All benchmarks (excluding CMRG) × all contexts × autosomes
    - v5q benchmarks only × all contexts × sexchrom

    Returns:
        List of file paths for context coverage outputs
    """
    contexts = ["TR", "TR10kb", "HP", "SD", "SD10kb", "MAP"]
    inputs = []

    for benchmark in config["benchmarksets"]:
        # Skip CMRG benchmarks (consistent with count_context_variants)
        if "cmrg" in benchmark.lower():
            continue

        for context in contexts:
            # All benchmarks get autosome coverage
            inputs.append(
                f"results/context_coverage/{benchmark}_{context}_autosomes_coverage.tsv"
            )

            # Only v5q benchmarks get sex chromosome coverage
            if benchmark.startswith("v5q_"):
                inputs.append(
                    f"results/context_coverage/{benchmark}_{context}_sexchrom_coverage.tsv"
                )

    return inputs
