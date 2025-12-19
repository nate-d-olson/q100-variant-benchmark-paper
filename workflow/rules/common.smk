"""
Common helper functions for the Q100 variant benchmark pipeline
"""

import os
from urllib.parse import urlparse


def get_filename_from_url(url):
    """Extract filename from URL."""
    return os.path.basename(urlparse(url).path)


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


def get_var_table_inputs(wildcards):
    """
    Generate list of variant table files for all benchmarks.

    Returns:
        List of file paths for variant table outputs
    """
    return [
        f"results/var_tables/{benchmark}/variants.tsv"
        for benchmark in config["benchmarksets"]
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


# ============================================================================
# Benchmark VCF and BED Helper Functions
# ============================================================================


def get_benchmark_vcf(benchmark):
    """
    Get VCF file path for a benchmark set.

    Args:
        benchmark: Name of the benchmark set

    Returns:
        Path to benchmark VCF file
    """
    vcf_entry = config["benchmarksets"].get(benchmark, {}).get("vcf")
    if isinstance(vcf_entry, dict):
        if url := vcf_entry.get("url"):
            return os.path.join("resources/benchmarksets", get_filename_from_url(url))
        if path := vcf_entry.get("path"):
            return path

    raise ValueError(f"No VCF file found for benchmark: {benchmark}")


def get_benchmark_bed(benchmark):
    """
    Get BED file path for a benchmark set based on url or path in config.

    Args:
        benchmark: Name of the benchmark set

    Returns:
        Path to benchmark BED file

    Raises:
        ValueError: If no BED file can be found for the benchmark
    """
    # Check config first
    bed_entry = config["benchmarksets"].get(benchmark, {}).get("bed")

    if isinstance(bed_entry, str):
        return bed_entry

    if isinstance(bed_entry, dict):
        if url := bed_entry.get("url"):
            return os.path.join("resources/benchmarksets", get_filename_from_url(url))
        if path := bed_entry.get("path"):
            return path

    raise ValueError(f"No BED file found for benchmark: {benchmark}")


def get_all_reference_files(wildcards):
    """
    Get all reference files defined in the config.

    Returns:
        List of file paths for all reference files
    """
    return [
        ref["path"] for ref in config.get("references", {}).values() if "path" in ref
    ]


def get_exclusion_table_inputs(wildcards):
    """
    Generate list of exclusion intersection tables for all benchmarks
    that have exclusions configured and files exist.

    Returns:
        List of file paths for exclusion table outputs
    """
    inputs = []
    for benchmark, conf in config["benchmarksets"].items():
        if "exclusions" in conf and conf["exclusions"]:
            # Only include if at least one exclusion has existing files
            exclusions = conf["exclusions"]
            has_any_files = False
            for excl in exclusions:
                if any(os.path.exists(f['path']) for f in excl.get("files", [])):
                    has_any_files = True
                    break
            if has_any_files:
                inputs.append(
                    f"results/exclusions/{benchmark}/exclusions_intersection_table.csv"
                )
    return inputs


# ============================================================================
# Exclusion Analysis Helper Functions
# ============================================================================


def get_exclusion_config(benchmark):
    """Get exclusion configuration for a benchmark set."""
    return config["benchmarksets"].get(benchmark, {}).get("exclusions", [])


def get_exclusion_items(wildcards):
    """Get list of exclusion names for a benchmark set."""
    exclusions = get_exclusion_config(wildcards.benchmark)
    return [item["name"] for item in exclusions]


def get_exclusion_entry(benchmark, exclusion_name):
    """Get a specific exclusion entry by name."""
    exclusions = get_exclusion_config(benchmark)
    for item in exclusions:
        if item["name"] == exclusion_name:
            return item
    raise ValueError(f"Exclusion {exclusion_name} not found for {benchmark}")


def get_exclusion_inputs(wildcards):
    """Get input file paths for an exclusion."""
    entry = get_exclusion_entry(wildcards.benchmark, wildcards.exclusion)
    return [f["path"] for f in entry["files"]]


def get_exclusion_type(wildcards):
    """Get type (single/pair) for an exclusion."""
    entry = get_exclusion_entry(wildcards.benchmark, wildcards.exclusion)
    return entry["type"]


def get_dip_bed_for_exclusions(wildcards):
    """Get the dip.bed path for a benchmark set."""
    return config["benchmarksets"][wildcards.benchmark]["dip_bed"]["path"]


def get_input_checksums(wildcards):
    """Get checksums for exclusion input files."""
    entry = get_exclusion_entry(wildcards.benchmark, wildcards.exclusion)
    return [f["sha256"] for f in entry["files"]]



def get_stratification_beds(wildcards):
    """Get list of stratification BED files with IDs."""
    ref = get_reference_for_benchmark(wildcards.benchmark)
    strats = config.get("stratifications", {})
    beds = []
    for name, strat in strats.items():
        # Use short ID in filename (e.g., CHM13_TR.bed.gz)
        path = f"resources/stratifications/{ref}_{name}.bed.gz"
        beds.append(f"{path}:{name}")
    return beds


def get_region_beds(wildcards):
    """Get benchmark region BED and exclusion BEDs (only if files exist)."""
    import os
    benchmark = wildcards.benchmark
    beds = []

    # Benchmark regions
    bed_path = get_benchmark_bed(benchmark)
    beds.append(f"{bed_path}:BMKREGIONS")

    # Exclusions (only for v5q benchmarks, and only if files exist)
    if benchmark.startswith("v5q_"):
        exclusions = get_exclusion_config(benchmark)
        for excl in exclusions:
            name = excl["name"].replace("-", "_").upper()
            for f in excl["files"]:
                # Only add if file exists
                if os.path.exists(f['path']):
                    beds.append(f"{f['path']}:EXCL_{name}")

    return beds


def get_strat_ids(wildcards):
    """Get list of stratification IDs."""
    return list(config.get("stratifications", {}).keys())


def get_region_ids(wildcards):
    """Get list of region IDs (benchmark + exclusions that exist)."""
    import os
    ids = ["BMKREGIONS"]
    if wildcards.benchmark.startswith("v5q_"):
        exclusions = get_exclusion_config(wildcards.benchmark)
        for excl in exclusions:
            name = excl["name"].replace("-", "_").upper()
            # Only add if at least one file exists for this exclusion
            if any(os.path.exists(f['path']) for f in excl["files"]):
                ids.append(f"EXCL_{name}")
    return ids

# Function to determine checksum type and return appropriate ensure() argument
def get_checksum_arg(checksum):
    if not checksum:
        return {}
    # MD5 is 32 hex chars, SHA256 is 64 hex chars
    if len(checksum) == 32:
        return {"md5": checksum}
    elif len(checksum) == 64:
        return {"sha256": checksum}
    else:
        # Default to md5 if length doesn't match standard sizes, or let Snakemake handle it
        return {"md5": checksum}
