"""
Common helper functions for the Q100 variant benchmark pipeline
"""

import os
from urllib.parse import urlparse


def get_var_table_inputs(wildcards):
    """
    Generate list of variant table files for all benchmarks.

    Returns:
        List of file paths for variant table outputs
    """
    return [
        f"results/variant_tables/{benchmark}/variants.tsv"
        for benchmark in config["benchmarksets"]
    ]


# ============================================================================
# Benchmark VCF and BED Helper Functions
# ============================================================================


def get_exclusion_file_path(benchmark, exclusion_name, file_idx):
    """
    Get the standardized local path for an exclusion file.

    Args:
        benchmark: Name of the benchmark set
        exclusion_name: Name of the exclusion category
        file_idx: Index of the file (0-based)

    Returns:
        Local path where the downloaded exclusion file will be stored
    """
    return f"resources/exclusions/{benchmark}/{exclusion_name}_{file_idx}.bed"


def get_exclusion_table_inputs(wildcards):
    """
    Generate list of exclusion intersection tables for all benchmarks
    that have exclusions configured.

    Returns:
        List of file paths for exclusion table outputs
    """
    inputs = []
    for benchmark, conf in config["benchmarksets"].items():
        if "exclusions" in conf and conf["exclusions"]:
            # Include all benchmarks with exclusions configured
            # Files will be downloaded by the download_exclusion rule
            inputs.append(
                f"results/exclusions/{benchmark}/exclusions_intersection_table.csv"
            )
    return inputs


def get_stratification_overlap_inputs(wildcards):
    """
    Generate list of stratification overlap tables for all benchmarks
    that have stratifications configured.

    Returns:
        List of file paths for stratification overlap table outputs
    """
    inputs = []
    for benchmark, conf in config["benchmarksets"].items():
        ref = conf.get("ref")
        if ref and ref in config["references"]:
            ref_config = config["references"][ref]
            if "stratifications" in ref_config and ref_config["stratifications"]:
                inputs.append(
                    f"results/stratification_overlap/{benchmark}/overlap_table.tsv"
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
    """
    Get input file paths for an exclusion.

    Returns standardized paths that match the download_exclusion rule outputs.
    """
    entry = get_exclusion_entry(wildcards.benchmark, wildcards.exclusion)
    return [
        get_exclusion_file_path(wildcards.benchmark, wildcards.exclusion, idx)
        for idx in range(len(entry["files"]))
    ]


def get_exclusion_type(wildcards):
    """Get type (single/pair) for an exclusion."""
    entry = get_exclusion_entry(wildcards.benchmark, wildcards.exclusion)
    return entry["type"]


def get_stratification_beds(wildcards):
    """Get list of stratification BED files with IDs."""
    ref = config["benchmarksets"][wildcards.benchmark].get("ref")
    strats = config["references"][ref].get("stratifications", {})
    beds = []
    for name, strat in strats.items():
        # Use short ID in filename (e.g., CHM13_TR.bed.gz)
        path = f"resources/stratifications/{ref}_{name}.bed.gz"
        beds.append(f"{path}:{name}")
    return beds


def get_region_beds(wildcards):
    """
    Get benchmark region BED and exclusion BEDs.

    Returns paths that match download rule outputs.
    """
    benchmark = wildcards.benchmark
    beds = []

    # Benchmark regions
    bed_path = f"resources/benchmarksets/{benchmark}_benchmark.bed"
    beds.append(f"{bed_path}:BMKREGIONS")

    # Exclusions (only for v5q benchmarks)
    if benchmark.startswith("v5q_"):
        exclusions = get_exclusion_config(benchmark)
        for excl in exclusions:
            name = excl["name"].replace("-", "_").upper()
            for idx, _ in enumerate(excl["files"]):
                path = get_exclusion_file_path(benchmark, excl["name"], idx)
                beds.append(f"{path}:EXCL_{name}")

    return beds


def get_strat_ids(wildcards):
    """Get list of stratification IDs."""
    return list(config.get("stratifications", {}).keys())


def get_region_ids(wildcards):
    """
    Get list of region IDs (benchmark + exclusions).

    All exclusions configured in config are included.
    """
    ids = ["BMKREGIONS"]
    if wildcards.benchmark.startswith("v5q_"):
        exclusions = get_exclusion_config(wildcards.benchmark)
        for excl in exclusions:
            name = excl["name"].replace("-", "_").upper()
            ids.append(f"EXCL_{name}")
    return ids


# ============================================================================
# Reference Download Helper Functions
# ============================================================================


def get_reference_checksum(ref_name):
    """
    Get checksum value for a reference genome.

    Args:
        ref_name: Name of the reference in config["references"]

    Returns:
        Checksum string (MD5 or SHA256)
    """
    ref_config = config["references"].get(ref_name, {})
    if "sha256" in ref_config:
        return ref_config["sha256"]
    if "md5" in ref_config:
        return ref_config["md5"]
    raise ValueError(f"No checksum found for reference: {ref_name}")


def get_reference_checksum_type(ref_name):
    """
    Determine checksum type for a reference genome.

    Args:
        ref_name: Name of the reference in config["references"]

    Returns:
        "md5" or "sha256"
    """
    ref_config = config["references"].get(ref_name, {})
    if "sha256" in ref_config:
        return "sha256"
    if "md5" in ref_config:
        return "md5"
    raise ValueError(f"No checksum found for reference: {ref_name}")


# ============================================================================
# Stratification Download Helper Functions
# ============================================================================


def get_strat_config(ref, strat_name):
    """
    Get stratification configuration for a given reference and stratification name.

    Args:
        ref: Reference genome name (GRCh37, GRCh38, CHM13)
        strat_name: Name of the stratification

    Returns:
        Stratification configuration dictionary
    """
    return (
        config.get("references", {})
        .get(ref, {})
        .get("stratifications", {})
        .get(strat_name, {})
    )


def get_stratification_url(wildcards):
    """
    Get URL for a stratification BED file.

    Handles GIAB stratification URL templating with {ref}@all pattern.

    Args:
        ref: Reference genome name (GRCh37, GRCh38, CHM13)
        strat_name: Name of the stratification

    Returns:
        Full URL for the stratification file
    """
    strat_config = get_strat_config(wildcards.ref, wildcards.strat_name)
    return strat_config.get("url", "")


def get_stratification_sha256(wildcards):
    """
    Get SHA256 checksum for a stratification BED file.

    Args:
        ref: Reference genome name (GRCh37, GRCh38, CHM13)
        strat_name: Name of the stratification

    Returns:
        SHA256 checksum string
    """
    strat_config = get_strat_config(wildcards.ref, wildcards.strat_name)
    return strat_config.get("sha256", "")


def get_stratification_bed_paths(wildcards):
    """
    Get list of stratification BED file paths for a benchmark.

    Args:
        wildcards: Must contain 'benchmark' attribute

    Returns:
        List of paths to stratification BED files
    """
    ref = config["benchmarksets"][wildcards.benchmark].get("ref")
    if not ref or ref not in config["references"]:
        return []

    strats = config["references"][ref].get("stratifications", {})
    paths = []
    for name in strats.keys():
        path = f"resources/stratifications/{ref}_{name}.bed.gz"
        paths.append(path)
    return paths


def get_stratification_names(wildcards):
    """
    Get list of stratification names for a benchmark.

    Args:
        wildcards: Must contain 'benchmark' attribute

    Returns:
        List of stratification names (e.g., ['TR', 'HP', 'SD', 'MAP'])
    """
    ref = config["benchmarksets"][wildcards.benchmark].get("ref")
    if not ref or ref not in config["references"]:
        return []

    strats = config["references"][ref].get("stratifications", {})
    return list(strats.keys())


# ============================================================================
# Exclusion Download Helper Functions
# ============================================================================


def get_exclusion_file_url(benchmark, exclusion_name, file_idx):
    """
    Get URL for an exclusion file by index.

    Args:
        benchmark: Name of the benchmark set
        exclusion_name: Name of the exclusion category
        file_idx: Index of the file in the files array (0-based)

    Returns:
        URL for the exclusion file
    """
    entry = get_exclusion_entry(benchmark, exclusion_name)
    files = entry.get("files", [])
    if file_idx >= len(files):
        raise ValueError(
            f"File index {file_idx} out of range for exclusion {exclusion_name} "
            f"in benchmark {benchmark} (only {len(files)} files)"
        )
    return files[file_idx]["url"]


def get_exclusion_file_checksum(benchmark, exclusion_name, file_idx):
    """
    Get SHA256 checksum for an exclusion file by index.

    Args:
        benchmark: Name of the benchmark set
        exclusion_name: Name of the exclusion category
        file_idx: Index of the file in the files array (0-based)

    Returns:
        SHA256 checksum for the exclusion file
    """
    entry = get_exclusion_entry(benchmark, exclusion_name)
    files = entry.get("files", [])
    if file_idx >= len(files):
        raise ValueError(
            f"File index {file_idx} out of range for exclusion {exclusion_name} "
            f"in benchmark {benchmark} (only {len(files)} files)"
        )
    return files[file_idx]["sha256"]
