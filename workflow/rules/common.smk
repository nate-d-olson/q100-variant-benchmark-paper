# AI Disclosure: This file was modified with assistance from Claude (Anthropic)
# to add exclusion analysis helper functions.
"""
Common helper functions for the Q100 variant benchmark pipeline
"""

from typing import List, Dict, Any, Tuple

# ============================================================================
# Configuration & Constraints
# ============================================================================

# Collect all unique stratification names from all references
_all_strat_names = set()
for ref_config in config.get("references", {}).values():
    _all_strat_names.update(ref_config.get("stratifications", {}).keys())

# Genomic context names (same sources as stratifications, used in updated analysis pipeline)
_all_genomic_context_names = _all_strat_names.copy()


wildcard_constraints:
    comp_id="[^/]+",
    strat_name="|".join(sorted(_all_strat_names)),
    genomic_context="|".join(sorted(_all_genomic_context_names)),
    region_type="regions|inverse",
    mode="complement|within",


# ============================================================================
# Benchmark Helpers
# ============================================================================


def get_bench_ids(wildcards) -> List[str]:
    """Get list of benchmark IDs."""
    return list(config.get("benchmarksets", {}).keys())


def get_var_table_inputs(wildcards) -> List[str]:
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
# Reference Helpers
# ============================================================================


def get_ref_ids(wildcards) -> List[str]:
    """Get list of reference IDs."""
    return list(config.get("references", {}).keys())


def get_chromosomes(wildcards) -> str:
    """
    Get list of main chromosomes for the reference genome.

    Args:
        wildcards: Snakemake wildcards with `ref` attribute

    Returns:
        Space-separated string of chromosomes
    """
    chroms = [str(i) for i in range(1, 23)] + ["X", "Y"]
    if wildcards.ref == "GRCh37":
        return " ".join(chroms)
    return " ".join([f"chr{c}" for c in chroms])


def get_reference_checksum(ref_name: str) -> str:
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


# ============================================================================
# Variant Table Helpers
# ============================================================================


def get_genomic_context_beds(wildcards) -> List[str]:
    """Get list of genomic context BED files with IDs."""
    ref = config["benchmarksets"][wildcards.benchmark].get("ref")
    contexts = config["references"][ref].get("stratifications", {})
    beds = []
    for name, _ in contexts.items():
        # Use short ID in filename (e.g., CHM13_TR.bed.gz)
        path = f"resources/stratifications/{ref}_{name}.bed.gz"
        beds.append(f"{path}:{name}")
    return beds


def get_region_beds(wildcards) -> List[str]:
    """
    Get benchmark region BED and exclusion BEDs.

    Returns paths that match download rule outputs.
    """
    benchmark = wildcards.benchmark
    beds = []

    # Benchmark regions
    bed_path = f"resources/benchmarksets/{benchmark}_benchmark.bed"
    beds.append(f"{bed_path}:BMKREGIONS")

    # Exclusions
    exclusions = get_exclusion_config(benchmark)
    for excl in exclusions:
        name = _format_exclusion_name(excl["name"])
        for idx, _ in enumerate(excl["files"]):
            path = get_exclusion_file_path(benchmark, excl["name"], idx)
            beds.append(f"{path}:EXCL_{name}")

    return beds


def get_region_ids(wildcards) -> List[str]:
    """
    Get list of region IDs (benchmark + exclusions).

    All exclusions configured in config are included.
    """
    ids = ["BMKREGIONS"]
    exclusions = get_exclusion_config(wildcards.benchmark)
    for excl in exclusions:
        name = _format_exclusion_name(excl["name"])
        ids.append(f"EXCL_{name}")
    return ids


# ============================================================================
# Comparison Helpers
# ============================================================================


def get_comparison_files(wildcards):
    comp = config["comparisons"][wildcards.comp_id]
    return {
        "new_vcf": f"resources/benchmarksets/{comp['new_benchmark']}_benchmark.vcf.gz",
        "new_vcfidx": f"resources/benchmarksets/{comp['new_benchmark']}_benchmark.vcf.gz.tbi",
        "new_bed": f"resources/benchmarksets/{comp['new_benchmark']}_benchmark.bed",
        "old_vcf": f"resources/benchmarksets/{comp['old_benchmark']}_benchmark.vcf.gz",
        "old_vcfidx": f"resources/benchmarksets/{comp['old_benchmark']}_benchmark.vcf.gz.tbi",
        "old_bed": f"resources/benchmarksets/{comp['old_benchmark']}_benchmark.bed",
        "ref": f"resources/references/{comp['ref']}.fa.gz",
    }


def get_stratifications_for_comp(wildcards):
    comp = config["comparisons"][wildcards.comp_id]
    ref = comp["ref"]
    strats = config["references"][ref].get("stratifications", {})
    return [f"resources/stratifications/{ref}_{s}.bed.gz" for s in strats]


def get_strat_inputs(wildcards):
    comp = config["comparisons"][wildcards.comp_id]
    ctype = comp["type"]

    inputs = {
        "new_bed": f"resources/benchmarksets/{comp['new_benchmark']}_benchmark.bed",
        "old_bed": f"resources/benchmarksets/{comp['old_benchmark']}_benchmark.bed",
        "strat_beds": get_stratifications_for_comp(wildcards),
    }

    if ctype == "smvar":
        base = f"results/comparisons/smvar/{wildcards.comp_id}"
        inputs.update(
            {
                "tp": f"{base}/tp.vcf.gz",
                "fp": f"{base}/fp.vcf.gz",
                "fn": f"{base}/fn.vcf.gz",
            }
        )
    else:
        # stvar (Truvari) - using refined outputs
        base = f"results/comparisons/stvar/{wildcards.comp_id}"
        return {
            "base": f"{base}/refine.base.vcf.gz",
            "comp": f"{base}/refine.comp.vcf.gz",
            "new_bed": f"resources/benchmarksets/{comp['new_benchmark']}_benchmark.bed",
            "old_bed": f"resources/benchmarksets/{comp['old_benchmark']}_benchmark.bed",
            "strat_beds": get_stratifications_for_comp(wildcards),
        }

    return inputs


def get_stratifications_for_ref(ref: str) -> List[str]:
    """Get list of stratification names for a reference."""
    return list(config["references"][ref].get("stratifications", {}).keys())


def get_stratification_ids(wildcards) -> List[str]:
    """
    Get list of stratification IDs for a benchmark.

    Args:
        wildcards: Snakemake wildcards with benchmark attribute

    Returns:
        List of stratification names for the benchmark's reference
    """
    benchmark = wildcards.benchmark
    ref = config["benchmarksets"][benchmark]["ref"]
    return get_stratifications_for_ref(ref)


def get_genomic_context_ids(wildcards) -> List[str]:
    """
    Get list of genomic context IDs for a benchmark.

    These are the same sources as stratifications but used in the updated analysis pipeline.

    Args:
        wildcards: Snakemake wildcards with benchmark attribute

    Returns:
        List of genomic context names for the benchmark's reference
    """
    benchmark = wildcards.benchmark
    ref = config["benchmarksets"][benchmark]["ref"]
    return get_stratifications_for_ref(ref)


def get_stratify_inputs(wildcards):
    """Get inputs for truvari stratify based on region type."""
    comp = config["comparisons"][wildcards.comp_id]
    ref = comp["ref"]

    inputs = {
        "bench_dir": f"results/comparisons/stvar/{comp}/",
        "summary_json": f"results/comparisons/stvar/{comp}/summary.json",
    }

    # Select bed file based on region_type
    if wildcards.region_type == "regions":
        inputs["bed"] = f"resources/stratifications/{ref}_{wildcards.strat_name}.bed.gz"
    else:  # inverse
        inputs["bed"] = (
            f"results/stratifications/inverse/{ref}_{wildcards.strat_name}_inverse.bed"
        )

    return inputs


def get_all_stratified_tsvs(wildcards):
    """Get all stratified TSV outputs for a comparison."""
    comp = config["comparisons"][wildcards.comp_id]
    ref = comp["ref"]
    strats = config["references"][ref].get("stratifications", {}).keys()

    tsvs = []
    for strat in strats:
        for region_type in ["regions", "inverse"]:
            for mode in ["overlap", "within"]:
                tsvs.append(
                    f"results/stratified_bench/{wildcards.comp_id}/"
                    f"{strat}_{region_type}_{mode}.tsv"
                )

    return tsvs


# ============================================================================
# Exclusion Helpers
# ============================================================================


def get_exclusion_config(benchmark: str) -> List[Dict[str, Any]]:
    """Get exclusion configuration for a benchmark set."""
    return config["benchmarksets"].get(benchmark, {}).get("exclusions", [])


def get_exclusion_entry(benchmark: str, exclusion_name: str) -> Dict[str, Any]:
    """Get a specific exclusion entry by name."""
    exclusions = get_exclusion_config(benchmark)
    for item in exclusions:
        if item["name"] == exclusion_name:
            return item
    raise ValueError(f"Exclusion {exclusion_name} not found for {benchmark}")


def _format_exclusion_name(name: str) -> str:
    """Format exclusion name for ID generation."""
    return name.replace("-", "_").upper()


def get_exclusion_file_path(benchmark: str, exclusion_name: str, file_idx: int) -> str:
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


def get_exclusion_table_inputs(wildcards) -> List[str]:
    """
    Generate list of exclusion intersection tables for all benchmarks
    that have exclusions configured.

    Returns:
        List of file paths for exclusion table outputs
    """
    inputs = []
    for benchmark, conf in config["benchmarksets"].items():
        if conf.get("exclusions"):
            # Include all benchmarks with exclusions configured
            inputs.append(
                f"results/exclusions/{benchmark}/exclusions_intersection_table.csv"
            )
    return inputs


def get_exclusion_items(wildcards) -> List[str]:
    """Get list of exclusion names for a benchmark set."""
    exclusions = get_exclusion_config(wildcards.benchmark)
    return [item["name"] for item in exclusions]


def get_exclusion_inputs(wildcards) -> List[str]:
    """
    Get input file paths for an exclusion.

    Returns standardized paths that match the download_exclusion rule outputs.
    """
    entry = get_exclusion_entry(wildcards.benchmark, wildcards.exclusion)
    return [
        get_exclusion_file_path(wildcards.benchmark, wildcards.exclusion, idx)
        for idx in range(len(entry["files"]))
    ]


def get_exclusion_type(wildcards) -> str:
    """Get type (single/pair) for an exclusion."""
    entry = get_exclusion_entry(wildcards.benchmark, wildcards.exclusion)
    return entry["type"]


def get_exclusion_file_url(benchmark: str, exclusion_name: str, file_idx: int) -> str:
    """
    Get URL for an exclusion file by index.

    Args:
        benchmark: Name of the benchmark set
        exclusion_name: Name of the exclusion category
        file_idx: Index of the file in the files array (0-based)
        ref_name: Reference name in config

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


def get_exclusion_file_checksum(benchmark: str, exclusion_name: str, file_idx: int) -> str:
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


def get_exclusion_name_mapping(benchmark: str) -> Dict[str, str]:
    """Map EXCL_* IDs (as in REGION_IDS) back to canonical exclusion names.

    Returns:
        Dict mapping e.g. "EXCL_FLANKS" -> "flanks"
    """
    exclusions = get_exclusion_config(benchmark)
    return {
        f"EXCL_{_format_exclusion_name(excl['name'])}": excl["name"] for excl in exclusions
    }


def get_exclusion_bed_paths(benchmark: str) -> List[str]:
    """Get list of materialized exclusion BED file paths for a benchmark."""
    exclusions = get_exclusion_config(benchmark)
    return [f"results/exclusions/{benchmark}/{excl['name']}.bed" for excl in exclusions]


def get_exclusion_impact_inputs(wildcards):
    """Get inputs for the compute_exclusion_impact rule."""
    benchmark = wildcards.benchmark
    excl_names = [e["name"] for e in get_exclusion_config(benchmark)]
    return {
        "variant_table": f"results/variant_tables/{benchmark}/variants.tsv",
        "exclusion_tsvs": expand(
            "results/exclusions/{benchmark}/coverage/{exclusion}.tsv",
            benchmark=benchmark,
            exclusion=excl_names,
        ),
    }


def get_exclusion_interaction_inputs(wildcards):
    """Get inputs for the compute_exclusion_interactions rule."""
    benchmark = wildcards.benchmark
    return {
        "dip_bed": f"resources/benchmarksets/{benchmark}_dip.bed",
        "benchmark_bed": f"resources/benchmarksets/{benchmark}_benchmark.bed",
        "exclusion_beds": get_exclusion_bed_paths(benchmark),
        "variant_table": f"results/variant_tables/{benchmark}/variants.tsv",
    }


def get_old_benchmark_analysis_inputs(wildcards):
    """Get inputs for the annotate_old_benchmark_status rule."""
    comp = config["comparisons"][wildcards.comp_id]
    new_bmk = comp["new_benchmark"]
    old_bmk = comp["old_benchmark"]
    return {
        "old_vcf": f"resources/benchmarksets/{old_bmk}_benchmark.vcf.gz",
        "old_vcfidx": f"resources/benchmarksets/{old_bmk}_benchmark.vcf.gz.tbi",
        "old_bed": f"resources/benchmarksets/{old_bmk}_benchmark.bed",
        "new_dip_bed": f"resources/benchmarksets/{new_bmk}_dip.bed",
        "new_benchmark_bed": f"resources/benchmarksets/{new_bmk}_benchmark.bed",
        "exclusion_beds": get_exclusion_bed_paths(new_bmk),
    }


def get_exclusion_impact_targets(wildcards) -> List[str]:
    """Get all exclusion impact table targets for rule all."""
    targets = []
    for benchmark, conf in config["benchmarksets"].items():
        if conf.get("exclusions"):
            targets.append(f"results/exclusions/{benchmark}/exclusion_impact.csv")
    return targets


def get_exclusion_interaction_targets(wildcards) -> List[str]:
    """Get all exclusion interaction table targets for rule all."""
    targets = []
    for benchmark, conf in config["benchmarksets"].items():
        if conf.get("exclusions"):
            targets.append(f"results/exclusions/{benchmark}/exclusion_interactions.csv")
    return targets


def get_reference_checksum_info(ref_name: str) -> Tuple[str, str]:
    """
    Get checksum value and type for a reference.
    Args:
        ref_name: Name of the reference in config["references"]
    Returns:
        Tuple of checksum_value and checksum_type
    Raises:
        KeyError: If reference not found
        ValueError: If no checksum found
    """
    ref_config = config["references"].get(ref_name)
    if ref_config is None:
        raise KeyError(f"Reference '{ref_name}' not found in config")

    if "sha256" in ref_config:
        return ref_config["sha256"], "sha256"
    if "md5" in ref_config:
        return ref_config["md5"], "md5"

    raise ValueError(f"No checksum (md5 or sha256) found for reference: {ref_name}")


def get_reference_checksum(ref_name: str) -> str:
    """Get checksum value for reference."""
    checksum, _ = get_reference_checksum_info(ref_name)
    return checksum


def get_reference_checksum_type(ref_name: str) -> str:
    """Get checksum type for reference."""
    _, checksum_type = get_reference_checksum_info(ref_name)
    return checksum_type


def get_chromosomes(wildcards) -> str:
    """
    Get list of main chromosomes for the reference genome.

    Args:
        wildcards: Snakemake wildcards with `ref` attribute

    Returns:
        Space-separated string of chromosomes
    """
    # Use string chromosome IDs so join works for GRCh37, then add chr prefix for others
    chroms = [str(i) for i in range(1, 23)] + ["X", "Y"]
    if wildcards.ref == "GRCh37":
        return " ".join(chroms)

        ## Adding chr prefix for non-GRCh37 references
    return " ".join([f"chr{i}" for i in chroms])


# ============================================================================
# Stratification Helpers
# ============================================================================


def get_strat_ids(wildcards):
    """Get list of stratification IDs."""
    return list(config.get("stratifications", {}).keys())


def get_strat_config(ref: str, strat_name: str) -> Dict[str, Any]:
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


def get_stratification_url(wildcards) -> str:
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


def get_stratification_sha256(wildcards) -> str:
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
