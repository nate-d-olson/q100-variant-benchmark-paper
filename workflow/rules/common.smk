# AI Disclosure: This file was modified with assistance from Claude (Anthropic)
# to add exclusion analysis helper functions.
"""
Common helper functions for the Q100 variant benchmark pipeline.

Organized by section:
1. Top-level variables and wildcard constraints
2. Chromosome helpers
3. Genomic context helpers
4. Comparison helpers
5. Exclusion helpers
6. Reference & stratification download helpers
7. Rule-all target generators
"""

from typing import List, Dict, Any, Tuple


# ============================================================================
# Top-Level Variables & Wildcard Constraints
# ============================================================================

# Collect all unique stratification/genomic context names across references
_all_strat_names = set()
for ref_config in config.get("references", {}).values():
    _all_strat_names.update(ref_config.get("stratifications", {}).keys())

# Benchmarks with exclusions configured (used for target generation)
BENCHMARKS_WITH_EXCLUSIONS = [
    bid
    for bid in config["benchmarksets"]
    if config["benchmarksets"][bid].get("exclusions")
]


wildcard_constraints:
    comp_id="[^/]+",
    strat_name="|".join(sorted(_all_strat_names)),
    genomic_context="|".join(sorted(_all_strat_names)),


# ============================================================================
# Chromosome Helpers
# ============================================================================


def get_chromosomes(wildcards) -> str:
    """Get list of chromosomes from config, with ref-specific prefixes."""
    autosomes = [str(c) for c in config["chromosomes"]["autosomes"]]
    sex_chroms = config["chromosomes"]["sex"]
    chroms = autosomes + sex_chroms
    if wildcards.ref == "GRCh37":
        return " ".join(chroms)
    return " ".join([f"chr{c}" for c in chroms])


# ============================================================================
# Genomic Context Helpers
# ============================================================================


def get_stratifications_for_ref(ref: str) -> List[str]:
    """Get list of stratification/genomic context names for a reference."""
    return list(config["references"][ref].get("stratifications", {}).keys())


def get_genomic_context_ids(wildcards) -> List[str]:
    """Get genomic context IDs for a benchmark's reference."""
    ref = config["benchmarksets"][wildcards.benchmark]["ref"]
    return get_stratifications_for_ref(ref)


def get_genomic_context_cov_beds(wildcards) -> List[str]:
    """Get coverage BED file paths for a benchmark's genomic contexts."""
    benchmark = wildcards.benchmark
    ref = config["benchmarksets"][benchmark]["ref"]
    contexts = get_stratifications_for_ref(ref)
    return [
        f"results/genomic_context/{benchmark}/coverage/{ctx}_cov.bed" for ctx in contexts
    ]


def get_genomic_context_bed_specs(wildcards) -> List[str]:
    """Get genomic context BED files with IDs (path:name format) for annotation."""
    ref = config["benchmarksets"][wildcards.benchmark]["ref"]
    contexts = config["references"][ref].get("stratifications", {})
    return [f"resources/stratifications/{ref}_{name}.bed.gz:{name}" for name in contexts]


def get_region_beds(wildcards) -> List[str]:
    """Get benchmark region BED and exclusion BEDs (path:ID format) for annotation."""
    benchmark = wildcards.benchmark
    beds = [f"resources/benchmarksets/{benchmark}_benchmark.bed:BMKREGIONS"]
    for excl in _get_exclusion_config(benchmark):
        excl_id = excl["name"].replace("-", "_").upper()
        for idx in range(len(excl["files"])):
            path = f"resources/exclusions/{benchmark}/{excl['name']}_{idx}.bed"
            beds.append(f"{path}:EXCL_{excl_id}")
    return beds


# ============================================================================
# Comparison Helpers
# ============================================================================


def get_comparison_files(wildcards):
    """Get all input files for a benchmark comparison."""
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


# ============================================================================
# Exclusion Helpers
# ============================================================================


def _get_exclusion_config(benchmark: str) -> List[Dict[str, Any]]:
    """Get exclusion configuration for a benchmark set."""
    return config["benchmarksets"].get(benchmark, {}).get("exclusions", [])


def _get_exclusion_entry(benchmark: str, exclusion_name: str) -> Dict[str, Any]:
    """Get a specific exclusion entry by name."""
    for item in _get_exclusion_config(benchmark):
        if item["name"] == exclusion_name:
            return item
    raise ValueError(f"Exclusion {exclusion_name} not found for {benchmark}")


def get_exclusion_inputs(wildcards) -> List[str]:
    """Get input file paths for an exclusion (matches download_exclusion outputs)."""
    entry = _get_exclusion_entry(wildcards.benchmark, wildcards.exclusion)
    return [
        f"resources/exclusions/{wildcards.benchmark}/{wildcards.exclusion}_{idx}.bed"
        for idx in range(len(entry["files"]))
    ]


def get_exclusion_type(wildcards) -> str:
    """Get type (single/pair) for an exclusion."""
    return _get_exclusion_entry(wildcards.benchmark, wildcards.exclusion)["type"]


def get_exclusion_file_url(benchmark: str, exclusion_name: str, file_idx: int) -> str:
    """Get URL for an exclusion file by index."""
    entry = _get_exclusion_entry(benchmark, exclusion_name)
    return entry["files"][file_idx]["url"]


def get_exclusion_file_checksum(benchmark: str, exclusion_name: str, file_idx: int) -> str:
    """Get SHA256 checksum for an exclusion file by index."""
    entry = _get_exclusion_entry(benchmark, exclusion_name)
    return entry["files"][file_idx]["sha256"]


def get_exclusion_name_mapping(benchmark: str) -> Dict[str, str]:
    """Map EXCL_* IDs back to canonical exclusion names (e.g. EXCL_FLANKS -> flanks)."""
    return {
        f"EXCL_{e['name'].replace('-', '_').upper()}": e["name"]
        for e in _get_exclusion_config(benchmark)
    }


def _get_exclusion_bed_paths(benchmark: str) -> List[str]:
    """Get materialized exclusion BED paths for a benchmark."""
    return [
        f"results/exclusions/{benchmark}/{e['name']}.bed"
        for e in _get_exclusion_config(benchmark)
    ]


def get_exclusion_impact_inputs(wildcards):
    """Get inputs for the compute_exclusion_impact rule."""
    benchmark = wildcards.benchmark
    excl_names = [e["name"] for e in _get_exclusion_config(benchmark)]
    return {
        "variant_table": f"results/variant_tables/{benchmark}/variants.parquet",
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
        "exclusion_beds": _get_exclusion_bed_paths(benchmark),
        "variant_table": f"results/variant_tables/{benchmark}/variants.parquet",
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
        "exclusion_beds": _get_exclusion_bed_paths(new_bmk),
    }


# ============================================================================
# Reference & Stratification Download Helpers
# ============================================================================


def get_reference_checksum(ref_name: str) -> str:
    """Get checksum value for a reference (sha256 preferred, then md5)."""
    ref_config = config["references"].get(ref_name)
    if ref_config is None:
        raise KeyError(f"Reference '{ref_name}' not found in config")
    return ref_config.get("sha256", ref_config.get("md5", ""))


def get_stratification_url(wildcards) -> str:
    """Get URL for a stratification BED file."""
    strat = config["references"][wildcards.ref]["stratifications"][wildcards.strat_name]
    return strat.get("url", "")


def get_stratification_sha256(wildcards) -> str:
    """Get SHA256 checksum for a stratification BED file."""
    strat = config["references"][wildcards.ref]["stratifications"][wildcards.strat_name]
    return strat.get("sha256", "")


# ============================================================================
# Rule-All Target Generators
# ============================================================================


def get_exclusion_impact_targets(wildcards) -> List[str]:
    """Get exclusion impact CSV targets for benchmarks with exclusions."""
    return expand(
        "results/exclusions/{benchmark}/exclusion_impact.csv",
        benchmark=BENCHMARKS_WITH_EXCLUSIONS,
    )


def get_exclusion_interaction_targets(wildcards) -> List[str]:
    """Get exclusion interaction CSV targets for benchmarks with exclusions."""
    return expand(
        "results/exclusions/{benchmark}/exclusion_interactions.csv",
        benchmark=BENCHMARKS_WITH_EXCLUSIONS,
    )


def get_chr8_synteny_targets(wildcards) -> List[str]:
    """Return chr8 synteny figure targets if chr8_synteny is configured."""
    if not config.get("chr8_synteny"):
        return []
    return [
        "results/chr8_synteny/chr8_figure.pdf",
        "results/chr8_synteny/chr8_figure.png",
    ]
