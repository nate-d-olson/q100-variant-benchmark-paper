# API Reference - Helper Functions

This document provides detailed reference for helper functions in the Q100 variant benchmark pipeline, primarily from `workflow/rules/common.smk`.

## Table of Contents

- [Benchmark Configuration Functions](#benchmark-configuration-functions)
- [Exclusion Management Functions](#exclusion-management-functions)
- [Stratification Management Functions](#stratification-management-functions)
- [Region Management Functions](#region-management-functions)
- [Reference Genome Functions](#reference-genome-functions)
- [Input Aggregation Functions](#input-aggregation-functions)

---

## Benchmark Configuration Functions

### `get_exclusion_config(benchmark)`

Get exclusion configuration for a benchmark set.

**Parameters:**
- `benchmark` (str): Benchmark set name (e.g., "v5q_GRCh38")

**Returns:**
- `list`: List of exclusion dictionaries with keys: `name`, `type`, `files`

**Example:**
```python
exclusions = get_exclusion_config("v5q_GRCh38")
# Returns: [
#     {"name": "consecutive-svs", "type": "single", "files": [{"url": "...", "sha256": "..."}]},
#     {"name": "lc-regions", "type": "single", "files": [...]},
#     ...
# ]
```

**Notes:**
- Returns empty list if benchmark has no exclusions configured
- Used by other exclusion functions to access config

---

## Exclusion Management Functions

### `get_exclusion_items(wildcards)`

Get list of exclusion names for a benchmark set.

**Parameters:**
- `wildcards`: Snakemake wildcards object with `benchmark` attribute

**Returns:**
- `list[str]`: List of exclusion names

**Example:**
```python
# In Snakemake rule:
rule process_exclusions:
    input:
        lambda wildcards: get_exclusion_items(wildcards)
# Returns: ["consecutive-svs", "lc-regions", "tandem-repeats-strs"]
```

---

### `get_exclusion_entry(benchmark, exclusion_name)`

Get a specific exclusion entry by name.

**Parameters:**
- `benchmark` (str): Benchmark set name
- `exclusion_name` (str): Name of the exclusion category

**Returns:**
- `dict`: Exclusion entry with keys: `name`, `type`, `files`

**Raises:**
- `ValueError`: If exclusion_name not found for the benchmark

**Example:**
```python
entry = get_exclusion_entry("v5q_GRCh38", "consecutive-svs")
# Returns: {"name": "consecutive-svs", "type": "single", "files": [...]}
```

---

### `get_exclusion_file_path(benchmark, exclusion_name, file_idx)`

Get the standardized local path for an exclusion file.

**Parameters:**
- `benchmark` (str): Benchmark set name
- `exclusion_name` (str): Name of the exclusion category
- `file_idx` (int): Index of the file (0-based)

**Returns:**
- `str`: Local file path where downloaded exclusion file will be stored

**Example:**
```python
path = get_exclusion_file_path("v5q_GRCh38", "consecutive-svs", 0)
# Returns: "resources/exclusions/v5q_GRCh38/consecutive-svs_0.bed"
```

**Notes:**
- Used by download rules and input functions
- Ensures consistent naming convention

---

### `get_exclusion_inputs(wildcards)`

Get input file paths for an exclusion.

**Parameters:**
- `wildcards`: Snakemake wildcards with `benchmark` and `exclusion` attributes

**Returns:**
- `list[str]`: List of file paths for exclusion BED files

**Example:**
```python
# For exclusion with 2 files:
paths = get_exclusion_inputs(wildcards)
# Returns: [
#     "resources/exclusions/v5q_GRCh38/consecutive-svs_0.bed",
#     "resources/exclusions/v5q_GRCh38/consecutive-svs_1.bed"
# ]
```

---

### `get_exclusion_type(wildcards)`

Get type (single/pair) for an exclusion.

**Parameters:**
- `wildcards`: Snakemake wildcards with `benchmark` and `exclusion` attributes

**Returns:**
- `str`: "single" or "pair"

**Example:**
```python
excl_type = get_exclusion_type(wildcards)
# Returns: "single" or "pair"
```

---

## Stratification Management Functions

### `get_stratification_beds(wildcards)`

Get list of stratification BED files with IDs.

**Parameters:**
- `wildcards`: Snakemake wildcards with `benchmark` attribute

**Returns:**
- `list[str]`: List of BED file paths with IDs in "path:ID" format

**Example:**
```python
beds = get_stratification_beds(wildcards)
# Returns: [
#     "resources/stratifications/GRCh38_TR.bed.gz:TR",
#     "resources/stratifications/GRCh38_HOMOPOLYMERS_7BP.bed.gz:HOMOPOLYMERS_7BP",
#     "resources/stratifications/GRCh38_SEGDUPS.bed.gz:SEGDUPS"
# ]
```

**Format:**
- Each entry is `{file_path}:{strat_id}`
- Used by annotation rules to add stratification IDs

---

### `get_strat_ids(wildcards)`

Get list of stratification IDs for a benchmark.

**Parameters:**
- `wildcards`: Snakemake wildcards with `benchmark` attribute

**Returns:**
- `list[str]`: List of stratification IDs (short names)

**Example:**
```python
ids = get_strat_ids(wildcards)
# Returns: ["TR", "HOMOPOLYMERS_7BP", "SEGDUPS", "LOWMAP", ...]
```

**Alias:**
- `get_stratification_ids(wildcards)` - Same functionality, different name

---

## Region Management Functions

### `get_region_beds(wildcards)`

Get benchmark region BED and exclusion BEDs with IDs.

**Parameters:**
- `wildcards`: Snakemake wildcards with `benchmark` attribute

**Returns:**
- `list[str]`: List of BED file paths with IDs in "path:ID" format

**Example:**
```python
beds = get_region_beds(wildcards)
# Returns: [
#     "resources/benchmarksets/v5q_GRCh38_benchmark.bed:BMKREGIONS",
#     "resources/exclusions/v5q_GRCh38/consecutive-svs_0.bed:EXCL_CONSECUTIVE_SVS",
#     "resources/exclusions/v5q_GRCh38/lc-regions_0.bed:EXCL_LC_REGIONS",
#     ...
# ]
```

**Format:**
- Benchmark regions always have ID "BMKREGIONS"
- Exclusions have ID format "EXCL_{EXCLUSION_NAME}"
- Exclusion names converted: "consecutive-svs" → "CONSECUTIVE_SVS"

**Logic:**
1. Always includes benchmark regions BED
2. For v5q benchmarks, adds all configured exclusion BEDs
3. Each exclusion file gets unique ID suffix if multiple files

**Notes:**
- Central function for combining benchmark and exclusion regions
- Used by annotation pipeline to track which regions variants fall in

---

### `get_region_ids(wildcards)`

Get list of region IDs (benchmark + exclusions).

**Parameters:**
- `wildcards`: Snakemake wildcards with `benchmark` attribute

**Returns:**
- `list[str]`: List of region ID strings

**Example:**
```python
ids = get_region_ids(wildcards)
# Returns: [
#     "BMKREGIONS",
#     "EXCL_CONSECUTIVE_SVS",
#     "EXCL_LC_REGIONS",
#     "EXCL_TANDEM_REPEATS_STRS"
# ]
```

**Notes:**
- Matches the IDs used in `get_region_beds()`
- Used by `expand_annotations.py` to create binary columns

---

## Reference Genome Functions

### `get_reference_checksum(ref_name)`

Get checksum value for a reference genome.

**Parameters:**
- `ref_name` (str): Name of the reference in `config["references"]`

**Returns:**
- `str`: Checksum string (MD5 or SHA256)

**Raises:**
- `ValueError`: If no checksum found for reference

**Example:**
```python
checksum = get_reference_checksum("GRCh38")
# Returns: "abc123def456..." (SHA256 or MD5 value)
```

---

### `get_reference_checksum_type(ref_name)`

Determine checksum type for a reference genome.

**Parameters:**
- `ref_name` (str): Name of the reference in `config["references"]`

**Returns:**
- `str`: "md5" or "sha256"

**Raises:**
- `ValueError`: If no checksum found for reference

**Example:**
```python
checksum_type = get_reference_checksum_type("GRCh38")
# Returns: "sha256" or "md5"
```

**Notes:**
- Checks for SHA256 first, then MD5
- Used by download rules to select appropriate hash algorithm

---

## Input Aggregation Functions

### `get_var_table_inputs(wildcards)`

Generate list of variant table files for all benchmarks.

**Parameters:**
- `wildcards`: Snakemake wildcards (unused but required by Snakemake)

**Returns:**
- `list[str]`: List of variant table file paths

**Example:**
```python
files = get_var_table_inputs(wildcards)
# Returns: [
#     "results/variant_tables/v5q_GRCh38/variants.tsv",
#     "results/variant_tables/v5q_GRCh37/variants.tsv",
#     "results/variant_tables/v421_GRCh38/variants.tsv",
#     ...
# ]
```

---

### `get_exclusion_table_inputs(wildcards)`

Generate list of exclusion intersection tables for all benchmarks that have exclusions configured.

**Parameters:**
- `wildcards`: Snakemake wildcards (unused but required by Snakemake)

**Returns:**
- `list[str]`: List of exclusion table file paths

**Example:**
```python
files = get_exclusion_table_inputs(wildcards)
# Returns: [
#     "results/exclusions/v5q_GRCh38/exclusions_intersection_table.csv",
#     "results/exclusions/v5q_GRCh37/exclusions_intersection_table.csv",
#     ...
# ]
```

**Notes:**
- Only includes benchmarks with `"exclusions"` in config
- Returns empty list if no benchmarks have exclusions

---

### `get_strat_metrics_inputs(wildcards)`

Generate list of stratification metrics table files for all benchmarks.

**Parameters:**
- `wildcards`: Snakemake wildcards (unused but required by Snakemake)

**Returns:**
- `list[str]`: List of stratification metrics file paths

**Example:**
```python
files = get_strat_metrics_inputs(wildcards)
# Returns: [
#     "results/strat_metrics/v5q_GRCh38/stratification_coverage_table.csv",
#     "results/strat_metrics/v5q_GRCh37/stratification_coverage_table.csv",
#     ...
# ]
```

**Notes:**
- Only includes benchmarks with `"dip_bed"` configured
- Stratification metrics require dip.bed file for computation

---

### `get_var_counts_inputs(wildcards)`

Generate list of variant count table files for all benchmarks.

**Parameters:**
- `wildcards`: Snakemake wildcards (unused but required by Snakemake)

**Returns:**
- `list[str]`: List of variant count file paths

**Example:**
```python
files = get_var_counts_inputs(wildcards)
# Returns: [
#     "results/var_counts/v5q_GRCh38/stratification_combined_metrics.csv",
#     "results/var_counts/v5q_GRCh37/stratification_combined_metrics.csv",
#     ...
# ]
```

**Notes:**
- Only includes benchmarks with `"dip_bed"` configured
- Variant counts depend on stratification metrics

---

## Common Patterns

### Path Construction Pattern

Many functions follow this pattern for building standardized paths:

```python
# General pattern:
"{base_dir}/{entity_type}/{identifier}/{file_name}"

# Examples:
"resources/benchmarksets/{benchmark}_benchmark.bed"
"resources/exclusions/{benchmark}/{exclusion}_{idx}.bed"
"resources/stratifications/{ref}_{strat}.bed.gz"
"results/var_tables/{benchmark}/{ref}/annotated.tsv"
```

### ID Transformation Pattern

Exclusion names undergo consistent transformation:

```python
# Pattern: dash-separated → underscore-uppercase
"consecutive-svs" → "CONSECUTIVE_SVS"
"lc-regions" → "LC_REGIONS"

# Implementation:
name.replace("-", "_").upper()
```

### Config Lookup Pattern

Helper functions consistently access config with safe defaults:

```python
# Safe access with get():
config["benchmarksets"].get(benchmark, {}).get("exclusions", [])

# Raises KeyError if required:
config["benchmarksets"][benchmark]["ref"]
```

---

## Usage in Snakemake Rules

### Input Function Example

```python
rule combine_region_beds:
    input:
        beds=get_region_beds  # Function reference (no parentheses)
    output:
        "results/{benchmark}/combined_regions.bed"
    shell:
        "cat {input.beds} > {output}"
```

### Expand with Helper Function

```python
rule aggregate:
    input:
        expand(
            "results/{benchmark}/output.txt",
            benchmark=[b for b in config["benchmarksets"]]
        )
```

### Lambda Wrapper for Complex Logic

```python
rule process:
    input:
        lambda wildcards: get_exclusion_inputs(wildcards)
    output:
        "results/{benchmark}/{exclusion}/processed.bed"
```

---

## Debugging Helper Functions

### Print Function Output

```python
# In Python script or interactive session:
from snakemake import Wildcards

wildcards = Wildcards(benchmark="v5q_GRCh38", exclusion="consecutive-svs")
print(get_exclusion_inputs(wildcards))
```

### Test Config Access

```python
# Load config and test functions:
import yaml
config = yaml.safe_load(open("config/config.yaml"))

# Test function:
print(get_region_beds(Wildcards(benchmark="v5q_GRCh38")))
```

### Dry-run to See Resolved Paths

```bash
# See what paths will be used:
snakemake --dry-run --printshellcmds target_rule

# See full DAG with resolved wildcards:
snakemake --dag target_rule | dot -Tpdf > dag.pdf
```

---

*Last Updated: 2026-01-13*
*Branch: feature/codebase-improvements*
