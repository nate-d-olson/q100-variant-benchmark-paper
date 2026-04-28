# API Reference — Snakemake Helper Functions

Helper functions exposed by `workflow/rules/common.smk`. These are imported into
all rule files via the include order in `workflow/Snakefile`.

For R-side loaders, see [pipeline-outputs.md](pipeline-outputs.md) and
[data-dictionary.md](data-dictionary.md).

## Top-level Variables

### `BENCHMARKS_WITH_EXCLUSIONS`

List of benchmark IDs that have an `exclusions:` block in `config.yaml`. Used by
target generators to skip benchmarks without exclusions.

```python
BENCHMARKS_WITH_EXCLUSIONS = [
    bid
    for bid in config["benchmarksets"]
    if config["benchmarksets"][bid].get("exclusions")
]
```

### `wildcard_constraints`

```python
wildcard_constraints:
    comp_id="[^/]+",
    strat_name="|".join(sorted(_all_strat_names)),
    genomic_context="|".join(sorted(_all_strat_names)),
```

`_all_strat_names` is built at config parse time from every reference's
`stratifications` map — when adding a new reference or stratification name, it
is automatically picked up.

## Chromosome Helpers

### `get_chromosomes(wildcards) -> str`

Space-separated chromosome list with the ref-specific prefix:

- `GRCh37` → `1 2 … 22 X Y` (no prefix)
- All other refs → `chr1 chr2 … chr22 chrX chrY`

Reads `config["chromosomes"]["autosomes"]` and `config["chromosomes"]["sex"]`.

## Genomic Context Helpers

### `get_stratifications_for_ref(ref: str) -> List[str]`

List of stratification (genomic context) names configured for a reference, e.g.
`["HP", "MAP", "SD", "SD10kb", "TR", "TR10kb"]`.

### `get_genomic_context_ids(wildcards) -> List[str]`

Same as above, looked up via the benchmark's reference. Used by rules that take
`{benchmark}` as a wildcard.

### `get_genomic_context_cov_beds(wildcards) -> List[str]`

Coverage BED paths for all of a benchmark's genomic contexts:

```
results/genomic_context/{benchmark}/coverage/HP_cov.bed
results/genomic_context/{benchmark}/coverage/MAP_cov.bed
…
```

### `get_genomic_context_bed_specs(wildcards) -> List[str]`

`path:name` specs (one per context) consumed by `combine_beds_with_id.py`. The
suffix after `:` becomes the value written to `INFO/CONTEXT_IDS`.

```
resources/stratifications/GRCh38_HP.bed.gz:HP
resources/stratifications/GRCh38_MAP.bed.gz:MAP
…
```

## Region Helpers

### `get_region_beds(wildcards) -> List[str]`

`path:ID` specs for the benchmark BED plus all configured exclusions. The
`BMKREGIONS` ID marks the benchmark itself; exclusion IDs are derived from the
config name (e.g. `consecutive-svs` → `EXCL_CONSECUTIVE_SVS`).

```
resources/benchmarksets/{benchmark}_benchmark.bed:BMKREGIONS
resources/exclusions/{benchmark}/consecutive-svs_0.bed:EXCL_CONSECUTIVE_SVS
…
```

## Comparison Helpers

### `get_comparison_files(wildcards) -> Dict[str, str]`

Dict of inputs for a `{comp_id}` benchmark comparison: paired VCFs+indexes,
BEDs, and the reference FASTA. Driven by `config["comparisons"][comp_id]` with
keys `new_benchmark`, `old_benchmark`, `ref`.

## Exclusion Helpers

### `get_exclusion_inputs(wildcards) -> List[str]`

Per-file BED paths for `{wildcards.exclusion}`, matching the `download_exclusion`
rule outputs. For multi-file exclusions (`type: pair`) returns one path per file.

### `get_exclusion_type(wildcards) -> str`

Returns `"single"` or `"pair"` from the config entry — controls how
`materialize_exclusion` combines source BEDs.

### `get_exclusion_file_url(benchmark, exclusion_name, file_idx) -> str`
### `get_exclusion_file_checksum(benchmark, exclusion_name, file_idx) -> str`

URL and SHA256 lookups for a specific exclusion file. Called from `download_exclusion`.

### `get_exclusion_name_mapping(benchmark: str) -> Dict[str, str]`

Maps `EXCL_*` IDs back to canonical exclusion names (e.g.
`EXCL_CONSECUTIVE_SVS` → `consecutive-svs`). Used by Python scripts that read
the annotated VCF and need to reconstruct config-level names.

### `get_exclusion_impact_inputs(wildcards) -> Dict[str, Any]`

Inputs for `compute_exclusion_impact`: variant Parquet + per-exclusion coverage
TSVs. Returned as a dict to be unpacked with `unpack(...)`.

### `get_exclusion_interaction_inputs(wildcards) -> Dict[str, Any]`

Inputs for `compute_exclusion_interactions`: dip.bed, benchmark.bed, all
materialized exclusion BEDs, and the variant Parquet. Returned as a dict.

### `get_old_benchmark_analysis_inputs(wildcards) -> Dict[str, Any]`

Inputs for `annotate_old_benchmark_status`: old benchmark VCF/BED plus the new
benchmark's dip.bed, benchmark.bed, and exclusion BEDs. Looks up
`config["comparisons"][comp_id]` for `new_benchmark` / `old_benchmark`.

### Private helpers

- `_get_exclusion_config(benchmark)` — list of exclusion entries from config
- `_get_exclusion_entry(benchmark, name)` — single exclusion entry by name
- `_get_exclusion_bed_paths(benchmark)` — materialized BED paths for all exclusions

Underscore prefix marks them as not part of the documented surface; rules should
use the public wrappers above.

## Reference & Stratification Download Helpers

### `get_reference_checksum(ref_name: str) -> str`

Returns SHA256 (preferred) or MD5 from the reference config. Used by
`download_reference`.

### `get_stratification_url(wildcards) -> str`
### `get_stratification_sha256(wildcards) -> str`

URL / SHA256 for a stratification BED. Used by `download_stratification`.

## Rule-All Target Generators

### `get_exclusion_impact_targets(wildcards) -> List[str]`

```
results/exclusions/{benchmark}/exclusion_impact.csv
```

for each benchmark in `BENCHMARKS_WITH_EXCLUSIONS`.

### `get_exclusion_interaction_targets(wildcards) -> List[str]`

```
results/exclusions/{benchmark}/exclusion_interactions.csv
```

for each benchmark in `BENCHMARKS_WITH_EXCLUSIONS`.

### `get_chr8_synteny_targets(wildcards) -> List[str]`

Returns `chr8_figure.pdf` + `chr8_figure.png` if `chr8_synteny:` is present in
config, else an empty list. Allows the chr8 pipeline to be opt-in.

## Adding a Helper

1. Add the function to the appropriate section header in `common.smk`
2. Use `_underscore_prefix` for private helpers
3. Add a one-line docstring (used by `snakemake --list` and IDE tooltips)
4. Document publicly-used helpers in this file
5. If the helper changes a wildcard pattern, update `wildcard_constraints` too

## Related

- [Architecture Overview](architecture.md) — module-level layout
- [Pipeline Outputs Reference](pipeline-outputs.md) — what the rules produce
- `workflow/rules/common.smk` — source of truth
