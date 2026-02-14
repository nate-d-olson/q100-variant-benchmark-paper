# Q100 Variant Benchmark — Project-Specific Guidance

## Project Overview

Quarto manuscript analyzing the GIAB Q100 HG002 variant benchmark. The Snakemake pipeline processes variant calls across multiple reference genomes (GRCh37, GRCh38, CHM13v2.0) and benchmark versions (v0.6, v4.2.1, v5.0q).

## Key Directories

- `R/` - R source files (data loading, schemas, caching)
- `analysis/` - Quarto notebooks and cached data
- `config/` - Pipeline configuration (config.yaml)
- `docs/` - Architecture docs, data dictionary, troubleshooting
- `tests/` - R tests (testthat) and Python tests (pytest)
- `workflow/` - Snakemake rules and Python scripts
- `results/` - Pipeline outputs (gitignored, must be generated)
- `resources/` - Benchmark BED files (gitignored, must be generated)

## Pipeline Organization

**Terminology:**

- **Genomic Context**: Difficult regions (HP, TR, SD, MAP, etc.)—used for variant analysis and metrics computation
- **Stratification**: Same regions as genomic context, but when used in GIAB comparison analysis or benchmarkset comparisons
- **Benchmark Regions**: The set of regions where variants are confident—separate from genomic contexts
- **Bench Type**: Classification for benchmark variants (smvar = small variants <50bp, stvar = structural variants >=50bp)
- **Exclusions**: Regions removed from dip.bed to produce the final v5.0q benchmark regions. Configured per-benchmark in `config.yaml`. Names must be consistent across all v5.0q benchmarks (e.g., use `gaps` not `gaps_slop`).
- **dip.bed**: Dipcall assembly regions — the starting point before exclusions are applied

**Key Rules by Layer:**

1. **Metrics** (`workflow/rules/genomic_context_metrics.smk`): Compute overlap statistics between genomic contexts and benchmarks
   - Outputs: `results/genomic_context/{benchmark}/genomic_context_coverage_table.csv`
2. **Variant Tables** (`workflow/rules/var_tables.smk`): Annotate VCFs with genomic context IDs and generate Parquet
   - Key rule: `annotate_vcf_genomic_contexts` → adds INFO/CONTEXT_IDS to VCF
   - `generate_variant_parquet` → reads annotated VCF with Truvari, classifies types, writes Parquet
   - Outputs: `results/variant_tables/{benchmark}/variants.parquet`
3. **Variant Counts** (`workflow/rules/var_counts.smk`): Count variants per genomic context
   - Reads Parquet variant table, explodes context_ids, groups by (context_name, var_type, szbin)
   - Outputs: `results/var_counts/{benchmark}/variants_by_genomic_context.parquet`
4. **Comparisons** (`workflow/rules/benchmark_comparisons.smk`): Truvari comparison between benchmark versions
   - Uses "stratification" terminology (not genomic_context) for GIAB comparison analysis
5. **Exclusions** (`workflow/rules/exclusions.smk`): Exclusion analysis for v5.0q benchmarks only
   - `materialize_exclusion`: Downloads/merges exclusion BED files
   - `compute_exclusion_metrics`: Per-exclusion BED overlap with dip.bed (uses `compute_bed_metrics.py`)
   - `compute_exclusion_impact`: Per-exclusion variant counts + BED metrics (Q1)
   - `compute_exclusion_interactions`: Upset-style exclusion combination analysis (Q2)
   - `annotate_old_benchmark_status`: Cross-version comparison of old vs v5.0q benchmarks (Q3)
   - Outputs: `results/exclusions/{benchmark}/exclusion_impact.csv`, `exclusion_interactions.csv`
   - Cross-version outputs: `results/exclusions/{comp_id}/old_only_*.csv`

**Output File Patterns:**

- Genomic context files: `results/genomic_context/{benchmark}/`
  - `coverage/` - BED files (symlinked `.bed.gz` and bedtools `_cov.bed`)
  - `metrics/` - Per-context metrics TSV files
  - `genomic_context_coverage_table.csv` - Aggregated coverage metrics
  - `sizes/` - Temporary size calculation files (auto-deleted by Snakemake)
- Variant tables: `results/variant_tables/{benchmark}/variants.parquet`
- Variant counts: `results/var_counts/{benchmark}/variants_by_genomic_context.parquet`
- Exclusions (v5.0q only): `results/exclusions/{benchmark}/*.csv`
- Cross-version exclusion analysis: `results/exclusions/{comp_id}/*.csv`

**Snakemake Environment:**

- Activate with: `micromamba activate q100-smk`
- Dry-run: `snakemake -n <target>`
- Snakemake version: 9.x (uses `--` before positional targets with flags like `-n`)

## Common Debugging Patterns

**bcftools annotation error**: "The INFO tag 'CONTEXT_IDS' is not defined"

- **Cause**: Header file was cached with old definitions before script changes
- **Root cause**: `generate_header_lines.py` was updated but `results/generate_annotation_headers/` had stale cache
- **Fix**: Delete affected output directories to force regeneration

  ```bash
  rm -rf results/generate_annotation_headers/
  rm -rf results/annotate_vcf_genomic_contexts/
  rm -rf results/annotate_vcf_regions/
  rm -rf results/combine_genomic_context_beds/
  rm -rf results/variant_tables/
  ```

- **Prevention**: Always clear downstream outputs when modifying scripts that generate headers/annotations

**Wildcard constraint errors**: Check `workflow/rules/common.smk` for `wildcard_constraints`

- Must match all wildcards used in affected rules
- Format: `strat_name="|".join(sorted(_all_strat_names))`
- For new genomic context rules: use `get_genomic_context_ids()` helper function

**Script changes not taking effect**:

- Snakemake does not automatically detect when a script changes
- Solution: Delete the output directory that depends on the script; re-run to regenerate
- Example: modify `combine_beds_with_id.py` → delete `results/combine_genomic_context_beds/`

**Empty variant counts with malformed context IDs**: "Contexts: ["'MAP'", "'MAP')", ...]"

- **Cause**: Truvari's `vcf_to_df()` converts multi-value VCF INFO fields into Python tuples/lists
- **Root cause**: When written to Parquet, these become string representations like `"('HP', 'MAP')"`
- **Fix**: Added `normalize_annotation()` function in `generate_variant_parquet.py` (line 66) for tuple-to-string conversion
- **Fixed**: February 2026
- **Verification**: Check logs show clean context names: `Contexts: ['HP', 'MAP', 'SD', ...]`
- **Impact**: Affected all variant count tables (genomic contexts and exclusions) across all benchmarks

## R Data Loading Infrastructure

### Data Dictionary

For a detailed description of the data objects returned by the loading functions and pipeline metric definitions, see:
[`docs/data-dictionary.md`](docs/data-dictionary.md)

> **Maintenance Note:** The data dictionary document must be updated whenever changes are made to `R/data_loading.R`, `R/schemas.R`, or pipeline output formats.

### File Organization

```
R/
├── schemas.R       # Arrow schema registry, factor levels, validation rules
├── cache.R         # Parquet caching (sources schemas.R)
└── data_loading.R  # All loading functions (sources schemas.R and cache.R)
```

Source order: `data_loading.R` sources `schemas.R` and `cache.R` at the top. Use `source(here::here("R/data_loading.R"))` to get everything.

### Caching System

- **Format:** Parquet with zstd compression (level 3) via Arrow
- **Location:** `analysis/cache/` (configurable via `getOption("q100.cache_dir")`)
- **Invalidation:** Cache key = `rlang::hash()` of dataset name + source file mtimes + params
- **Metadata:** Pipeline metadata (R version, package versions, config) stored as JSON in Parquet file-level key-value metadata under key `q100_pipeline_metadata`
- **Factors:** Stripped to character before Parquet write, restored on read from schema registry

### Adding a New Cached Dataset

1. `R/schemas.R`: Add entries to `get_arrow_schema()`, `get_factor_levels()`, `get_validation_rules()`
2. `R/data_loading.R`: Add `use_cache`/`force_refresh` params, call `read_cache()` before loading, `write_cache()` after

### Testing

```bash
Rscript -e 'testthat::test_file("tests/test_cache.R")'      # Schema + cache tests (45 tests)
Rscript -e 'testthat::test_file("tests/test_data_loading.R")'  # Data loading tests (has pre-existing failures)
```

Cache tests use `withr::local_options(q100.cache_dir = tempdir)` for isolation.

### Known Issues

- `parse_benchmark_id()` regex captures `"5.0q"` not `"v5.0q"` -- tests expect `"v"` prefix
- `test_data_loading.R` references old column names (`strat_name`, `strat_filter`) from pre-refactor
- Tests that load real data fail without `results/` directory (pipeline outputs not in repo)

## Code Quality

- **Formatter:** `air` (config in `air.toml`: 100-char lines, 2-space indent)
- **Linter:** `lintr` (config in `.lintr`: allows `dotted.case` for internal helpers)
- **Pipe Operator:** Use `%>%` (magrittr pipe) instead of `|>` (native pipe) for consistency with existing codebase

## Important Patterns

- `fs::dir_ls` glob matches against **full path** via `glob2rx` -- use `grepl()` on `fs::path_file()` for filename-only matching
- Factor levels are centralized in `R/schemas.R` constants (e.g., `BENCH_VERSION_LEVELS`, `CHROM_LEVELS`)
- All loading functions have `use_cache = TRUE` and `force_refresh = FALSE` defaults; cache writes are wrapped in `tryCatch` so failures don't break loading

**Quarto notebooks:**

- `analysis/benchmarkset_characterization.qmd` — Primary analysis; loads `variants_df` and `genomic_context_metrics_df`
- `analysis/benchmark_difficult.qmd` — Coverage analysis; loads `diff_cov_df`
- `analysis/external_evaluation.qmd` — External benchmark comparisons

## Column Naming Conventions

**Genomic Context Files** (CSV/TSV):

- `context_name` — Genomic context identifier (HP, TR, SD, MAP, etc.)
- `context_bp` — Total size of genomic context
- `intersect_bp` — Overlap with benchmark
- `pct_of_context` — % of genomic context covered by benchmark
- `pct_of_bench` — % of benchmark within this genomic context
- `total_variants` — Total variant count
- `variant_density_per_mb` — Variants per Mb of benchmark-covered region

**VCF Annotation Fields** (after bcftools annotate):

- INFO/CONTEXT_IDS — Comma-separated genomic context region IDs overlapping variant
- INFO/REGION_IDS — Comma-separated benchmark + exclusion region IDs

**Variant Table** (`variants.parquet`):

- Generated by `generate_variant_parquet.py` using Truvari's VariantRecord API
- Columns: bench_version, ref, bench_type, chrom, pos, end, gt, var_type, var_size, szbin, ref_len, alt_len, qual, filter, is_pass, context_ids, region_ids

## Benchmark ID Format

Pattern: `{bench_version}_{ref}_{bench_type}`

Examples:

- `v5.0q_GRCh38_smvar` — v5.0q small variants on GRCh38
- `v5.0q_GRCh37_stvar` — v5.0q structural variants on GRCh37
- `v4.2.1_GRCh38_smvar` — v4.2.1 small variants on GRCh38

File paths encode this pattern:

```
results/genomic_context/v5.0q_GRCh38_smvar/
results/variant_tables/v5.0q_GRCh38_smvar/variants.parquet
results/var_counts/v5.0q_GRCh38_smvar/variants_by_genomic_context.parquet
results/exclusions/v5.0q_GRCh38_smvar/  # only if exclusions configured
```

## Variant Type Classification

The pipeline uses Truvari's VariantRecord API in `workflow/scripts/generate_variant_parquet.py`:

**Variant Type (`var_type`):**

- Uses `VariantRecord.var_type()` which returns an `SV` enum object
- Enum values: SNP, DEL, INS, DUP, INV, BND, NON, UNK
- Extracted to string using `.name` (e.g., `SV.DEL` → `"DEL"`)
- NON (non-variant) records are filtered out

**Genotype (`gt`):**

- Uses `VariantRecord.gt()` which returns a tuple (e.g., `(0, 1)`)
- Converted to `GT` enum using `truvari.get_gt()`
- Enum values: REF, HET, HOM, NON, UNK
- Extracted to string using `.name` (e.g., `GT.HET` → `"HET"`)

**Size Bin (`szbin`):**

- Uses `truvari.get_sizebin(var_size)` which returns a string directly
- Values: `"SNP"`, `"[1,5)"`, `"[50,100)"`, `">=5k"`, etc.

**Size Filtering:**

- smvar: keeps variants with abs(var_size) < 50bp
- stvar: keeps variants with abs(var_size) >= 50bp

**Pipeline Flow:**

1. `annotate_vcf_genomic_contexts` — bcftools annotate adds CONTEXT_IDS/REGION_IDS to VCF
2. `generate_variant_parquet` — Truvari VariantRecord reads annotated VCF, classifies types, writes Parquet
3. `count_variants_by_genomic_context` — Reads Parquet, counts by (context_name, var_type, szbin)
