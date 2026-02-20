# Q100 Variant Benchmark — Project-Specific Guidance

## Project Overview

Quarto manuscript analyzing the GIAB Q100 HG002 variant benchmark. The Snakemake pipeline processes variant calls across multiple reference genomes (GRCh37, GRCh38, CHM13v2.0) and benchmark versions (v0.6, v4.2.1, v5.0q).

## Key Directories

- `R/` - R source files (data loading, schemas, caching)
- `analysis/` - Quarto notebooks and cached data
- `config/` - Pipeline configuration (config.yaml)
- `docs/` - Architecture docs, data dictionary, troubleshooting
- `manuscript/` - Quarto manuscript chapters (introduction, methods, results, discussion)
- `scripts/` - Utility scripts (create_grch38_debug_subset.py, happy_giab.R)
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

1. **Genomic Context Analysis** (`workflow/rules/genomic_context_analysis.smk`): Coverage computation, metrics, variant table generation, and variant counting
   - `genomic_context_coverage`: bedtools coverage per context
   - `compute_genomic_context_coverage_table`: uses `compute_coverage_table.py` to process all contexts
   - `generate_variant_parquet`: reads annotated VCF with Truvari, classifies types, writes Parquet
   - `count_variants_by_genomic_context`: counts by (context_name, var_type, szbin)
   - Outputs: `results/genomic_context/{benchmark}/genomic_context_coverage_table.csv`, `results/variant_tables/{benchmark}/variants.parquet`, `results/genomic_context/{benchmark}/variants_by_genomic_context.parquet`
2. **Annotation** (`workflow/rules/annotation.smk`): Annotate VCFs with genomic context and region IDs
   - `combine_genomic_context_beds` / `combine_region_beds`: merge BEDs for bcftools
   - `generate_annotation_headers`: inline heredoc for VCF INFO headers
   - `annotate_vcf_genomic_contexts` → adds INFO/CONTEXT_IDS to VCF
   - `annotate_vcf_regions` → adds INFO/REGION_IDS to VCF
3. **Comparisons** (`workflow/rules/benchmark_comparisons.smk`): Truvari comparison between benchmark versions
   - Uses "stratification" terminology (not genomic_context) for GIAB comparison analysis
4. **Exclusions** (`workflow/rules/exclusions.smk`): Exclusion analysis for v5.0q benchmarks only
   - `materialize_exclusion`: Sort/merge exclusion BED files (NOT symlinks—performs data processing)
   - `compute_exclusion_metrics`: Per-exclusion BED overlap with dip.bed (inline bedtools shell commands)
   - `compute_exclusion_impact`: Per-exclusion variant counts + BED metrics (Q1)
   - `compute_exclusion_interactions`: Upset-style exclusion combination analysis (Q2)
   - `annotate_old_benchmark_status`: Cross-version comparison of old vs v5.0q benchmarks (Q3)
   - Outputs: `results/exclusions/{benchmark}/exclusion_impact.csv`, `exclusion_interactions.csv`
   - Cross-version outputs: `results/exclusions/{comp_id}/old_only_*.csv`
   - **Note:** Exclusions pipeline has different architecture than genomic contexts (see `docs/exclusions-pipeline-evaluation.md`)

**Output File Patterns:**

- Genomic context files: `results/genomic_context/{benchmark}/`
  - `coverage/` - bedtools coverage output (`_cov.bed`) from `genomic_context_analysis.smk`
  - `genomic_context_coverage_table.csv` - Coverage metrics (from `compute_coverage_table.py`)
- Variant tables: `results/variant_tables/{benchmark}/variants.parquet`
- Variant counts: `results/genomic_context/{benchmark}/variants_by_genomic_context.parquet`
- Exclusions (v5.0q only): `results/exclusions/{benchmark}/*.csv`
- Cross-version exclusion analysis: `results/exclusions/{comp_id}/*.csv`

**Snakemake Environment:**

- Activate with: `micromamba activate q100-smk`
- Dry-run: `snakemake -n <target>`
- Snakemake version: 8.x (`min_version("8.0")` in Snakefile)

**Makefile Shortcuts** (preferred over raw snakemake commands):

```bash
make dry-run       # Validate workflow (snakemake -n --quiet)
make lint          # Snakemake + Python (ruff) + R (lintr) linting
make format        # Format Python (ruff), Snakemake (snakefmt), Markdown, R (styler)
make format-check  # Check formatting without modifying files
make test          # lint + format-check + dry-run
make run           # Execute pipeline (--cores 14 --sdm conda)
make dag           # Generate pipeline DAG PDF (results/dag/pipeline.pdf)
make clean         # Remove logs/, .snakemake/, results/, analysis/cache/
```

**Conda Channel Configuration:**

Anaconda CDN (`conda.anaconda.org`, `repo.anaconda.com`) is blocked by the organization firewall.
Channels are redirected to prefix.dev via `~/.condarc`:

```yaml
custom_multichannels:
  conda-forge:
    - https://prefix.dev/conda-forge
  bioconda:
    - https://prefix.dev/bioconda
default_channels: []
```

This remaps `conda-forge` and `bioconda` channel names used in all `workflow/envs/*.yaml` files.
**Note**: `channel_alias` does not work for prefix.dev because its URL path (`/conda-forge`) differs
from the conda channel alias convention (`/channels/conda-forge`). Use `custom_multichannels` instead.

**GitHub CLI TLS issue:**

The `gh` CLI intermittently fails with `tls: failed to verify certificate: x509: OSStatus -26276` due to the organization network proxy. When this occurs, use the MCP GitHub tools (e.g., `mcp__plugin_github_github__merge_pull_request`) as a fallback for PR operations. Git push/pull over SSH is unaffected.

## Chr8 Synteny Figure Pipeline

The chr8 synteny pipeline (`workflow/rules/chr8_synteny.smk`) produces a multi-panel PDF/PNG figure
showing GRCh38 chr8 vs HG002 maternal and paternal haplotypes, with a zoom panel on the largest
inversion. It is configured under `chr8_synteny:` in `config/config.yaml`.

**Pipeline steps:**
1. `chr8_extract_contig` — extract chr8 from each full-assembly FASTA (samtools faidx)
2. `chr8_index_fasta` — index extracted FASTA and write `.cl` (chromosome-length) file for plotsr
3. `chr8_align` — minimap2 asm5 alignment; wildcards `{ref_samp}_{qry_samp}` (e.g. `ref_mat`, `mat_pat`)
4. `chr8_syri` — SyRI structural rearrangement; same wildcard pattern as align
5. `chr8_find_inversion` — parse `ref_patsyri.out` for largest PAT inversion, write `inversion_coords.json`
6. `chr8_make_figure` — plotsr + matplotlib multi-panel figure

**Required SyRI runs** (plotsr needs consecutive-genome pairs for a 3-genome [REF, MAT, PAT] layout):
- `ref_matsyri.out` — REF↔MAT (first pair)
- `mat_patsyri.out` — MAT↔PAT (second pair; **not** `ref_patsyri.out`)
- `ref_patsyri.out` — REF↔PAT (used only by `chr8_find_inversion` to detect PAT inversions in REF coords)

**Key outputs:**
- `results/chr8_synteny/syri/{ref_samp}_{qry_samp}syri.out`
- `results/chr8_synteny/inversion_coords.json`
- `results/chr8_synteny/chr8_figure.{pdf,png}`

**Convenience target:** `snakemake chr8_synteny`

## CI / GitHub Actions

- **Workflow file:** `.github/workflows/main.yml` — runs on push to `main` and on PRs; also supports `workflow_dispatch` for manual triggering from the Actions UI
- **Jobs:** `r-format-suggest` (air formatter), `py-format` (ruff), `snk-linting` (snakemake --lint), `snk-dry-run` (snakemake --dry-run)
- **snakemake-github-action@v2 defaults:** This action defaults `directory` to `.test` and `snakefile` to `Snakefile`. This repo requires explicit `directory: .` and `snakefile: workflow/Snakefile`.
- **Workflow disabling:** GitHub silently disables workflows after repeated consecutive failures. Re-enable with `gh workflow enable main.yml`. The `workflow_dispatch` trigger allows manual runs as a workaround.

## Common Debugging Patterns

**bcftools annotation error**: "The INFO tag 'CONTEXT_IDS' is not defined"

- **Cause**: Header file was cached with old definitions before rule changes
- **Root cause**: `generate_annotation_headers` rule (inline heredoc in `annotation.smk`) was updated but `results/generate_annotation_headers/` had stale cache
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
- For new genomic context rules: use `get_stratifications_for_ref()` helper function
- Internal helper functions use underscore prefix convention (e.g., `_get_exclusion_config`)
- `BENCHMARKS_WITH_EXCLUSIONS` top-level variable filters benchmarks with exclusions configured

**Chromosome configuration**: Chromosome lists are defined in `config.yaml` under `chromosomes:`, not hardcoded

- `get_chromosomes()` reads from config and applies ref-specific prefixes (GRCh37 vs GRCh38)

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

**SyRI crash**: `ValueError: buffer source array is read-only` in `synsearchFunctions.pyx`

- **Cause**: pandas 2.0 Copy-on-Write makes DataFrame slice arrays non-writeable; SyRI's Cython code
  tries to create a writable memoryview from them
- **Fix**: Pin `pandas<2.0` in `workflow/envs/plotsr.yaml`; delete `.snakemake/conda/<hash>/` to force
  env rebuild, or delete `results/chr8_synteny/syri/` so the rule re-runs in the rebuilt env

**SyRI `--prefix` deprecation warning** (may cause crashes with some SyRI versions):

- **Symptom**: "For specifying output folder use --dir, use --prefix for modifying the output file names"
- **Fix**: Use `--dir <directory> --prefix <basename>` instead of `--prefix <full/path>`
- Already fixed in `chr8_syri` rule (uses `params.outdir` + `params.prefix`)

**SyRI `.out` column indices** (`find_chr8_inversion.py`):

- 0-indexed columns: `[0]`=refChr, `[1]`=refStart, `[2]`=refEnd, `[5]`=**qryChr (string)**, `[6]`=qryStart, `[7]`=qryEnd, `[10]`=type
- Query coordinates are **cols 6–7**, not 5–6. Reading col 5 as int raises `ValueError: int("chr8")`

## Data Loading (R/Quarto)

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
Rscript -e 'testthat::test_file("tests/test_cache.R")'           # Schema + cache tests (45 tests)
Rscript -e 'testthat::test_file("tests/test_data_loading.R")'    # Data loading tests (has pre-existing failures)
Rscript -e 'testthat::test_file("tests/test_exclusion_loading.R")'  # Exclusion loading tests
Rscript -e 'testthat::test_file("tests/test_schema_update.R")'   # Schema update tests
pytest tests/unit/  # Python unit tests (common helpers)
```

Cache tests use `withr::local_options(q100.cache_dir = tempdir)` for isolation.

### Known Issues

- `parse_benchmark_id()` regex captures `"5.0q"` not `"v5.0q"` -- tests expect `"v"` prefix
- `test_data_loading.R` references old column names (`strat_name`, `strat_filter`) from pre-refactor
- Tests that load real data fail without `results/` directory (pipeline outputs not in repo)

## Code Quality

- **Formatter:** `air` (config in `air.toml`: 100-char lines, 2-space indent)
- **Linter:** `lintr` (config in `.lintr`: allows `dotted.case` for internal helpers; excludes `R/data_loading.R` line 7 for `%>%` assignment false positive)
- **lintr/styler in renv:** Both are recorded in `renv.lock` and available when renv is active
- **Pipe Operator:** Use `%>%` (magrittr pipe) instead of `|>` (native pipe) for consistency with existing codebase

### renv Package Management

- **Adding dev-only packages:** Use `renv::record("pkg")` to add packages to `renv.lock` without disturbing existing entries. **Do NOT use `renv::snapshot()`** — its implicit mode strips packages not imported in project code, which will remove most of the lockfile (lintr, styler, and other dev tools are not imported).
- **renv vs conda:** `.Rprofile` activates renv which isolates R library paths. Packages installed via conda are not visible when renv is active. Dev tools (lintr, styler) must be added to renv directly.

## Plot Themes & Styling

The `R/plot_themes.R` module provides consistent ggplot2 and gt table styling for all manuscript figures and tables (Cell Genomics submission requirements). Key functions support flexible parameter overrides via `...`:

### Theme Functions

**`theme_manuscript(...)`**
- Base manuscript theme with journal-appropriate sizing (8-11pt fonts, 85-180mm figure width)
- Accepts `...` to override any theme element
- Example: `theme_manuscript(axis.title = element_text(size = 12, face = "italic"))`
- Returns composition: `theme_minimal() + base_theme + theme(...)`

**`theme_gt_manuscript(gt_object, striped = TRUE, ...)`**
- gt table styling with consistent fonts, borders, and spacing
- Accepts `...` to override any gt table option via `utils::modifyList()`
- Example: `theme_gt_manuscript(gt_table, table.font.size = "10pt")`
- User options take precedence over defaults

### Scale Functions

All scale functions accept `...` to pass additional parameters to underlying ggplot2 functions:

**`scale_benchmark_version(aesthetic, name, guide, ...)`**
- Colors for benchmark versions (v0.6, v4.2.1, v5.0q, PP)
- Example: `scale_benchmark_version(fill = "color", limits = c("v5.0q"))`

**`scale_bench_type(aesthetic, name, guide, ...)`**
- Colors for small variants (smvar) vs structural variants (stvar)
- Auto-applies readable labels: "Small Variants", "Structural Variants"

**`scale_genomic_context(aesthetic, name, guide, ...)`**
- Colors for genomic context types (HP, MAP, SD, SD10kb, TR, TR10kb)
- Auto-applies readable labels: "Homopolymers", "Low Mappability", etc.

### Color Palettes

`get_color_palettes()` returns all available palettes:
- `bench_version` — distinct colors for benchmark versions
- `ref` — colors for reference genomes (GRCh37, GRCh38, CHM13v2.0)
- `bench_type` — colors for smvar/stvar
- `context_name` — colors for genomic contexts
- `chrom_type` — colors for autosomes vs sex chromosomes
- `binary` — colors for TRUE/FALSE and Yes/No

All palettes are colorblind-friendly and print-friendly.

## Important Patterns

- `fs::dir_ls` glob matches against **full path** via `glob2rx` -- use `grepl()` on `fs::path_file()` for filename-only matching
- Factor levels are centralized in `R/schemas.R` constants (e.g., `BENCH_VERSION_LEVELS`, `CHROM_LEVELS`)
- All loading functions have `use_cache = TRUE` and `force_refresh = FALSE` defaults; cache writes are wrapped in `tryCatch` so failures don't break loading

**Quarto notebooks:**

- `analysis/benchmarkset_characterization.qmd` — Primary analysis; loads `variants_df` and `genomic_context_metrics_df`
- `analysis/benchmark_difficult.qmd` — Coverage analysis; loads `diff_cov_df`
- `analysis/benchmark_exclusions.qmd` — Exclusion region analysis (v5.0q only)
- `analysis/genomic_context_analysis.qmd` — Genomic context overlap analysis
- `analysis/benchmark_interval_size_distributions.qmd` — Interval size distributions
- `analysis/benchmark_unique_regions.qmd` — Unique region analysis across benchmark versions
- `analysis/external_evaluation.qmd` — External benchmark comparisons
- `analysis/_notebook_setup.R` — `analysis_setup()` helper: loads tidyverse/patchwork, sources `R/data_loading.R` and `R/plot_themes.R`; call at top of each notebook

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
results/genomic_context/v5.0q_GRCh38_smvar/variants_by_genomic_context.parquet
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
