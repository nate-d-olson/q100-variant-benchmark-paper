# Pipeline Architecture

Snakemake-based workflow for analyzing the GIAB Q100 HG002 variant benchmark across
GRCh37, GRCh38, CHM13v2.0 and benchmark versions v0.6, v4.2.1, v5.0q. The pipeline
produces aggregated metrics consumed by Quarto notebooks in `analysis/` and the
`manuscript/` chapters.

## Top-level layout

```
config/                   # config.yaml + JSON schema
workflow/
├── Snakefile             # min_version("8.0"); entrypoint
├── rules/                # 9 rule files (see below)
├── scripts/              # 10 Python scripts (see below)
└── envs/                 # 6 conda environments
R/                        # data loading + plot themes (Quarto consumers)
analysis/                 # 8 Quarto notebooks
manuscript/               # Quarto manuscript chapters
results/                  # gitignored pipeline outputs
resources/                # gitignored downloaded benchmarks/refs/strats
```

## Rule files (`workflow/rules/`)

| File | Lines | Purpose |
|---|---|---|
| `common.smk` | 274 | Helper functions, wildcard constraints, BENCHMARKS_WITH_EXCLUSIONS |
| `downloads.smk` | 333 | SHA256-verified downloads (benchmarks, refs, stratifications, exclusions) |
| `ref_processing.smk` | 31 | Reference indexing (`samtools faidx`, `seqkit stats`) |
| `vcf_processing.smk` | 79 | VCF normalization and indexing |
| `annotation.smk` | 173 | Combine BEDs + annotate VCFs with `INFO/CONTEXT_IDS`, `INFO/REGION_IDS` |
| `genomic_context_analysis.smk` | 129 | bedtools coverage, coverage tables, variant Parquet, per-context counts |
| `exclusions.smk` | 258 | v5.0q exclusion impact, interactions, cross-version analysis |
| `benchmark_comparisons.smk` | 123 | Truvari comparison between benchmark versions |
| `chr8_synteny.smk` | 247 | Chr8 multi-panel SyRI/plotsr figure pipeline |
| `use_case_evaluation.smk` | 53 | SV use-case metric extraction (wraps `scripts/extract_sv_metrics.py`); opt-in target `snakemake use_case_evaluation` |

`common.smk` reads `config.yaml`, builds `BENCHMARKS_WITH_EXCLUSIONS` from
benchmarks that have exclusions configured, and exposes `get_*` helpers used by
all rule files.

## Python scripts (`workflow/scripts/`)

| Script | Used by | Purpose |
|---|---|---|
| `logging_config.py` | all | Centralized logger setup |
| `combine_beds_with_id.py` | annotation.smk | Merge BEDs and inject the ID column for `INFO/*_IDS` annotations |
| `compute_coverage_table.py` | genomic_context_analysis.smk | Per-context BED overlap → CSV |
| `generate_variant_parquet.py` | genomic_context_analysis.smk | Annotated VCF → Parquet via Truvari `VariantRecord` API; includes `normalize_annotation()` for tuple-to-string conversion |
| `count_variants_by_genomic_context.py` | genomic_context_analysis.smk | Group variants by `(context_name, var_type, szbin)` |
| `count_exclusion_variants.py` | exclusions.smk | Variant counts per exclusion |
| `compute_exclusion_interactions.py` | exclusions.smk | Upset-style decomposition of exclusion overlaps |
| `annotate_old_benchmark_status.py` | exclusions.smk | Cross-version v4.2.1/v0.6 vs v5.0q annotation |
| `find_chr8_inversion.py` | chr8_synteny.smk | Parse SyRI `.out` for largest PAT inversion |
| `make_chr8_figure.py` | chr8_synteny.smk | Multi-panel matplotlib + plotsr figure |

## Conda environments (`workflow/envs/`)

Consolidated from 8 → 6 envs on 2026-02-23. See `workflow/envs/README.md` for
rationale.

| Env | Purpose | Notable pins |
|---|---|---|
| `biotools.yaml` | Core CLI (bcftools, rtg-tools) | — |
| `python-biotools.yaml` | Python data processing + bcftools/bedtools/tabix | python=3.11 |
| `samtools.yaml` | Sequence handling (samtools, seqkit) | — |
| `downloads.yaml` | wget for downloads | — |
| `plotsr.yaml` | Chr8 synteny (minimap2, syri, plotsr) | **pandas<2.0** (SyRI Cython bug) |
| `truvari.yaml` | Variant analysis | Truvari==5.4.0 (pip), bcftools=1.20 |

## Data Flow

```
config.yaml
    │
    ▼
Downloads (SHA256-verified)
    │  benchmark VCFs/BEDs, references, stratifications, exclusions
    ▼
VCF normalization + indexing            Reference indexing
    │                                       │
    └──────────────┬────────────────────────┘
                   ▼
         Combine BEDs with ID columns
         (combine_beds_with_id.py)
                   │
                   ▼
         Annotate VCFs (bcftools annotate)
         INFO/CONTEXT_IDS, INFO/REGION_IDS
                   │
        ┌──────────┼──────────┬───────────────┐
        ▼          ▼          ▼               ▼
  bedtools     Truvari     Exclusion     Truvari benchmark
  coverage     Parquet     impact        comparisons
  per-context  generation  + interactions
        │          │          │               │
        └──────────┼──────────┘               │
                   ▼                          │
         Per-context counts                   │
         (count_variants_by_genomic_context)  │
                   │                          │
                   ▼                          ▼
         results/genomic_context/   results/comparisons/
         results/variant_tables/
         results/exclusions/
                   │
                   ▼
         R loaders (R/data_loading.R) → Parquet cache (analysis/cache/)
                   │
                   ▼
         Quarto notebooks → manuscript figures/tables
```

## Annotation Model

Variants are annotated by `bcftools annotate` with two multi-value INFO fields:

| Field | Source BED | Loader column |
|---|---|---|
| `CONTEXT_IDS` | `combine_genomic_context_beds` (HP, MAP, SD, …) | `context_ids` |
| `REGION_IDS` | `combine_region_beds` (benchmark + per-exclusion BEDs) | `region_ids` |

Truvari's `vcf_to_df()` converts these to Python tuples; `normalize_annotation()`
in `generate_variant_parquet.py` converts them back to comma-separated strings
before Parquet write (fix shipped in commit `78747aa`).

## R Data Layer

```
R/
├── schemas.R       # Arrow schema, factor levels, validation rules
├── cache.R         # Parquet caching infrastructure
├── data_loading.R  # 11 load_* functions (sources schemas + cache)
├── plot_themes.R   # ggplot2 themes + flextable helpers (gt legacy)
└── bed_helpers.R   # Interval arithmetic for benchmark_unique_regions notebook
```

### Schemas (`schemas.R`)

`get_arrow_schema(dataset)` returns Arrow schema for: `variant_table`,
`diff_coverage`, `benchmark_regions`, `platinum_pedigree_regions`. Factor levels
in `get_factor_levels()` restore factors after Parquet read (Arrow stores them as
strings).

### Cache (`cache.R`)

Parquet + zstd in `analysis/cache/` (override via `getOption("q100.cache_dir")`).
Cache key: `rlang::hash()` of dataset name + source-file mtimes + parameters.
Pipeline metadata (R version, package versions, config summary) embedded as
JSON in Parquet file-level kv-metadata under `q100_pipeline_metadata`.

Cached loaders: `load_variant_table()`, `load_genomic_context_coverage()`,
`load_benchmark_regions()`, `load_platinum_pedigree_regions()`. Small-dataset
loaders (metrics, exclusions, reference sizes) read directly each call.

## Notebooks (`analysis/`)

| Notebook | Loads |
|---|---|
| `benchmarkset_characterization.qmd` | `variants_df`, `genomic_context_metrics_df` |
| `benchmark_difficult.qmd` | `diff_cov_df` |
| `benchmark_exclusions.qmd` | `exclusion_metrics_df`, `exclusion_interactions_df` |
| `benchmark_unique_regions.qmd` | cross-version `old_only_*.csv`, BED helpers |
| `benchmark_interval_size_distributions.qmd` | `bench_regions_df` |
| `genomic_context_analysis.qmd` | per-context overlaps |
| `external_evaluation.qmd` | external benchmark comparisons |
| `use_case_evaluation.qmd` | hap.py + Truvari outputs (manual; **not yet wired into Snakemake** — see TODO.md) |

`analysis/_notebook_setup.R` provides `analysis_setup()` — loads tidyverse,
sources `R/data_loading.R` and `R/plot_themes.R`. Call at the top of each notebook.

## Tooling

- **Workflow:** Snakemake 8.x (`min_version("8.0")` in Snakefile); execution via
  `make run` (`--cores 14 --sdm conda`)
- **Linting:** `make lint` runs `snakemake --lint`, `ruff check`, `lintr`
- **Formatting:** `make format` runs `ruff format`, `snakefmt`, `prettier`, `styler`
- **Pre-merge gate:** `make test` (lint + format-check + dry-run)
- **CI:** `.github/workflows/main.yml` runs format checks, snakemake lint, dry-run.
  No unit tests run in CI yet — `pytest` and `Rscript testthat` are manual.

## Adding a New Cached Dataset

1. `R/schemas.R` — add to `get_arrow_schema()`, `get_factor_levels()`,
   `get_validation_rules()`
2. `R/data_loading.R` — add `use_cache` and `force_refresh` parameters; call
   `read_cache()` before loading and `write_cache()` after

`R/cache.R` is dataset-agnostic — no changes needed.

## Related Docs

- [Pipeline Outputs Reference](pipeline-outputs.md)
- [Data Dictionary](data-dictionary.md)
- [Troubleshooting](troubleshooting.md)
- [Plot Themes Guide](plot_themes_guide.md)
- `workflow/envs/README.md` — env consolidation rationale
- `CLAUDE.md` — project-specific conventions and known issues
