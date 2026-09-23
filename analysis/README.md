# Analysis Notebooks

This directory contains Quarto notebooks for downstream analysis and figure generation from workflow outputs.

## Current Notebooks

Manuscript figures and tables (full map in the top-level `README.md`):

- `benchmarkset_characterization.qmd`: Table 1, Fig 2
- `external_evaluation.qmd`: Fig 5
- `benchmark_interval_size_distributions.qmd`: Fig 6
- `benchmark_exclusions.qmd`: Table 4
- `manuscript_value_verification.qmd`: checks values quoted in the manuscript text

Supporting analyses (no manuscript figure or table):

- `benchmark_difficult.qmd`
- `benchmark_unique_regions.qmd`
- `genomic_context_analysis.qmd`
- `use_case_evaluation.qmd`

Shared setup helper:

- `_notebook_setup.R`

## Typical Inputs

Notebooks primarily consume:

- `results/variant_tables/*/variants.parquet`
- `results/genomic_context/*/genomic_context_coverage_table.csv`
- `results/genomic_context/*/variants_by_genomic_context.parquet`
- `results/exclusions/**`
- `resources/benchmarksets/**`
- `resources/stratifications/**`
- `config/config.yaml`
- `data/external-evaluations/*` (for `external_evaluation.qmd`)

## Typical Outputs

Most rendered figures/tables are written to:

- `figures/`

Manuscript figures produced by notebooks:

- `figures/variant_size_genomic_context.{pdf,png}` (Fig 2)
- `figures/combined_eval_strata.{pdf,png}` (Fig 5)
- `figures/benchmark_intervals.{pdf,png}` (Fig 6)

Tables 1 and 4 are rendered inline in the notebook HTML.

## Running

From repository root:

```bash
# Render one notebook
quarto render analysis/benchmarkset_characterization.qmd

# Render all notebooks in this directory
quarto render analysis/
```

For consistent setup in notebook chunks:

```r
source(here::here("analysis/_notebook_setup.R"))
analysis_setup(load_plot_themes = TRUE)
```

## Cache

Local cache artifacts are stored in `analysis/cache/`.

- Populated by data-loading helpers in `R/data_loading.R`
- Intended as local acceleration artifacts
- Excluded from git

Useful helpers:

```r
source(here::here("R/data_loading.R"))
cache_info()
invalidate_cache("variant_table")
clear_cache()
```

## Rendered HTML Artifacts

Rendered notebook HTML files in this folder are local build artifacts and are gitignored.
