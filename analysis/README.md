# Analysis Notebooks

This directory contains Quarto notebooks for downstream analysis and figure generation from workflow outputs.

## Current Notebooks

- `benchmarkset_characterization.qmd`
- `benchmark_unique_regions.qmd`
- `benchmark_interval_size_distributions.qmd`
- `benchmark_exclusions.qmd`
- `benchmark_difficult.qmd`
- `benchmark_difficult_var.qmd`
- `external_evaluation.qmd`
- `genomic_context_analysis.qmd`

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

- `manuscript/figs/`
- `manuscript/tables/`

Examples currently produced by notebooks include:

- `manuscript/figs/benchmark_intervals.png`
- `manuscript/figs/benchmark_region_size_change.png`
- `manuscript/figs/genomic_context_variant_counts.png`
- `manuscript/figs/genomic_context_variant_counts_tall.png`
- `manuscript/figs/smvar_strat_flow.png`
- `manuscript/figs/stvar_eval_strata_supplemental.png`
- `manuscript/figs/combined_eval.png`

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
analysis_setup(load_plot_themes = TRUE, load_gt = TRUE)
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
