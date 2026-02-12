# Analysis Notebooks

This directory contains Quarto notebooks for analyzing the GIAB v5q benchmark set data.

## Notebooks

### benchmarkset_characterization.qmd

Comprehensive characterization of the v5q benchmark set variants and regions, including comparisons to historical benchmark sets.

**Input Files:**

| File | Source | Description |
|------|--------|-------------|
| `results/variant_tables/{benchmark}/variants.tsv` | Snakemake pipeline | Annotated variant tables for each benchmark set |
| `resources/benchmarksets/{benchmark}_benchmark.bed` | Snakemake pipeline | Benchmark region BED files |
| `results/context_coverage/all_context_coverage.tsv` | Snakemake pipeline | Context coverage summary |
| `data/20250117_v0.020_HG002Q100v1.1/draft_benchmarksets/*/exclusion_intersection_summary.csv` | DeFrABB run | Exclusion intersection summaries |
| `config/config.yaml` | Repository | Pipeline configuration with genome sizes |

**Output Files:**

| File | Description |
|------|-------------|
| `manuscript/figs/sv_length_distribution.png` | SV length distribution histogram |
| `manuscript/figs/diff_variant_counts.png` | Variant counts in difficult genomic contexts |
| `manuscript/figs/benchmark_region_size_change.png` | Coverage change between benchmark versions |

**Key Analyses:**

- Variant counts by type (SNV, INDEL, INS, DEL) across benchmark versions
- Size distributions for indels and structural variants
- Variants in difficult genomic contexts (tandem repeats, homopolymers, segmental duplications, low mappability)
- Benchmark region coverage comparisons
- Inclusion of difficult regions in benchmark
- Exclusion region analysis

---

### benchmark_unique_regions.qmd

Focused comparison of benchmark regions and variants unique to v5 or previous benchmark versions for:

- small variants on `GRCh38` (`v5.0q` vs `v4.2.1`)
- structural variants on `GRCh37` (`v5.0q` vs `v0.6`)

**Input Files:**

| File | Source | Description |
|------|--------|-------------|
| `resources/benchmarksets/{benchmark}_benchmark.bed` | Snakemake pipeline | Benchmark region BED files |
| `resources/benchmarksets/{benchmark}_dip.bed` | Snakemake pipeline | v5 diploid regions |
| `resources/exclusions/{benchmark}/*.bed` | Snakemake pipeline | v5 exclusion BED files |
| `results/exclusions/{comp_id}/old_only_variants.tsv` | Snakemake pipeline | Old-only variant status with exclusion overlap |
| `results/variant_tables/{benchmark}/variants.tsv` | Snakemake pipeline | Variant tables for v5-only region variant summaries |
| `resources/stratifications/{ref}_{context}.bed.gz` | Snakemake pipeline | Genomic context BED files |

**Key Analyses:**

- unique base and variant counts for v5-only vs previous-only benchmark regions
- previous-only region status (`not_in_dipbed`, `excluded`, `in_v5_dipbed_not_excluded`)
- exclusion-overlap reasons and genomic context for excluded previous-only regions
- v5-only variant type/size breakdown and genomic context for new regions

---

### benchmark_interval_size_distributions.qmd

Benchmark interval size distribution analysis, combining benchmark BED regions with Platinum Pedigree regions.

**Input Files:**

| File | Source | Description |
|------|--------|-------------|
| `resources/benchmarksets/{benchmark}_benchmark.bed` | Snakemake pipeline | Benchmark region BED files |
| `resources/platinum-pedigree-data/truthset_v1.2/NA12878_hq_v1.2.smallvar.bed.gz` | Public S3 via `load_platinum_pedigree_regions()` | Platinum Pedigree small variant regions |
| `resources/platinum-pedigree-data/truthset_v1.2/NA12878_hq_v1.2.svs.bed.gz` | Public S3 via `load_platinum_pedigree_regions()` | Platinum Pedigree SV regions |

**Output Files:**

| File | Description |
|------|-------------|
| `manuscript/figs/benchmark_intervals.png` | Benchmark interval size distributions |

**Key Analyses:**

- Interval size density distributions across benchmark sets
- Interval count comparisons by benchmark version
- Combined benchmark + Platinum Pedigree interval summary table

---

### external_evaluation.qmd

Analysis of external evaluations from variant calling pipelines compared against the v5q benchmark.

**Input Files:**

| File | Source | Description |
|------|--------|-------------|
| `data/external-evaluations/evaluator-curations/*Miqa.csv` | External evaluators | Curation results from evaluators |
| NA | `https://docs.google.com/spreadsheets/d/1evd6q8jqF51IM_uX6XrTlc-3f7dDUvurlGI-\_vcV6t0/edit?usp=sharing` | JZ re-curations and notes for external evaluations |

**Output Files:**

| File | Description |
|------|-------------|
| `manuscript/figs/smvar_eval_strata_supplemental.png` | Small variant evaluations by strata |
| `manuscript/figs/stvar_eval_strata_supplemental.png` | Structural variant evaluations by strata |
| `manuscript/figs/combined_eval_callset.png` | Combined evaluation summary by callset |

**Key Analyses:**

- Small variant evaluation summaries by callset and strata
- Structural variant evaluation summaries by callset and strata
- Identification of variants marked as "not correct" or "unsure" in benchmark
- Visualization of evaluation results across different genomic contexts

---

## Running the Notebooks

### Shared Notebook Setup

Analysis notebooks now use `analysis/_notebook_setup.R` to reduce repeated setup code.
Use this helper in notebook setup chunks to consistently load packages and shared
R modules:

```r
source(here::here("analysis/_notebook_setup.R"))
analysis_setup(load_plot_themes = TRUE, load_gt = TRUE)
```

The data-loading layer also includes `load_primary_analysis_data()` to load
validated core data frames in one call.

### Prerequisites

Ensure the Snakemake pipeline has been run to generate input files:

```bash
# From project root
snakemake --cores 4 --sdm conda
```

### Rendering Notebooks

```bash
# Render a single notebook
quarto render analysis/benchmarkset_characterization.qmd

# Render unique-region comparison notebook
quarto render analysis/benchmark_unique_regions.qmd

# Render interval size distribution notebook
quarto render analysis/benchmark_interval_size_distributions.qmd

# Render all analysis notebooks
quarto render analysis/
```

### Required R Packages

The notebooks use the following R packages:

- `tidyverse` - Data manipulation and visualization
- `arrow` - Parquet caching and Arrow schema support
- `jsonlite` - Pipeline metadata serialization
- `DT` - Interactive data tables
- `here` - Project-relative file paths
- `assertthat` - Data validation assertions
- `patchwork` - Plot composition
- `gt` - Publication-quality tables
- `vroom` - Fast file reading
- `furrr` - Parallel processing
- `yaml` - YAML file parsing

Install packages via:

```r
install.packages(c("tidyverse", "arrow", "jsonlite", "DT", "here",
                   "assertthat", "patchwork", "gt", "vroom", "furrr", "yaml"))
```

### Parquet Data Cache

Large datasets (variant tables, coverage files, benchmark regions) are cached as Parquet files in `analysis/cache/` for fast reloading. The cache is:

- Automatically populated on first load via `R/data_loading.R` functions
- Invalidated when source files change (based on modification times)
- Validated during data loading before cache generation, with schema checks
  re-applied at cache write time for defense in depth
- Excluded from git (see `.gitignore`)

Cache management:

```r
source(here::here("R/data_loading.R"))

cache_info()                        # List cached files
invalidate_cache("variant_table")   # Clear one dataset
clear_cache()                       # Clear everything
```

## Data Dependencies

### From Snakemake Pipeline

The analysis notebooks depend on outputs from the Snakemake pipeline:

```
results/
├── variant_tables/           # Per-benchmark variant tables
│   ├── v5q_chm13_smvar/variants.tsv
│   ├── v5q_grch38_smvar/variants.tsv
│   ├── v421_grch38_smvar/variants.tsv
│   └── ...
├── context_coverage/         # Context coverage summaries
│   └── all_context_coverage.tsv
└── exclusions/               # Exclusion analysis
    └── {benchmark}/exclusions_intersection_table.csv

resources/
└── benchmarksets/            # Benchmark BED files
    ├── v5q_chm13_smvar_benchmark.bed
    └── ...
```

### External Data

Some notebooks require data downloaded separately:

- **Platinum Pedigree**: loaded via `load_platinum_pedigree_regions()`, which
  downloads from public S3 into `resources/platinum-pedigree-data/truthset_v1.2/`
  when required files are missing.
- **External Evaluations**: Curation files from evaluators in `data/external-evaluations/`

## Output Figures

Generated figures are saved to `manuscript/figs/` for inclusion in the manuscript.
