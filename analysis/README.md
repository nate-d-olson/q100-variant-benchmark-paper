# Analysis Notebooks

This directory contains Quarto notebooks for analyzing the GIAB v5q benchmark set data.

## Notebooks

### benchmarkset-characterization.qmd

Comprehensive characterization of the v5q benchmark set variants and regions, including comparisons to historical benchmark sets.

**Input Files:**

| File | Source | Description |
|------|--------|-------------|
| `results/variant_tables/{benchmark}/variants.tsv` | Snakemake pipeline | Annotated variant tables for each benchmark set |
| `resources/benchmarksets/{benchmark}_benchmark.bed` | Snakemake pipeline | Benchmark region BED files |
| `data/platinum-pedigree-data/NA12878_hq_v1.2.smallvar.bed.gz` | External download | Platinum Pedigree small variant regions |
| `data/platinum-pedigree-data/NA12878_hq_v1.2.svs.bed.gz` | External download | Platinum Pedigree SV regions |
| `results/context_coverage/all_context_coverage.tsv` | Snakemake pipeline | Context coverage summary |
| `data/20250117_v0.020_HG002Q100v1.1/draft_benchmarksets/*/exclusion_intersection_summary.csv` | DeFrABB run | Exclusion intersection summaries |
| `config/config.yaml` | Repository | Pipeline configuration with genome sizes |

**Output Files:**

| File | Description |
|------|-------------|
| `manuscript/figs/sv_length_distribution.png` | SV length distribution histogram |
| `manuscript/figs/diff_variant_counts.png` | Variant counts in difficult genomic contexts |
| `manuscript/figs/benchmark_region_size_change.png` | Coverage change between benchmark versions |
| `manuscript/figs/benchmark_intervals.png` | Benchmark interval size distributions |

**Key Analyses:**

- Variant counts by type (SNP, INDEL, INS, DEL) across benchmark versions
- Size distributions for indels and structural variants
- Variants in difficult genomic contexts (tandem repeats, homopolymers, segmental duplications, low mappability)
- Benchmark region coverage comparisons
- Inclusion of difficult regions in benchmark
- Interval size distributions
- Exclusion region analysis

---

### external_evaluation.qmd

Analysis of external evaluations from variant calling pipelines compared against the v5q benchmark.

**Input Files:**

| File | Source | Description |
|------|--------|-------------|
| `data/external-evaluations/evaluator-curations/*Miqa.csv` | External evaluators | Curation results from evaluators |

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

### Prerequisites

Ensure the Snakemake pipeline has been run to generate input files:

```bash
# From project root
snakemake --cores 4 --sdm conda
```

### Rendering Notebooks

```bash
# Render a single notebook
quarto render analysis/benchmarkset-characterization.qmd

# Render all analysis notebooks
quarto render analysis/
```

### Required R Packages

The notebooks use the following R packages:

- `tidyverse` - Data manipulation and visualization
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
install.packages(c("tidyverse", "DT", "here", "assertthat",
                   "patchwork", "gt", "vroom", "furrr", "yaml"))
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

- **Platinum Pedigree**: Download from AWS S3
  ```bash
  aws s3 sync --no-sign-request s3://platinum-pedigree-data/truthset_v1.2/ data/platinum-pedigree-data/
  ```

- **External Evaluations**: Curation files from evaluators in `data/external-evaluations/`

## Output Figures

Generated figures are saved to `manuscript/figs/` for inclusion in the manuscript.
