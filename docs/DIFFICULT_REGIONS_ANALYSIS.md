# Difficult Regions Analysis Guide

This document explains how to analyze variant counts and benchmark region coverage for difficult genomic regions (stratifications) using the Q100 variant benchmark pipeline.

## Overview

The pipeline provides two complementary analyses for difficult regions:

1. **Stratification Coverage Metrics** - Quantifies how much of each difficult region is covered by benchmark regions
2. **Variant Counts by Stratification** - Counts benchmark variants within each difficult region

These analyses enable:
- Understanding benchmark representativeness across difficult genomic contexts
- Generating cumulative plots showing progressive coverage and variant inclusion
- Comparing benchmarks on difficult region characteristics

## Pipeline Outputs

### 1. Stratification Coverage Tables

**Location:** `results/strat_metrics/{benchmark}/stratification_coverage_table.csv`

**Columns:**
- `strat_name`: Stratification region name (TR, HP, SD, MAP, etc.)
- `strat_bp`: Total bases in the stratification region (after merging overlaps)
- `intersect_bp`: Bases overlapping between stratification and benchmark regions
- `pct_of_strat`: Percentage of stratification covered by benchmark
- `pct_of_dip`: Percentage of benchmark regions covered by stratification

**Example:**
```csv
strat_name,strat_bp,intersect_bp,pct_of_strat,pct_of_dip
TR,245691234,123456789,50.25,4.12
HP,89234567,45123456,50.58,1.51
SD,134567890,67890123,50.45,2.27
```

### 2. Variant Count Tables

**Location:** `results/var_counts/{benchmark}/variants_by_stratification.csv`

**Columns:**
- `strat_name`: Stratification region name
- `var_type`: Variant type (SNP, INDEL, DEL, INS, COMPLEX)
- `count`: Number of variants

**Example:**
```csv
strat_name,var_type,count
TR,SNP,45678
TR,INDEL,12345
HP,SNP,23456
HP,INDEL,6789
```

### 3. Summary Tables

**Location:** `results/var_counts/{benchmark}/stratification_summary.csv`

**Columns:**
- `strat_name`: Stratification region name
- `total_variants`: Total variant count across all types
- `snp_count`, `indel_count`, `del_count`, `ins_count`, `complex_count`, `other_count`: Type-specific counts

### 4. Combined Metrics

**Location:** `results/var_counts/{benchmark}/stratification_combined_metrics.csv`

**Columns:**
- All columns from stratification coverage table
- All columns from summary table
- `variant_density_per_mb`: Variants per megabase of benchmark-covered stratification

**This is the primary file for generating cumulative plots and comprehensive analyses.**

## Running the Pipeline

### Generate All Outputs

```bash
# From the project root directory
snakemake --cores 8
```

This will generate:
- Variant tables (existing)
- Exclusion tables (existing)
- **Stratification coverage metrics (NEW)**
- **Variant counts by stratification (NEW)**
- Reference genome sizes (existing)

### Generate Specific Outputs

```bash
# Only stratification metrics
snakemake --cores 8 results/strat_metrics/v5q_grch38_smvar/stratification_coverage_table.csv

# Only variant counts
snakemake --cores 8 results/var_counts/v5q_grch38_smvar/stratification_combined_metrics.csv

# All metrics for one benchmark
snakemake --cores 8 \
  results/strat_metrics/v5q_grch38_smvar/stratification_coverage_table.csv \
  results/var_counts/v5q_grch38_smvar/stratification_combined_metrics.csv
```

## Generating Cumulative Plots

### Using the Provided R Script

```bash
# Plot for a specific benchmark
Rscript scripts/plot_difficult_region_cumulative.R v5q_grch38_smvar

# Default (v5q_grch38_smvar)
Rscript scripts/plot_difficult_region_cumulative.R
```

**Outputs:**
- `results/plots/{benchmark}/difficult_regions_cumulative.png` - Combined plot
- `results/plots/{benchmark}/difficult_regions_summary.csv` - Summary table

### Custom Analysis in R

```r
library(tidyverse)
library(here)

# Load combined metrics
metrics <- read_csv(here(
  "results/var_counts/v5q_grch38_smvar/stratification_combined_metrics.csv"
))

# Sort by benchmark coverage percentage
coverage_sorted <- metrics %>%
  arrange(desc(pct_of_dip)) %>%
  mutate(
    cumulative_pct = cumsum(pct_of_dip),
    cumulative_bp = cumsum(intersect_bp)
  )

# Cumulative coverage plot
ggplot(coverage_sorted, aes(x = seq_along(strat_name))) +
  geom_step(aes(y = cumulative_pct), linewidth = 1.2) +
  geom_point(aes(y = cumulative_pct), size = 3) +
  scale_x_continuous(breaks = seq_along(coverage_sorted$strat_name),
                    labels = coverage_sorted$strat_name) +
  labs(
    title = "Cumulative Benchmark Coverage by Difficult Regions",
    x = "Stratification (ordered by coverage %)",
    y = "Cumulative Coverage (%)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

### Custom Analysis in Quarto

Add to your Quarto document (e.g., `analysis/difficult_regions.qmd`):

```{r}
# Load all benchmarks
benchmark_files <- list.files(
  path = here("results/var_counts"),
  pattern = "stratification_combined_metrics.csv",
  recursive = TRUE,
  full.names = TRUE
)

combined_metrics <- benchmark_files %>%
  set_names(str_extract(., "(?<=var_counts/).*(?=/strat)")) %>%
  map_dfr(read_csv, show_col_types = FALSE, .id = "benchmark")

# Compare coverage across benchmarks
combined_metrics %>%
  separate(benchmark, into = c("version", "ref", "type"), sep = "_") %>%
  ggplot(aes(x = strat_name, y = pct_of_dip, fill = version)) +
  geom_col(position = "dodge") +
  facet_wrap(~ type, scales = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

## Understanding the Metrics

### Stratification Coverage Interpretation

**`pct_of_strat`** (% of stratification in benchmark):
- High value (>50%): Benchmark well-represents this difficult region
- Low value (<20%): Benchmark may under-represent this difficult region
- Use this to assess **representativeness**

**`pct_of_dip`** (% of benchmark in stratification):
- High value (>10%): Large portion of benchmark is in this difficult region
- Low value (<1%): Small portion of benchmark is in this difficult region
- Use this to assess **difficulty** of benchmark

### Variant Density Interpretation

**`variant_density_per_mb`**:
- Calculated as: `total_variants / (intersect_bp / 1,000,000)`
- Higher density → More variants per unit of benchmark-covered difficult region
- Useful for comparing variant richness across stratifications

### Example: Tandem Repeats (TR)

```
strat_name: TR
strat_bp: 245,691,234
intersect_bp: 123,456,789
pct_of_strat: 50.25
pct_of_dip: 4.12
total_variants: 45,678
variant_density_per_mb: 369.8
```

**Interpretation:**
- 50% of all tandem repeat regions are covered by the benchmark (good representativeness)
- 4% of the benchmark consists of tandem repeats (moderate difficulty)
- 370 variants per Mb in TR regions (high variant density)

## Comparing Benchmarks

### Example: v5.0q vs v4.2.1 Small Variants on GRCh38

```r
# Load both benchmarks
v5q <- read_csv("results/var_counts/v5q_grch38_smvar/stratification_combined_metrics.csv")
v421 <- read_csv("results/var_counts/v421_grch38_smvar/stratification_combined_metrics.csv")

# Compare coverage
comparison <- v5q %>%
  select(strat_name, pct_of_dip_v5q = pct_of_dip, vars_v5q = total_variants) %>%
  left_join(
    v421 %>% select(strat_name, pct_of_dip_v421 = pct_of_dip, vars_v421 = total_variants),
    by = "strat_name"
  ) %>%
  mutate(
    coverage_change = pct_of_dip_v5q - pct_of_dip_v421,
    variant_change = vars_v5q - vars_v421
  )

# Plot changes
ggplot(comparison, aes(x = coverage_change, y = variant_change, label = strat_name)) +
  geom_point(size = 3) +
  geom_text_repel() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    title = "Difficult Region Coverage and Variant Changes: v5.0q vs v4.2.1",
    x = "Change in Benchmark Coverage (%)",
    y = "Change in Variant Count"
  )
```

## Integration with Existing Analysis

The difficult regions analysis integrates seamlessly with the existing Quarto documents in `analysis/`:

### Update `benchmarkset-characterization.qmd`

Replace the TODO section (lines 762-896) with:

```{r}
# Load stratification metrics
strat_metrics_files <- list.files(
  path = here("results/var_counts"),
  pattern = "stratification_combined_metrics.csv",
  recursive = TRUE,
  full.names = TRUE
) %>%
  set_names(str_extract(., "(?<=var_counts/).*(?=/strat)"))

strat_metrics_df <- strat_metrics_files %>%
  map_dfr(read_csv, show_col_types = FALSE, .id = "benchmarkset") %>%
  separate(
    col = benchmarkset,
    into = c("bench_version", "ref", "bench_type"),
    sep = "_"
  )

# Visualize difficult region inclusion
diff_inclusion_plt <- strat_metrics_df %>%
  filter(bench_version %in% c("v5q", "v421", "v06")) %>%
  ggplot(aes(x = strat_name, y = pct_of_dip, fill = bench_version)) +
  geom_col(position = "dodge") +
  facet_wrap(~ bench_type, scales = "free") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(
    title = "Difficult Region Coverage by Benchmark Version",
    x = "Stratification",
    y = "% of Benchmark in Difficult Region",
    fill = "Version"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

diff_inclusion_plt
```

## Troubleshooting

### Issue: Missing stratification metrics

**Error:** `FileNotFoundError: results/strat_metrics/{benchmark}/stratification_coverage_table.csv`

**Solution:**
```bash
# Ensure stratifications are downloaded
snakemake --cores 1 resources/stratifications/GRCh38_TR.bed.gz

# Run metrics generation
snakemake --cores 8 results/strat_metrics/v5q_grch38_smvar/stratification_coverage_table.csv
```

### Issue: Empty variant counts

**Symptom:** `total_variants = 0` for all stratifications

**Causes:**
1. Variant table missing STRAT_IDS annotations
2. Variants filtered out (not in benchmark regions)

**Solution:**
```bash
# Verify variant table has annotations
zcat results/variant_tables/v5q_grch38_smvar/variants.tsv.gz | head -1 | grep STRAT_IDS

# Regenerate variant table if needed
snakemake --cores 8 --forcerun generate_var_table results/variant_tables/v5q_grch38_smvar/variants.tsv
```

### Issue: Reference mismatch

**Error:** `Reference mismatch: benchmark uses GRCh38, stratification uses GRCh37`

**Solution:**
- Ensure config.yaml correctly specifies `ref:` for each benchmark
- Stratifications are automatically matched to benchmark reference

## Pipeline Architecture

```
Variant Table Generation (existing)
    ↓
    ├→ annotate_vcf_stratifications (adds STRAT_IDS)
    ├→ annotate_vcf_regions (adds REGION_IDS)
    └→ generate_var_table (TSV output)

New: Stratification Metrics
    ↓
    ├→ materialize_stratification (link to downloaded BED)
    ├→ compute_stratification_size (total stratification size)
    ├→ compute_stratification_metrics (overlap with benchmark)
    └→ aggregate_stratification_metrics (CSV table)

New: Variant Counts
    ↓
    ├→ count_variants_by_stratification (parse STRAT_IDS)
    ├→ summarize_variant_counts (aggregate by stratification)
    └→ combine_metrics_and_counts (unified table with density)
```

## Related Files

- **Snakemake rules:**
  - `workflow/rules/strat_metrics.smk` - Stratification coverage metrics
  - `workflow/rules/var_counts.smk` - Variant counting
  - `workflow/rules/var_tables.smk` - Variant table generation (existing)
  - `workflow/rules/exclusions.smk` - Exclusion metrics (template for strat_metrics)

- **Python scripts:**
  - `workflow/scripts/count_variants_by_strat.py` - Count variants by stratification
  - `workflow/scripts/summarize_var_counts.py` - Aggregate counts
  - `workflow/scripts/combine_metrics_counts.py` - Merge coverage and counts

- **R scripts:**
  - `scripts/plot_difficult_region_cumulative.R` - Generate cumulative plots

- **Analysis documents:**
  - `analysis/benchmarkset-characterization.qmd` - Main analysis (update with metrics)

## Further Reading

- GIAB Genome Stratifications: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/
- Benchmark methodology: See manuscript files in `manuscript/`
