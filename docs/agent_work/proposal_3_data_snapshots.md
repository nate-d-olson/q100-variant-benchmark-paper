# Proposal 3: Analysis-Ready Data Snapshots with Version Control

## Overview

Generate curated, analysis-ready data snapshots from Snakemake pipeline outputs and store them as versioned, compressed R objects. This proposal creates a middle tier between raw pipeline outputs and analysis code, providing stable data artifacts that are optimized for manuscript figure generation.

## Key Principles

- **Stable artifacts** - Data snapshots don't change unless explicitly regenerated
- **Version control friendly** - Compressed binary formats with version tags
- **Optimized for analysis** - Pre-joined, pre-factored, ready for plotting
- **Explicit refresh** - Analysts control when to update data

## Technical Approach

### 1. Snapshot Generation Framework

Create a script that generates analysis-ready snapshots from pipeline outputs:

```r
# R/generate_snapshots.R (NEW FILE)

#' Generate Analysis-Ready Data Snapshots
#' 
#' Creates compressed, versioned snapshots of pipeline outputs optimized
#' for manuscript analysis. Snapshots are stored in analysis/snapshots/
#' and loaded automatically by notebooks.
#' 
#' @param snapshot_version Version tag (default: auto-generated from date)
#' @param results_dir Pipeline results directory
#' @param output_dir Snapshot output directory
#' 
#' @export
generate_analysis_snapshots <- function(
  snapshot_version = NULL,
  results_dir = NULL,
  output_dir = NULL
) {
  # Setup paths
  if (is.null(results_dir)) {
    results_dir <- here::here("results")
  }
  
  if (is.null(output_dir)) {
    output_dir <- here::here("analysis/snapshots")
  }
  
  if (is.null(snapshot_version)) {
    snapshot_version <- format(Sys.Date(), "%Y%m%d")
  }
  
  fs::dir_create(output_dir)
  
  message("=" %R% 60)
  message("Generating Analysis Snapshots v", snapshot_version)
  message("=" %R% 60)
  
  # Load and prepare all data
  snapshots <- list()
  
  ## 1. Primary Metrics ----
  message("\n[1/6] Loading genomic context metrics...")
  snapshots$genomic_metrics <- load_genomic_context_metrics(
    results_dir = results_dir
  ) %>%
    # Add derived columns useful for analysis
    dplyr::mutate(
      benchmark_id = paste(bench_version, ref, var_type, sep = "_"),
      # Proportion instead of percentage for some calculations
      prop_of_context = pct_of_context / 100,
      prop_of_bench = pct_of_bench / 100,
      # Variant rate per kilobase for finer resolution
      variant_density_per_kb = variant_density_per_mb / 1000
    ) %>%
    # Ensure factors are set
    dplyr::mutate(
      bench_version = factor(bench_version, levels = c("v0.6", "v4.2.1", "v5.0q")),
      ref = factor(ref, levels = c("GRCh37", "GRCh38", "CHM13v2.0")),
      var_type = factor(var_type, levels = c("smvar", "stvar")),
      context_name = factor(context_name, 
                           levels = c("HP", "MAP", "SD", "SD10kb", "TR", "TR10kb"))
    )
  
  ## 2. Benchmark Regions ----
  message("[2/6] Loading benchmark regions...")
  snapshots$benchmark_regions <- load_benchmark_regions() %>%
    # Add chromosome ordering
    dplyr::mutate(
      benchmark_id = paste(bench_version, ref, bench_type, sep = "_"),
      chrom_num = as.integer(stringr::str_remove(chrom, "chr"))
    ) %>%
    dplyr::arrange(bench_version, ref, chrom_num, start)
  
  ## 3. Reference Sizes ----
  message("[3/6] Loading reference genome sizes...")
  snapshots$reference_sizes <- load_reference_sizes(
    results_dir = results_dir
  ) %>%
    dplyr::mutate(
      ref = factor(ref, levels = c("GRCh37", "GRCh38", "CHM13v2.0"))
    )
  
  ## 4. Reference Summaries ----
  message("[4/6] Computing reference summaries...")
  snapshots$reference_summary <- snapshots$reference_sizes %>%
    dplyr::group_by(ref) %>%
    dplyr::summarise(
      total_length_bp = sum(length),
      total_ns = sum(ns),
      total_asm_bp = sum(asm_bp),
      n_chromosomes = dplyr::n(),
      pct_ns = (total_ns / total_length_bp) * 100,
      .groups = "drop"
    )
  
  ## 5. Benchmark Summaries ----
  message("[5/6] Computing benchmark summaries...")
  snapshots$benchmark_summary <- snapshots$genomic_metrics %>%
    dplyr::group_by(benchmark_id, bench_version, ref, var_type) %>%
    dplyr::summarise(
      n_contexts = dplyr::n(),
      total_variants = sum(total_variants),
      total_snv = sum(snv_count),
      total_indel = sum(indel_count),
      total_del = sum(del_count),
      total_ins = sum(ins_count),
      total_complex = sum(complex_count),
      total_other = sum(other_count),
      mean_variant_density = mean(variant_density_per_mb),
      .groups = "drop"
    )
  
  ## 6. Exclusions (v5.0q only) ----
  message("[6/6] Loading exclusion metrics...")
  snapshots$exclusions <- tryCatch(
    {
      load_exclusion_metrics(results_dir = results_dir) %>%
        dplyr::mutate(
          benchmark_id = paste(bench_version, ref, var_type, sep = "_")
        )
    },
    error = function(e) {
      message("Note: No exclusion data (only available for v5.0q)")
      tibble::tibble()
    }
  )
  
  ## 7. Pipeline Metadata ----
  message("\nAdding pipeline metadata...")
  snapshots$metadata <- list(
    snapshot_version = snapshot_version,
    snapshot_date = Sys.time(),
    pipeline_config = parse_pipeline_config(),
    r_version = R.version.string,
    package_versions = sessioninfo::package_info() %>%
      dplyr::filter(package %in% c(
        "tidyverse", "here", "vroom", "glue", "assertthat"
      )) %>%
      dplyr::select(package, ondiskversion)
  )
  
  # Save snapshot
  snapshot_file <- fs::path(
    output_dir,
    glue::glue("analysis_snapshot_v{snapshot_version}.rds")
  )
  
  message("\nSaving snapshot to:", snapshot_file)
  saveRDS(snapshots, snapshot_file, compress = "xz")
  
  # Create symlink to latest
  latest_link <- fs::path(output_dir, "latest.rds")
  if (fs::file_exists(latest_link)) {
    fs::file_delete(latest_link)
  }
  fs::link_create(snapshot_file, latest_link)
  
  # Generate snapshot manifest
  manifest <- tibble::tibble(
    dataset = names(snapshots),
    n_rows = purrr::map_int(snapshots, ~ {
      if (is.data.frame(.x)) nrow(.x) else NA_integer_
    }),
    n_cols = purrr::map_int(snapshots, ~ {
      if (is.data.frame(.x)) ncol(.x) else NA_integer_
    }),
    size_bytes = purrr::map_dbl(snapshots, ~ {
      object.size(.x)
    })
  )
  
  manifest_file <- fs::path(
    output_dir,
    glue::glue("manifest_v{snapshot_version}.csv")
  )
  vroom::vroom_write(manifest, manifest_file)
  
  message("\n" %R% 60)
  message("✓ Snapshot generated successfully!")
  message("  Version:", snapshot_version)
  message("  Files:", length(snapshots), "datasets")
  message("  Size:", format(fs::file_size(snapshot_file), units = "auto"))
  message("=" %R% 60)
  
  return(invisible(snapshot_file))
}
```

### 2. Snapshot Loading Interface

Create simple functions to load snapshots in notebooks:

```r
# R/load_snapshots.R (NEW FILE)

#' Load Analysis Snapshot
#' 
#' Loads a versioned data snapshot. By default, loads the latest snapshot.
#' 
#' @param version Snapshot version (default: "latest")
#' @param snapshot_dir Snapshot directory
#' 
#' @return Named list with all snapshot datasets
#' 
#' @examples
#' \dontrun{
#' # Load latest snapshot
#' data <- load_snapshot()
#' 
#' # Access datasets
#' data$genomic_metrics
#' data$benchmark_summary
#' 
#' # Load specific version
#' data <- load_snapshot("20250115")
#' }
#' 
#' @export
load_snapshot <- function(version = "latest", snapshot_dir = NULL) {
  if (is.null(snapshot_dir)) {
    snapshot_dir <- here::here("analysis/snapshots")
  }
  
  # Resolve version
  if (version == "latest") {
    snapshot_file <- fs::path(snapshot_dir, "latest.rds")
  } else {
    snapshot_file <- fs::path(
      snapshot_dir,
      glue::glue("analysis_snapshot_v{version}.rds")
    )
  }
  
  # Check file exists
  if (!fs::file_exists(snapshot_file)) {
    available <- fs::dir_ls(snapshot_dir, glob = "*.rds") %>%
      fs::path_file() %>%
      stringr::str_remove("analysis_snapshot_v") %>%
      stringr::str_remove(".rds")
    
    stop(
      glue::glue(
        "Snapshot version '{version}' not found.\n",
        "Available versions: {paste(available, collapse = ', ')}\n",
        "Run generate_analysis_snapshots() to create a new snapshot."
      ),
      call. = FALSE
    )
  }
  
  # Load snapshot
  message(glue::glue("Loading snapshot from: {snapshot_file}"))
  snapshot <- readRDS(snapshot_file)
  
  # Display metadata
  metadata <- snapshot$metadata
  message(glue::glue(
    "Snapshot version: {metadata$snapshot_version}\n",
    "Generated: {metadata$snapshot_date}\n",
    "Datasets: {length(snapshot) - 1}"  # Exclude metadata
  ))
  
  return(snapshot)
}

#' List Available Snapshots
#' 
#' @param snapshot_dir Snapshot directory
#' 
#' @return Tibble with snapshot information
#' 
#' @export
list_snapshots <- function(snapshot_dir = NULL) {
  if (is.null(snapshot_dir)) {
    snapshot_dir <- here::here("analysis/snapshots")
  }
  
  if (!fs::dir_exists(snapshot_dir)) {
    message("No snapshots directory found. Run generate_analysis_snapshots() first.")
    return(tibble::tibble())
  }
  
  # Find all snapshot files
  snapshot_files <- fs::dir_ls(snapshot_dir, glob = "*analysis_snapshot*.rds")
  
  if (length(snapshot_files) == 0) {
    message("No snapshots found. Run generate_analysis_snapshots() to create one.")
    return(tibble::tibble())
  }
  
  # Extract metadata from each
  snapshot_info <- snapshot_files %>%
    purrr::map_dfr(function(file) {
      meta <- readRDS(file)$metadata
      tibble::tibble(
        version = meta$snapshot_version,
        created = meta$snapshot_date,
        file = fs::path_file(file),
        size_mb = fs::file_size(file) / 1024^2,
        n_benchmarks = meta$pipeline_config$num_benchmarks
      )
    }) %>%
    dplyr::arrange(dplyr::desc(created))
  
  return(snapshot_info)
}
```

### 3. Makefile Integration

Add snapshot generation to the Makefile:

```makefile
# Add to existing Makefile

.PHONY: snapshot
snapshot: ## Generate analysis-ready data snapshot
	@echo "Generating analysis snapshots from pipeline outputs..."
	Rscript -e "source('R/generate_snapshots.R'); generate_analysis_snapshots()"

.PHONY: list-snapshots
list-snapshots: ## List available data snapshots
	@echo "Available data snapshots:"
	Rscript -e "source('R/load_snapshots.R'); print(list_snapshots())"
```

### 4. Notebook Template

Standardized setup for all Quarto notebooks:

```r
# Standard notebook setup (add to all *.qmd files)

# Load snapshot (not raw pipeline outputs)
source(here::here("R/load_snapshots.R"))

# Load latest snapshot
snapshot <- load_snapshot()
# Loading snapshot from: analysis/snapshots/latest.rds
# Snapshot version: 20250207
# Generated: 2025-02-07 10:30:45
# Datasets: 6

# Extract datasets for analysis
genomic_metrics <- snapshot$genomic_metrics
benchmark_summary <- snapshot$benchmark_summary
benchmark_regions <- snapshot$benchmark_regions
reference_sizes <- snapshot$reference_sizes
exclusions <- snapshot$exclusions

# All data is pre-processed and ready for analysis!
```

## Implementation Steps

### Phase 1: Snapshot Generation (Week 1)
1. Create `R/generate_snapshots.R` with generation logic
2. Create `R/load_snapshots.R` with loading functions
3. Add snapshot generation tests
4. Document snapshot structure

### Phase 2: Initial Snapshot (Week 2)
1. Run pipeline to generate fresh outputs
2. Generate first official snapshot
3. Validate snapshot completeness
4. Archive snapshot with version tag

### Phase 3: Notebook Migration (Week 3)
1. Update `benchmarkset_characterization.qmd` to use snapshots
2. Update `benchmark_difficult.qmd` to use snapshots
3. Update `external_evaluation.qmd` (external data separate)
4. Verify all figures reproduce correctly

### Phase 4: Workflow Integration (Week 4)
1. Add `make snapshot` target
2. Create snapshot refresh workflow
3. Document snapshot lifecycle
4. Add snapshot versioning guide

## Benefits

### Stability and Reproducibility
- **Frozen artifacts** - Data doesn't change between analysis runs
- **Version control** - Track which snapshot produced each manuscript version
- **Portability** - Single file contains all analysis data
- **Rollback capability** - Can always return to previous snapshot

### Performance
- **Instant loading** - Pre-processed data loads in <1 second
- **Optimized format** - Compressed RDS with factors pre-set
- **No parsing** - Skip CSV reading and type conversion
- **Pre-joined data** - Common joins already performed

### Collaboration
- **Easy sharing** - Send single snapshot file to collaborators
- **Consistent results** - Everyone uses same data version
- **Reduced pipeline dependencies** - Analysts don't need to run full pipeline
- **Clear versioning** - Explicit snapshot versions in filenames

### Analysis-Ready
- **Pre-factored** - All categorical variables with correct levels
- **Pre-joined** - Common joins already performed
- **Derived columns** - Useful calculations pre-computed
- **Metadata included** - Pipeline config and package versions

## Drawbacks and Mitigations

### Stale Data Risk
- **Risk**: Analysts may forget to regenerate snapshots
- **Mitigation**:
  - Display snapshot age when loading
  - Add `make snapshot` to workflow documentation
  - Include snapshot generation in CI/CD
  - Warn if snapshot is >7 days old

### Storage Requirements
- **Risk**: Multiple snapshot versions consume disk space
- **Mitigation**:
  - Compress with `xz` compression (typically 80-90% reduction)
  - Keep only last 3 versions by default
  - Document cleanup procedure
  - Snapshots are much smaller than raw pipeline outputs

### Two-Step Workflow
- **Risk**: Requires running pipeline then generating snapshot
- **Mitigation**:
  - Integrate snapshot generation into pipeline (`onsuccess` hook)
  - Provide `make pipeline-and-snapshot` combined target
  - Document workflow clearly

### Binary Format
- **Risk**: RDS files are not human-readable
- **Mitigation**:
  - Include manifest CSV showing snapshot contents
  - Provide `snapshot_info()` function to inspect
  - Keep raw CSVs in results/ for debugging

## File Structure Changes

```
q100-variant-benchmark-paper/
├── R/
│   ├── data_loading.R          # Existing (still used by snapshot generation)
│   ├── generate_snapshots.R    # NEW: Snapshot generation
│   └── load_snapshots.R        # NEW: Snapshot loading
├── analysis/
│   ├── snapshots/              # NEW: Versioned snapshots
│   │   ├── analysis_snapshot_v20250207.rds
│   │   ├── manifest_v20250207.csv
│   │   ├── latest.rds          # Symlink to latest
│   │   └── README.md           # Snapshot documentation
│   └── *.qmd                   # Updated to load snapshots
├── Makefile                    # Add snapshot targets
└── .gitignore                  # Add analysis/snapshots/*.rds
```

## Example Workflow

### Analyst Workflow (Typical)

```bash
# 1. Load latest snapshot in R notebook
# (no pipeline run needed)
data <- load_snapshot()

# 2. Create figures and analysis
# (all data pre-loaded and ready)

# 3. When pipeline updates:
make snapshot  # Generate new snapshot
# Then re-run notebooks to get latest data
```

### Full Workflow (When Pipeline Changes)

```bash
# 1. Run Snakemake pipeline
snakemake --cores 4 --sdm conda

# 2. Generate fresh snapshot
make snapshot
# Generating Analysis Snapshots v20250207
# [1/6] Loading genomic context metrics...
# [2/6] Loading benchmark regions...
# ✓ Snapshot generated successfully!

# 3. Verify snapshot
make list-snapshots
# version     created              size_mb  n_benchmarks
# 20250207    2025-02-07 10:30:45  12.4     15

# 4. Run analysis notebooks
quarto render analysis/*.qmd
```

### Snapshot Management

```r
# List available snapshots
list_snapshots()
# version     created              size_mb  n_benchmarks
# 20250207    2025-02-07 10:30:45  12.4     15
# 20250115    2025-01-15 14:22:10  11.8     12
# 20241220    2024-12-20 09:15:33  10.2     10

# Load specific version
data_jan <- load_snapshot("20250115")

# Compare snapshots
current <- load_snapshot("20250207")
previous <- load_snapshot("20250115")

# How many more variants in new benchmarks?
sum(current$benchmark_summary$total_variants) -
  sum(previous$benchmark_summary$total_variants)
```

## Comparison to Other Proposals

**Advantages over Proposal 1 (Caching):**
- Version control and explicit versioning
- Pre-processed and optimized for analysis
- Portable and shareable
- Stable across sessions

**Advantages over Proposal 2 (Data Access Layer):**
- Zero runtime overhead (data pre-loaded)
- Simple interface (single load function)
- Faster implementation
- No validation runtime cost

**Disadvantages:**
- Requires explicit refresh step
- Additional storage for snapshots
- Binary format less transparent

## Recommendation

**Best suited for:** Teams in manuscript preparation phase who need stable, reproducible data artifacts. Ideal when notebooks are being run repeatedly to finalize figures and collaborators need consistent data versions.

**Priority Level:** HIGH - Provides excellent balance of simplicity, performance, and reproducibility for manuscript finalization.

## Additional Considerations

### Integration with Git LFS

For large snapshots (>50 MB), consider using Git LFS:

```bash
# Track snapshot files with Git LFS
git lfs track "analysis/snapshots/*.rds"
git add .gitattributes
```

### Snapshot Naming Convention

Use semantic versioning for major snapshot milestones:

- `v1.0.0` - Initial submission snapshot
- `v1.1.0` - Revision with new benchmarks
- `v2.0.0` - Major pipeline update

For development snapshots, use dates: `20250207`, `20250208`, etc.
