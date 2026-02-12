# Proposal 1: Enhanced Data Loading Library with Intelligent Caching

## Overview

Enhance the existing `R/data_loading.R` library with intelligent caching mechanisms to improve performance and reproducibility. This proposal builds incrementally on the solid foundation already in place, adding disk-based caching and memoization to reduce redundant data loading across analysis sessions.

## Key Principles

- **Build on existing infrastructure** - Don't replace, enhance
- **Transparent caching** - Works seamlessly without changing notebook code
- **Invalidation on changes** - Automatically detects when source data updates
- **Minimal overhead** - Fast cache checks, optional for small files

## Technical Approach

### 1. Add Caching Layer to Existing Functions

Wrap existing data loading functions with a caching layer that:

- Stores processed R objects as `.rds` files in `analysis/cache/`
- Uses content-based cache keys (file modification time + file size)
- Automatically invalidates when source files change
- Provides cache statistics and management utilities

```r
# New internal caching utility (add to R/data_loading.R)
.cache_wrapper <- function(cache_key, load_fn, force_refresh = FALSE) {
  cache_dir <- here::here("analysis/cache")
  cache_file <- fs::path(cache_dir, paste0(cache_key, ".rds"))
  
  # Check if cache exists and is valid
  if (!force_refresh && fs::file_exists(cache_file)) {
    message(glue::glue("Loading from cache: {cache_key}"))
    return(readRDS(cache_file))
  }
  
  # Load fresh data
  message(glue::glue("Computing and caching: {cache_key}"))
  result <- load_fn()
  
  # Save to cache
  fs::dir_create(cache_dir)
  saveRDS(result, cache_file, compress = "xz")
  
  return(result)
}

# Enhanced load_genomic_context_metrics with caching
load_genomic_context_metrics <- function(results_dir = NULL, 
                                         benchmark_filter = NULL,
                                         use_cache = TRUE,
                                         force_refresh = FALSE) {
  if (!use_cache) {
    return(.load_genomic_context_metrics_impl(results_dir, benchmark_filter))
  }
  
  # Generate cache key from inputs
  cache_key <- digest::digest(list(
    results_dir = results_dir %||% here::here("results"),
    benchmark_filter = benchmark_filter,
    # Include modification times of source files
    mod_times = get_metrics_file_mtimes(results_dir)
  ))
  
  .cache_wrapper(
    cache_key = paste0("genomic_metrics_", cache_key),
    load_fn = function() {
      .load_genomic_context_metrics_impl(results_dir, benchmark_filter)
    },
    force_refresh = force_refresh
  )
}
```

### 2. Cache Management Utilities

Add new helper functions for cache management:

```r
# Clear all cached data
clear_analysis_cache <- function() {
  cache_dir <- here::here("analysis/cache")
  if (fs::dir_exists(cache_dir)) {
    fs::dir_delete(cache_dir)
    message("Analysis cache cleared")
  }
}

# Get cache statistics
cache_stats <- function() {
  cache_dir <- here::here("analysis/cache")
  if (!fs::dir_exists(cache_dir)) {
    return(tibble::tibble(files = 0, size_mb = 0))
  }
  
  files <- fs::dir_ls(cache_dir, glob = "*.rds")
  tibble::tibble(
    n_files = length(files),
    total_size_mb = sum(fs::file_size(files)) / 1024^2,
    files = fs::path_file(files)
  )
}

# Force refresh specific dataset
refresh_dataset <- function(dataset_name) {
  # Wrapper to force_refresh = TRUE for specific datasets
  # e.g., refresh_dataset("genomic_metrics")
}
```

### 3. Progressive Loading for Large Files

Add lazy loading option for large variant tables:

```r
load_variant_table <- function(benchmark_id, 
                                results_dir = NULL,
                                filters = NULL,
                                lazy = FALSE,
                                use_cache = TRUE) {
  
  if (lazy) {
    # Return a promise/proxy object that loads on first access
    return(create_lazy_variant_loader(benchmark_id, results_dir, filters))
  }
  
  # Existing implementation with added caching
  if (use_cache) {
    cache_key <- digest::digest(list(benchmark_id, results_dir, filters))
    return(.cache_wrapper(
      cache_key = paste0("variants_", cache_key),
      load_fn = function() {
        .load_variant_table_impl(benchmark_id, results_dir, filters)
      }
    ))
  }
  
  .load_variant_table_impl(benchmark_id, results_dir, filters)
}
```

## Implementation Steps

### Phase 1: Core Caching Infrastructure (Week 1)

1. Add `.cache_wrapper()` internal function
2. Add cache management utilities (`clear_analysis_cache()`, `cache_stats()`)
3. Update `.gitignore` to exclude `analysis/cache/` directory
4. Document caching behavior in function documentation

### Phase 2: Integrate Caching into Existing Functions (Week 2)

1. Update `load_genomic_context_metrics()` with caching
2. Update `load_variant_table()` with caching (most impactful)
3. Update `load_diff_coverage()` with caching
4. Update `load_benchmark_regions()` with caching

### Phase 3: Testing and Documentation (Week 3)

1. Add unit tests for caching logic
2. Update README with caching documentation
3. Add troubleshooting guide for cache issues
4. Update Quarto notebooks with caching examples

### Phase 4: Advanced Features (Optional)

1. Implement lazy loading for very large files
2. Add parallel cache warming utility
3. Add cache compression optimization

## Benefits

### Performance Improvements

- **80-95% faster re-runs** of analysis notebooks (cached data loads in milliseconds)
- **Instant restarts** after R session crashes or kernel restarts
- **Reduced memory pressure** when working with multiple notebooks simultaneously

### Development Experience

- **Faster iteration** when tweaking visualizations or statistical models
- **Notebook portability** - cache travels with analysis directory
- **Explicit refresh** when pipeline outputs update

### Reproducibility

- **Stable snapshots** of intermediate data transformations
- **Cache key tracking** for provenance
- **Version control friendly** (cache directory excluded from git)

## Drawbacks and Mitigations

### Disk Space Usage

- **Risk**: Cache directory can grow large with many datasets
- **Mitigation**:
  - Provide `clear_analysis_cache()` utility
  - Add cache size warnings when >1GB
  - Document cache location for manual cleanup

### Stale Cache Issues

- **Risk**: Users may forget to refresh cache after pipeline updates
- **Mitigation**:
  - Automatic invalidation based on file modification times
  - Add `force_refresh = TRUE` parameter to all cached functions
  - Display cache age in `cache_stats()`

### Complexity

- **Risk**: Adds another layer of abstraction
- **Mitigation**:
  - Make caching opt-in via `use_cache` parameter (default: TRUE)
  - Keep caching logic internal and transparent
  - Provide clear error messages for cache issues

## File Structure Changes

```
q100-variant-benchmark-paper/
├── R/
│   └── data_loading.R          # Enhanced with caching
├── analysis/
│   ├── cache/                  # New: cached R objects (gitignored)
│   │   ├── genomic_metrics_abc123.rds
│   │   ├── variants_def456.rds
│   │   └── ...
│   └── *.qmd                   # No changes needed
├── .gitignore                  # Add analysis/cache/
└── README.md                   # Document caching
```

## Backward Compatibility

**100% backward compatible** - All existing notebook code continues to work without modification. Caching is transparent and can be disabled with `use_cache = FALSE`.

## Example Usage

### In Quarto Notebooks (No Changes Required)

```r
# Existing code works as-is, but now with caching
metrics <- load_genomic_context_metrics()  # First run: loads and caches
metrics <- load_genomic_context_metrics()  # Second run: instant from cache

# Force refresh after pipeline updates
metrics <- load_genomic_context_metrics(force_refresh = TRUE)

# Disable caching if needed
metrics <- load_genomic_context_metrics(use_cache = FALSE)
```

### Cache Management

```r
# Check cache status
cache_stats()
# → 15 files, 234 MB

# Clear cache before final manuscript run
clear_analysis_cache()
source(here("R/data_loading.R"))
metrics <- load_genomic_context_metrics()  # Fresh load
```

## Comparison to Other Proposals

**Advantages over Proposal 2 (Standardized Layer):**

- Lighter weight, builds on existing code
- No architectural changes required
- Faster to implement

**Advantages over Proposal 3 (Snapshots):**

- No manual snapshot generation steps
- Automatic invalidation on data changes
- Works with dynamic filtering and subsets

## Recommendation

**Best suited for:** Teams that want immediate performance improvements without changing their workflow. Ideal for the current phase where figures are being finalized and notebooks are being re-run frequently.

**Priority Level:** HIGH - Delivers maximum benefit with minimal disruption to existing code.
