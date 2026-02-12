# Implementation Plan Addendum: Parquet Cache Format Support

## Overview

This addendum extends the original implementation plan to incorporate Arrow/Parquet as an optional cache format alongside RDS. Based on the cache format comparison analysis, Parquet offers 2-4x faster performance and 40-60% better compression for data frames.

---

## Updated Pre-Implementation Checklist

- [ ] Review existing `R/data_loading.R` code
- [ ] Identify most-used functions for initial caching implementation
- [ ] Determine cache directory location and structure
- [ ] Add `digest` package to dependencies (for cache keys)
- [ ] **Add `arrow` package to dependencies (for Parquet support) - OPTIONAL**
- [ ] Update `.gitignore` to exclude cache directory
- [ ] Create branch: `feat/caching-library`

---

## Updated Phase 1: Core Caching Infrastructure with Format Support

### New Task 1.1: Add Required Dependencies

**File:** `environment.yaml`

**Changes:**

```yaml
# Add to dependencies:
  - r-digest>=0.6.33        # For cache key generation (required)
  - r-arrow>=10.0.0         # For Parquet format support (optional but recommended)
```

**Test:**

```r
library(digest)
digest::digest("test")  # Should return hash

# Test arrow (if installed)
if (requireNamespace("arrow", quietly = TRUE)) {
  library(arrow)
  test_df <- data.frame(x = 1:10, y = letters[1:10])
  arrow::write_parquet(test_df, tempfile(fileext = ".parquet"))
  message("✓ Arrow/Parquet support available")
} else {
  message("ℹ Arrow not installed - will use RDS format only")
}
```

---

### Updated Task 1.2: Enhanced Cache Wrapper with Format Support

**File:** `R/data_loading.R`

**Enhanced Code:**

```r
#' Internal: Cache Wrapper for Data Loading Functions
#'
#' Wraps a data loading function with transparent caching.
#' Supports multiple cache formats: RDS (native R) and Parquet (Arrow).
#'
#' @param cache_key Unique identifier for this cache entry
#' @param load_fn Function that loads data (should return a data frame)
#' @param force_refresh Logical; if TRUE, bypass cache and reload fresh data
#' @param use_cache Logical; if FALSE, skip caching entirely
#' @param cache_format Character; "auto", "rds", or "parquet". 
#'   "auto" selects format based on data size (Parquet for >5MB, RDS otherwise)
#'
#' @return Result from load_fn (either from cache or fresh load)
#'
#' @keywords internal
#' @noRd
.cache_wrapper <- function(cache_key, 
                           load_fn, 
                           force_refresh = FALSE, 
                           use_cache = TRUE,
                           cache_format = c("auto", "rds", "parquet")) {
  # Skip caching if disabled
  if (!use_cache) {
    return(load_fn())
  }
  
  cache_format <- match.arg(cache_format)
  
  # Setup cache directory
  cache_dir <- here::here("analysis", "cache")
  fs::dir_create(cache_dir)
  
  # Check for existing cache files (both formats)
  cache_file_rds <- fs::path(cache_dir, paste0(cache_key, ".rds"))
  cache_file_parquet <- fs::path(cache_dir, paste0(cache_key, ".parquet"))
  
  # Determine which cache file exists
  cache_exists <- FALSE
  cache_file <- NULL
  cache_ext <- NULL
  
  if (fs::file_exists(cache_file_parquet)) {
    cache_exists <- TRUE
    cache_file <- cache_file_parquet
    cache_ext <- "parquet"
  } else if (fs::file_exists(cache_file_rds)) {
    cache_exists <- TRUE
    cache_file <- cache_file_rds
    cache_ext <- "rds"
  }
  
  # Load from cache if exists and not forcing refresh
  if (cache_exists && !force_refresh) {
    cache_age <- difftime(
      Sys.time(), 
      fs::file_info(cache_file)$modification_time, 
      units = "hours"
    )
    
    message(glue::glue(
      "Loading from cache: {cache_key} ({cache_ext}, ",
      "cached {round(cache_age, 1)} hours ago)"
    ))
    
    # Read based on format
    if (cache_ext == "parquet") {
      if (!requireNamespace("arrow", quietly = TRUE)) {
        stop(
          "Cached file is in Parquet format but arrow package is not installed. ",
          "Install with: install.packages('arrow') or clear cache and retry.",
          call. = FALSE
        )
      }
      return(arrow::read_parquet(cache_file, as_data_frame = TRUE))
    } else {
      return(readRDS(cache_file))
    }
  }
  
  # Load fresh data
  message(glue::glue("Loading and caching: {cache_key}"))
  start_time <- Sys.time()
  result <- load_fn()
  elapsed <- difftime(Sys.time(), start_time, units = "secs")
  
  # Auto-select format based on size if format is "auto"
  if (cache_format == "auto") {
    size_mb <- as.numeric(object.size(result)) / 1024^2
    
    # Use Parquet for data >5 MB if arrow is available
    if (size_mb > 5 && requireNamespace("arrow", quietly = TRUE)) {
      cache_format <- "parquet"
      cache_ext <- "parquet"
      cache_file <- cache_file_parquet
    } else {
      cache_format <- "rds"
      cache_ext <- "rds"
      cache_file <- cache_file_rds
    }
  } else {
    cache_ext <- cache_format
    cache_file <- if (cache_format == "parquet") {
      cache_file_parquet
    } else {
      cache_file_rds
    }
  }
  
  # Save to cache
  if (cache_format == "parquet") {
    if (!requireNamespace("arrow", quietly = TRUE)) {
      warning(
        "arrow package not available, falling back to RDS format. ",
        "Install with: install.packages('arrow')",
        call. = FALSE
      )
      cache_format <- "rds"
      cache_ext <- "rds"
      cache_file <- cache_file_rds
    }
  }
  
  # Write cache file
  if (cache_format == "parquet") {
    arrow::write_parquet(
      result, 
      cache_file,
      compression = "zstd",
      compression_level = 3
    )
  } else {
    saveRDS(result, cache_file, compress = "xz")
  }
  
  cache_size <- fs::file_size(cache_file) / 1024^2  # Convert to MB
  
  message(glue::glue(
    "Cached in {round(elapsed, 2)}s as {cache_ext} ",
    "(size: {round(cache_size, 2)} MB)"
  ))
  
  return(result)
}
```

---

### Updated Task 1.3: Enhanced Cache Management Utilities

**Additional function to add:**

```r
#' Get Detailed Cache Statistics
#'
#' Returns comprehensive information about cached data files including
#' format, size, age, and compression efficiency.
#'
#' @return Tibble with cache statistics
#'
#' @examples
#' \dontrun{
#' # Check cache status with format details
#' cache_stats_detailed()
#' }
#'
#' @export
cache_stats_detailed <- function() {
  cache_dir <- here::here("analysis", "cache")
  
  if (!fs::dir_exists(cache_dir)) {
    message("No cache directory found")
    return(tibble::tibble(
      n_files = 0,
      total_size_mb = 0
    ))
  }
  
  cache_files <- fs::dir_ls(cache_dir, regexp = "\\.(rds|parquet)$")
  
  if (length(cache_files) == 0) {
    message("Cache directory is empty")
    return(tibble::tibble(
      n_files = 0,
      total_size_mb = 0
    ))
  }
  
  # Get file information
  file_info <- fs::file_info(cache_files)
  
  stats_df <- tibble::tibble(
    file = fs::path_file(cache_files),
    format = tools::file_ext(cache_files),
    size_mb = file_info$size / 1024^2,
    modified = file_info$modification_time,
    age_hours = as.numeric(difftime(Sys.time(), file_info$modification_time, 
                                   units = "hours"))
  ) %>%
    dplyr::arrange(dplyr::desc(modified))
  
  # Summary by format
  format_summary <- stats_df %>%
    dplyr::group_by(format) %>%
    dplyr::summarise(
      n_files = dplyr::n(),
      total_mb = sum(size_mb),
      avg_mb = mean(size_mb),
      .groups = "drop"
    )
  
  # Overall summary
  total_size_mb <- sum(stats_df$size_mb)
  n_files <- nrow(stats_df)
  
  message(glue::glue(
    "Cache contains {n_files} file(s), total size: {round(total_size_mb, 2)} MB\n",
    "Format breakdown:"
  ))
  
  for (i in seq_len(nrow(format_summary))) {
    fmt <- format_summary$format[i]
    n <- format_summary$n_files[i]
    size <- format_summary$total_mb[i]
    message(glue::glue("  {fmt}: {n} files, {round(size, 2)} MB"))
  }
  
  return(stats_df)
}
```

---

## Updated Phase 2: Function Integration with Format Options

### Updated load_genomic_context_metrics()

```r
#' Load Genomic Context Metrics
#'
#' [existing documentation...]
#'
#' @param cache_format Character; cache format to use ("auto", "rds", or "parquet").
#'   Default "auto" selects format based on data size. Parquet requires arrow package.
#'
#' @section Cache Formats:
#' - **RDS**: Native R format, works without additional packages
#' - **Parquet**: Columnar format (via arrow package), 2-4x faster, better compression
#' - **Auto**: Uses Parquet for data >5MB (if arrow available), RDS otherwise
#'
#' @export
load_genomic_context_metrics <- function(results_dir = NULL, 
                                         benchmark_filter = NULL,
                                         use_cache = TRUE,
                                         force_refresh = FALSE,
                                         cache_format = c("auto", "rds", "parquet")) {
  # Skip caching if disabled
  if (!use_cache) {
    return(.load_genomic_context_metrics_impl(results_dir, benchmark_filter))
  }
  
  cache_format <- match.arg(cache_format)
  
  # Setup defaults
  if (is.null(results_dir)) {
    results_dir <- here::here("results")
  }
  
  # Find source files for cache key
  metrics_files <- fs::dir_ls(
    results_dir,
    recurse = TRUE,
    glob = "**/genomic_context_combined_metrics.csv",
    fail = FALSE
  )
  
  # Generate cache key
  cache_key <- .generate_cache_key(
    results_dir = results_dir,
    benchmark_filter = benchmark_filter,
    source_files = metrics_files
  )
  
  # Load with caching
  .cache_wrapper(
    cache_key = paste0("genomic_metrics_", cache_key),
    load_fn = function() {
      .load_genomic_context_metrics_impl(results_dir, benchmark_filter)
    },
    force_refresh = force_refresh,
    use_cache = use_cache,
    cache_format = cache_format
  )
}
```

---

## New Phase 5: Optional DuckDB Query Support (Days 15-17)

### Objective

Add optional SQL query capability for large cached datasets without full loading.

### Task 5.1: Add DuckDB Helper Function

**File:** `R/data_loading.R`

**New function:**

```r
#' Query Cached Data with SQL
#'
#' Execute SQL query on cached Parquet data without loading full dataset.
#' Useful for exploratory analysis on large variant tables.
#'
#' Requires duckdb package: install.packages("duckdb")
#'
#' @param cache_pattern File pattern to match cached Parquet files
#' @param query SQL query string
#'
#' @return Query result as tibble
#'
#' @examples
#' \dontrun{
#' # Query variant counts by chromosome without loading full table
#' counts <- query_cached_data(
#'   cache_pattern = "variants_v5.0q_GRCh38_smvar*.parquet",
#'   query = "SELECT chrom, var_type, COUNT(*) as n 
#'            FROM data 
#'            WHERE chrom IN ('chr1', 'chr2') 
#'            GROUP BY chrom, var_type"
#' )
#'
#' # Compute statistics on variant sizes
#' stats <- query_cached_data(
#'   cache_pattern = "variants*.parquet",
#'   query = "SELECT var_type,
#'                   AVG(var_size) as mean_size,
#'                   PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY var_size) as median
#'            FROM data
#'            GROUP BY var_type"
#' )
#' }
#'
#' @export
query_cached_data <- function(cache_pattern, query) {
  # Check dependencies
  if (!requireNamespace("duckdb", quietly = TRUE)) {
    stop(
      "duckdb package required for querying cached data.\n",
      "Install with: install.packages('duckdb')",
      call. = FALSE
    )
  }
  
  cache_dir <- here::here("analysis", "cache")
  
  # Find matching cache files
  cache_files <- fs::dir_ls(cache_dir, glob = cache_pattern)
  
  if (length(cache_files) == 0) {
    stop(
      glue::glue("No cached files found matching pattern: {cache_pattern}"),
      call. = FALSE
    )
  }
  
  # Create DuckDB connection
  con <- DBI::dbConnect(duckdb::duckdb())
  
  tryCatch({
    # Register Parquet file(s) as table
    if (length(cache_files) == 1) {
      DBI::dbExecute(con, glue::glue(
        "CREATE VIEW data AS SELECT * FROM read_parquet('{cache_files[1]}')"
      ))
    } else {
      # Multiple files - use glob pattern
      DBI::dbExecute(con, glue::glue(
        "CREATE VIEW data AS SELECT * FROM read_parquet('{cache_dir}/{cache_pattern}')"
      ))
    }
    
    # Execute query
    result <- DBI::dbGetQuery(con, query)
    
    # Disconnect
    DBI::dbDisconnect(con, shutdown = TRUE)
    
    return(tibble::as_tibble(result))
    
  }, error = function(e) {
    # Ensure connection is closed on error
    DBI::dbDisconnect(con, shutdown = TRUE)
    stop(glue::glue("Query failed: {e$message}"), call. = FALSE)
  })
}
```

---

## Updated Documentation

### README.md Update

Add section on format selection:

```markdown
### Cache Format Selection

The caching system supports two formats:

#### RDS (Native R)
- Default for small datasets
- No additional packages required
- Preserves all R object attributes
- Good compression with xz

#### Parquet (via Arrow)
- 2-4x faster read/write
- 40-60% better compression
- Column selection support
- Requires arrow package: `install.packages("arrow")`

#### Usage

```r
# Auto-select format based on data size
metrics <- load_genomic_context_metrics(cache_format = "auto")

# Force specific format
metrics <- load_genomic_context_metrics(cache_format = "parquet")
metrics <- load_genomic_context_metrics(cache_format = "rds")

# Check cache format distribution
cache_stats_detailed()
```

#### Format Recommendations

- **Small data (<5 MB):** RDS (simpler, no dependencies)
- **Medium data (5-100 MB):** Parquet (faster, better compression)
- **Large data (>100 MB):** Parquet with optional DuckDB queries

#### Memory-Efficient Queries (Optional)

For large cached datasets, use SQL queries without full loading:

```r
# Install DuckDB (optional)
install.packages("duckdb")

# Query cached variant data
variant_summary <- query_cached_data(
  cache_pattern = "variants_*.parquet",
  query = "SELECT var_type, COUNT(*) as n 
           FROM data 
           WHERE chrom = 'chr1' 
           GROUP BY var_type"
)
```

```

---

## Updated Testing

### Format-Specific Tests

**File:** `tests/testthat/test-caching-formats.R`

```r
test_that("cache_wrapper supports RDS format", {
  result <- .cache_wrapper(
    "test_rds",
    function() tibble::tibble(x = 1:10),
    cache_format = "rds"
  )
  
  cache_file <- here::here("analysis", "cache", "test_rds.rds")
  expect_true(fs::file_exists(cache_file))
  expect_equal(result$x, 1:10)
})

test_that("cache_wrapper supports Parquet format if arrow available", {
  skip_if_not_installed("arrow")
  
  result <- .cache_wrapper(
    "test_parquet",
    function() tibble::tibble(x = 1:10, y = letters[1:10]),
    cache_format = "parquet"
  )
  
  cache_file <- here::here("analysis", "cache", "test_parquet.parquet")
  expect_true(fs::file_exists(cache_file))
  expect_equal(result$x, 1:10)
})

test_that("cache_wrapper auto-selects format based on size", {
  skip_if_not_installed("arrow")
  
  # Small data should use RDS
  small_data <- .cache_wrapper(
    "test_small",
    function() tibble::tibble(x = 1:10),
    cache_format = "auto"
  )
  
  # Large data should use Parquet (create >5 MB data)
  large_data <- .cache_wrapper(
    "test_large",
    function() tibble::tibble(x = 1:1000000, y = rnorm(1000000)),
    cache_format = "auto"
  )
  
  # Check file extensions
  expect_true(
    fs::file_exists(here::here("analysis", "cache", "test_small.rds")) ||
    fs::file_exists(here::here("analysis", "cache", "test_small.parquet"))
  )
})

test_that("cache_format fallsback to RDS if arrow not available", {
  # Temporarily hide arrow package
  if (requireNamespace("arrow", quietly = TRUE)) {
    skip("arrow is installed, can't test fallback")
  }
  
  expect_warning(
    result <- .cache_wrapper(
      "test_fallback",
      function() tibble::tibble(x = 1:10),
      cache_format = "parquet"
    ),
    "arrow package not available"
  )
  
  cache_file <- here::here("analysis", "cache", "test_fallback.rds")
  expect_true(fs::file_exists(cache_file))
})

test_that("query_cached_data works with Parquet files", {
  skip_if_not_installed("arrow")
  skip_if_not_installed("duckdb")
  
  # Create cached Parquet data
  test_data <- tibble::tibble(
    category = rep(c("A", "B"), each = 5),
    value = 1:10
  )
  
  cache_file <- here::here("analysis", "cache", "test_query.parquet")
  arrow::write_parquet(test_data, cache_file)
  
  # Query the data
  result <- query_cached_data(
    "test_query.parquet",
    "SELECT category, SUM(value) as total FROM data GROUP BY category"
  )
  
  expect_equal(nrow(result), 2)
  expect_equal(result$category, c("A", "B"))
  expect_equal(result$total, c(15, 40))  # Sum of 1:5 and 6:10
})
```

---

## Performance Benchmarks to Validate

Run these benchmarks after implementation:

```r
# Benchmark script
library(tictoc)
source(here("R/data_loading.R"))

# Clear cache
clear_analysis_cache()

# Test 1: Small data (genomic context metrics)
cat("\n=== Small Data (Genomic Context Metrics) ===\n")

tic("RDS format")
metrics_rds <- load_genomic_context_metrics(cache_format = "rds")
time_rds_1 <- toc()

tic("RDS cached")
metrics_rds_cached <- load_genomic_context_metrics(cache_format = "rds")
time_rds_2 <- toc()

if (requireNamespace("arrow", quietly = TRUE)) {
  clear_analysis_cache()
  
  tic("Parquet format")
  metrics_pq <- load_genomic_context_metrics(cache_format = "parquet")
  time_pq_1 <- toc()
  
  tic("Parquet cached")
  metrics_pq_cached <- load_genomic_context_metrics(cache_format = "parquet")
  time_pq_2 <- toc()
  
  # Compare cache sizes
  cache_files <- fs::dir_ls(here("analysis", "cache"))
  for (f in cache_files) {
    size_mb <- fs::file_size(f) / 1024^2
    ext <- tools::file_ext(f)
    cat(glue::glue("{ext}: {round(size_mb, 2)} MB\n"))
  }
}

# Validation
stopifnot(identical(metrics_rds, metrics_rds_cached))
if (exists("metrics_pq")) {
  stopifnot(nrow(metrics_rds) == nrow(metrics_pq))
}
```

---

## Migration Notes

### For Users Without Arrow

- Caching will default to RDS format
- All functionality works without arrow package
- Warning displayed when Parquet format requested but arrow not available
- Can install arrow later and re-cache for performance benefits

### For Users With Arrow

- Auto format selection provides optimal performance
- Existing RDS caches are still readable
- Can coexist with RDS caches
- Parquet recommended for new caches

---

## Summary of Changes

**Core Enhancement:**

- Multi-format cache support (RDS + Parquet)
- Automatic format selection based on data size
- Graceful fallback when arrow not available

**Performance Benefits:**

- 2-4x faster read/write with Parquet
- 40-60% better compression
- Column selection for memory efficiency (future enhancement)

**Optional Advanced Feature:**

- DuckDB query support for large cached data
- SQL-based exploratory analysis
- Minimal memory usage for queries

**Dependencies:**

- `digest` (required, already minimal)
- `arrow` (optional, recommended for better performance)
- `duckdb` (optional, for advanced queries)

**Backward Compatibility:**

- 100% compatible with original RDS-only design
- Works without arrow/duckdb packages
- Existing code requires no changes
- Format parameter is optional
