# Proposal 2: Standardized Data Access Layer with Validation Framework

## Overview

Create a comprehensive data access layer that standardizes how Quarto notebooks interact with pipeline outputs. This proposal introduces a unified interface with built-in validation, automatic schema checking, and consistent error handling across all data loading operations.

## Key Principles

- **Single source of truth** - One canonical way to access each data type
- **Fail fast** - Validate data immediately on load
- **Self-documenting** - Schema definitions serve as documentation
- **Type safety** - Explicit column types and value ranges

## Technical Approach

### 1. Data Schema Registry

Define explicit schemas for all pipeline outputs in a central registry:

```r
# R/schemas.R (NEW FILE)

#' Data Schema Registry
#' 
#' Defines expected structure for all pipeline outputs
PIPELINE_SCHEMAS <- list(
  
  genomic_context_metrics = list(
    required_columns = c(
      "bench_version", "ref", "var_type", "context_name",
      "context_bp", "intersect_bp", "pct_of_context", "pct_of_bench",
      "total_variants", "snp_count", "indel_count", 
      "del_count", "ins_count", "complex_count", "other_count",
      "variant_density_per_mb"
    ),
    column_types = list(
      bench_version = "character",
      ref = "character",
      var_type = "character",
      context_name = "character",
      context_bp = "numeric",
      intersect_bp = "numeric",
      pct_of_context = "numeric",
      pct_of_bench = "numeric",
      total_variants = "integer",
      snp_count = "integer",
      indel_count = "integer",
      del_count = "integer",
      ins_count = "integer",
      complex_count = "integer",
      other_count = "integer",
      variant_density_per_mb = "numeric"
    ),
    factor_levels = list(
      bench_version = c("v0.6", "v4.2.1", "v5.0q"),
      ref = c("GRCh37", "GRCh38", "CHM13v2.0"),
      var_type = c("smvar", "stvar"),
      context_name = c("HP", "MAP", "SD", "SD10kb", "TR", "TR10kb")
    ),
    validation_rules = list(
      # Custom validation functions
      context_bp = function(x) all(x > 0),
      pct_of_context = function(x) all(x >= 0 & x <= 100),
      pct_of_bench = function(x) all(x >= 0 & x <= 100),
      total_variants = function(x) all(x >= 0)
    ),
    description = "Per-genomic-context metrics with variant counts"
  ),
  
  variant_table = list(
    required_columns = c(
      "chrom", "pos", "end", "gt", "vkx", "var_type",
      "len_ref", "len_alt", "var_size", "region_ids"
    ),
    column_types = list(
      chrom = "character",
      pos = "integer",
      end = "integer",
      gt = "character",
      vkx = "character",
      var_type = "character",
      len_ref = "integer",
      len_alt = "integer",
      var_size = "integer",
      region_ids = "character"
    ),
    validation_rules = list(
      pos = function(x) all(x > 0),
      end = function(x) all(x > 0),
      var_size = function(x) TRUE  # Can be negative for deletions
    ),
    description = "Full variant-level data"
  ),
  
  benchmark_regions = list(
    required_columns = c("bench_version", "ref", "bench_type", 
                        "chrom", "start", "end", "interval_size"),
    column_types = list(
      bench_version = "character",
      ref = "character",
      bench_type = "character",
      chrom = "character",
      start = "integer",
      end = "integer",
      interval_size = "integer"
    ),
    validation_rules = list(
      start = function(x) all(x >= 0),
      end = function(x) all(x > 0),
      interval_size = function(x) all(x > 0)
    ),
    description = "Benchmark region BED intervals"
  )
  
  # Additional schemas for other data types...
)
```

### 2. Validation Engine

Create a validation engine that checks data against schemas:

```r
# R/validation.R (NEW FILE)

#' Validate DataFrame Against Schema
#' 
#' @param df Data frame to validate
#' @param schema_name Name of schema in PIPELINE_SCHEMAS
#' @param strict If TRUE, fail on extra columns; if FALSE, warn only
#' 
#' @return Validated data frame (throws error if validation fails)
validate_data <- function(df, schema_name, strict = FALSE) {
  schema <- PIPELINE_SCHEMAS[[schema_name]]
  
  if (is.null(schema)) {
    stop(glue::glue("Unknown schema: {schema_name}"), call. = FALSE)
  }
  
  errors <- c()
  warnings <- c()
  
  # Check required columns
  missing_cols <- setdiff(schema$required_columns, colnames(df))
  if (length(missing_cols) > 0) {
    errors <- c(errors, glue::glue(
      "Missing required columns: {paste(missing_cols, collapse = ', ')}"
    ))
  }
  
  # Check for extra columns
  extra_cols <- setdiff(colnames(df), schema$required_columns)
  if (length(extra_cols) > 0) {
    msg <- glue::glue(
      "Extra columns found: {paste(extra_cols, collapse = ', ')}"
    )
    if (strict) {
      errors <- c(errors, msg)
    } else {
      warnings <- c(warnings, msg)
    }
  }
  
  # Validate column types
  for (col in names(schema$column_types)) {
    if (col %in% colnames(df)) {
      expected_type <- schema$column_types[[col]]
      actual_type <- typeof(df[[col]])
      
      if (!types_compatible(actual_type, expected_type)) {
        errors <- c(errors, glue::glue(
          "Column '{col}' has type '{actual_type}' but expected '{expected_type}'"
        ))
      }
    }
  }
  
  # Run custom validation rules
  for (col in names(schema$validation_rules)) {
    if (col %in% colnames(df)) {
      validator <- schema$validation_rules[[col]]
      if (!validator(df[[col]])) {
        errors <- c(errors, glue::glue(
          "Validation failed for column '{col}'"
        ))
      }
    }
  }
  
  # Report results
  if (length(warnings) > 0) {
    warning(paste(warnings, collapse = "\n"), call. = FALSE)
  }
  
  if (length(errors) > 0) {
    stop(
      glue::glue(
        "Data validation failed for schema '{schema_name}':\n",
        "{paste(errors, collapse = '\n')}"
      ),
      call. = FALSE
    )
  }
  
  message(glue::glue("✓ Data validated against schema: {schema_name}"))
  return(df)
}

#' Check Type Compatibility
types_compatible <- function(actual, expected) {
  type_map <- list(
    character = c("character", "factor"),
    numeric = c("double", "numeric", "integer"),
    integer = c("integer", "numeric")
  )
  
  actual %in% type_map[[expected]]
}
```

### 3. Unified Data Access Interface

Wrap existing loading functions with validation layer:

```r
# R/data_access.R (NEW FILE)

#' Unified Data Access Layer
#' 
#' Provides validated, type-safe access to all pipeline outputs

#' Load Genomic Context Metrics (Validated)
#' 
#' @inheritParams load_genomic_context_metrics
#' @param validate Enable schema validation (default: TRUE)
#' 
#' @return Validated tibble
data_access_genomic_metrics <- function(results_dir = NULL,
                                        benchmark_filter = NULL,
                                        validate = TRUE) {
  # Load using existing function
  df <- load_genomic_context_metrics(
    results_dir = results_dir,
    benchmark_filter = benchmark_filter
  )
  
  # Validate against schema
  if (validate) {
    df <- validate_data(df, "genomic_context_metrics")
  }
  
  return(df)
}

#' Load Variant Table (Validated)
data_access_variants <- function(benchmark_id,
                                  results_dir = NULL,
                                  filters = NULL,
                                  validate = TRUE) {
  # Load using existing function
  df <- load_variant_table(
    benchmark_id = benchmark_id,
    results_dir = results_dir,
    filters = filters
  )
  
  # Validate against schema
  if (validate) {
    df <- validate_data(df, "variant_table")
  }
  
  return(df)
}

#' Load Benchmark Regions (Validated)
data_access_benchmark_regions <- function(resources_dir = NULL,
                                          validate = TRUE) {
  df <- load_benchmark_regions(resources_dir = resources_dir)
  
  if (validate) {
    df <- validate_data(df, "benchmark_regions")
  }
  
  return(df)
}

# Convenience function: Load All Primary Data
data_access_primary <- function(results_dir = NULL,
                                 resources_dir = NULL,
                                 validate = TRUE) {
  list(
    genomic_metrics = data_access_genomic_metrics(
      results_dir = results_dir,
      validate = validate
    ),
    benchmark_regions = data_access_benchmark_regions(
      resources_dir = resources_dir,
      validate = validate
    ),
    ref_sizes = load_reference_sizes(results_dir = results_dir)
  )
}
```

### 4. External Data Management

Add validation for external data sources (Google Sheets, CSV files):

```r
# R/external_data.R (NEW FILE)

#' Load External Evaluation Data
#' 
#' Validates and standardizes external evaluation datasets
#' 
#' @param data_dir Directory containing external CSV files
#' @param validate Enable schema validation
#' 
#' @return List with smvar_evals and stvar_evals tibbles
load_external_evaluations <- function(data_dir = NULL,
                                      validate = TRUE) {
  if (is.null(data_dir)) {
    data_dir <- here::here("data/external-evaluations")
  }
  
  # Find all evaluation CSV files
  eval_csvs <- fs::dir_ls(
    data_dir,
    regexp = "Miqa\\.csv$",
    recurse = TRUE
  )
  
  if (length(eval_csvs) == 0) {
    stop(glue::glue("No evaluation CSV files found in {data_dir}"), 
         call. = FALSE)
  }
  
  # Load and classify by variant type
  all_evals <- eval_csvs %>%
    purrr::map_dfr(
      ~ vroom::vroom(.x, show_col_types = FALSE),
      .id = "callset"
    ) %>%
    dplyr::mutate(
      callset = stringr::str_extract(
        callset,
        "(?<=evaluations/).*(?= On)"
      )
    )
  
  # Split by variant type (based on row count heuristic)
  nrows_by_file <- eval_csvs %>%
    purrr::map_int(~ nrow(vroom::vroom(.x, show_col_types = FALSE)))
  
  smvar_evals <- all_evals %>%
    dplyr::filter(callset %in% names(nrows_by_file)[nrows_by_file > 50])
  
  stvar_evals <- all_evals %>%
    dplyr::filter(callset %in% names(nrows_by_file)[nrows_by_file <= 50])
  
  if (validate) {
    # Add schema validation for external data
    validate_external_evals(smvar_evals, "smvar")
    validate_external_evals(stvar_evals, "stvar")
  }
  
  list(
    smvar_evals = smvar_evals,
    stvar_evals = stvar_evals
  )
}

#' Validate External Evaluation Data
validate_external_evals <- function(df, var_type) {
  required_cols <- c(
    "callset", "GRCh38_chr", "GRCh38_start", "GRCh38_end",
    "correct.in.benchmark", "var_type", "strata", "unique_to"
  )
  
  missing <- setdiff(required_cols, colnames(df))
  if (length(missing) > 0) {
    stop(glue::glue(
      "Missing columns in {var_type} evaluations: {paste(missing, collapse = ', ')}"
    ))
  }
  
  message(glue::glue("✓ External {var_type} data validated"))
}
```

## Implementation Steps

### Phase 1: Schema Definition (Week 1)

1. Create `R/schemas.R` with complete schema registry
2. Document all pipeline output structures
3. Add schema for external data sources

### Phase 2: Validation Engine (Week 2)

1. Implement `R/validation.R` with validation logic
2. Add comprehensive test suite for validation
3. Create validation report generator

### Phase 3: Data Access Layer (Week 3)

1. Create `R/data_access.R` with validated wrappers
2. Create `R/external_data.R` for external data
3. Update documentation with new functions

### Phase 4: Notebook Migration (Week 4)

1. Update `benchmarkset_characterization.qmd` to use new layer
2. Update `external_evaluation.qmd` to use `load_external_evaluations()`
3. Update `benchmark_difficult.qmd`
4. Add migration guide for future notebooks

## Benefits

### Data Quality Assurance

- **Early detection** of pipeline output issues
- **Type safety** prevents silent type coercion bugs
- **Schema documentation** serves as data dictionary
- **Consistent validation** across all notebooks

### Developer Experience

- **Self-documenting** - schemas explain data structure
- **Autocomplete support** - clear function signatures
- **Error messages** pinpoint exact validation failures
- **Standardized interface** - one way to load each data type

### Maintainability

- **Centralized** - schema changes update all notebooks
- **Versioned** - schemas can evolve with pipeline
- **Testable** - validation logic has comprehensive tests
- **Future-proof** - easy to add new data sources

## Drawbacks and Mitigations

### Additional Complexity

- **Risk**: More files and abstraction layers
- **Mitigation**:
  - Keep existing `data_loading.R` functions
  - Make validation optional (`validate = FALSE`)
  - Provide migration guide
  - Document with examples

### Performance Overhead

- **Risk**: Validation adds runtime cost
- **Mitigation**:
  - Validation is fast (milliseconds for most datasets)
  - Can be disabled for production runs
  - Cache validation results

### Migration Effort

- **Risk**: Requires updating existing notebooks
- **Mitigation**:
  - Both old and new interfaces work simultaneously
  - Gradual migration, one notebook at a time
  - Backwards compatibility maintained

## File Structure Changes

```
q100-variant-benchmark-paper/
├── R/
│   ├── data_loading.R          # Existing functions (unchanged)
│   ├── schemas.R               # NEW: Data schema registry
│   ├── validation.R            # NEW: Validation engine
│   ├── data_access.R           # NEW: Validated wrappers
│   └── external_data.R         # NEW: External data loaders
├── analysis/
│   └── *.qmd                   # Update to use new functions
└── tests/
    └── test_validation.R       # NEW: Validation tests
```

## Example Usage

### Before (Current Approach)

```r
# benchmarkset_characterization.qmd
source(here("R/data_loading.R"))

metrics <- load_genomic_context_metrics()
regions <- load_benchmark_regions()

# Hope the data is correct...
```

### After (Standardized Layer)

```r
# benchmarkset_characterization.qmd
source(here("R/data_access.R"))

# Load with automatic validation
metrics <- data_access_genomic_metrics()
# ✓ Data validated against schema: genomic_context_metrics

regions <- data_access_benchmark_regions()
# ✓ Data validated against schema: benchmark_regions

# Or load all primary data at once
primary_data <- data_access_primary()
# ✓ Data validated against schema: genomic_context_metrics
# ✓ Data validated against schema: benchmark_regions
```

### Validation Error Example

```r
# If pipeline output is corrupted or schema changes:
metrics <- data_access_genomic_metrics()
# Error: Data validation failed for schema 'genomic_context_metrics':
# - Missing required columns: variant_density_per_mb
# - Column 'total_variants' has type 'character' but expected 'integer'
```

## Comparison to Other Proposals

**Advantages over Proposal 1 (Caching):**

- Guarantees data quality
- Self-documenting schemas
- Catches pipeline bugs early

**Advantages over Proposal 3 (Snapshots):**

- Real-time validation on current data
- No manual snapshot generation
- Works with dynamic queries

**Disadvantages compared to Proposal 1:**

- Higher implementation effort
- More files to maintain
- Requires notebook migration

## Recommendation

**Best suited for:** Teams prioritizing data quality and long-term maintainability. Ideal for production pipelines where data validation is critical.

**Priority Level:** MEDIUM-HIGH - Provides strong foundation but requires upfront investment.
