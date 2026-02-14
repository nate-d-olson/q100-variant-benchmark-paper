#' Data Schema Registry
#'
#' Centralized Arrow schema definitions, factor levels, and validation rules
#' for cached datasets. Adding a new dataset requires entries in
#' get_arrow_schema(), get_factor_levels(), and get_validation_rules().

# ---- Factor level constants ------------------------------------------------
BENCH_VERSION_LEVELS <- c("v0.6", "v4.2.1", "v5.0q")
REF_LEVELS <- c("GRCh37", "GRCh38", "CHM13v2.0")
BENCH_TYPE_LEVELS <- c("smvar", "stvar")
CONTEXT_NAME_LEVELS <- c("HP", "MAP", "SD", "SD10kb", "TR", "TR10kb")
CHROM_LEVELS <- paste0("chr", c(1:22, "X", "Y"))
AUTOSOME_CHROM_LEVELS <- paste0("chr", 1:22)


#' Get Arrow Schema for a Dataset
#'
#' Returns an arrow::Schema used for type enforcement at Parquet write time.
#'
#' @param dataset_name Character string identifying the dataset
#'
#' @return An arrow::Schema object
#'
#' @examples
#' \dontrun{
#' schema <- get_arrow_schema("variant_table")
#' }
#'
#' @export
get_arrow_schema <- function(dataset_name) {
  schemas <- list(
    variant_table = arrow::schema(
      bench_version = arrow::utf8(),
      ref = arrow::utf8(),
      bench_type = arrow::utf8(),
      chrom = arrow::utf8(),
      pos = arrow::int64(),
      end = arrow::int64(),
      gt = arrow::utf8(),
      var_type = arrow::utf8(),
      var_size = arrow::int64(),
      szbin = arrow::utf8(),
      ref_len = arrow::int32(),
      alt_len = arrow::int32(),
      qual = arrow::float64(),
      filter = arrow::utf8(),
      is_pass = arrow::boolean(),
      ## Genomic context boolean columns (one per context)
      HP = arrow::boolean(),
      MAP = arrow::boolean(),
      SD = arrow::boolean(),
      SD10kb = arrow::boolean(),
      TR = arrow::boolean(),
      TR10kb = arrow::boolean(),
      ## Region membership
      in_benchmark = arrow::boolean()
      ## Note: excl_* exclusion columns are dynamic per benchmark and
      ## not included in the fixed schema. They are present in v5.0q
      ## Parquet files but absent from v4.2.1 and v0.6 files.
    ),
    diff_coverage = arrow::schema(
      context_name = arrow::utf8(),
      chrom = arrow::utf8(),
      start = arrow::int64(),
      end = arrow::int64(),
      n_overlap = arrow::int32(),
      bases_cov = arrow::int32(),
      ivl_len = arrow::int32(),
      frac_cov = arrow::float64()
    ),
    benchmark_regions = arrow::schema(
      bench_version = arrow::utf8(),
      ref = arrow::utf8(),
      bench_type = arrow::utf8(),
      chrom = arrow::utf8(),
      start = arrow::int64(),
      end = arrow::int64(),
      interval_size = arrow::int64()
    ),
    platinum_pedigree_regions = arrow::schema(
      bench_version = arrow::utf8(),
      ref = arrow::utf8(),
      bench_type = arrow::utf8(),
      chrom = arrow::utf8(),
      start = arrow::int64(),
      end = arrow::int64(),
      interval_size = arrow::int64()
    )
  )

  if (!dataset_name %in% names(schemas)) {
    stop(
      glue::glue(
        "Unknown dataset: '{dataset_name}'. ",
        "Available datasets: {paste(names(schemas), collapse = ', ')}"
      ),
      call. = FALSE
    )
  }

  schemas[[dataset_name]]
}


#' Get Factor Levels for a Dataset
#'
#' Returns a named list of factor levels for columns that should be converted
#' back to factors after reading from Parquet (which stores them as strings).
#'
#' @param dataset_name Character string identifying the dataset
#'
#' @return Named list where names are column names and values are character
#'   vectors of factor levels
#'
#' @examples
#' \dontrun{
#' levels <- get_factor_levels("benchmark_regions")
#' }
#'
#' @export
get_factor_levels <- function(dataset_name) {
  factor_defs <- list(
    variant_table = list(
      bench_version = BENCH_VERSION_LEVELS,
      ref = REF_LEVELS,
      bench_type = BENCH_TYPE_LEVELS,
      chrom = CHROM_LEVELS
    ),
    diff_coverage = list(
      context_name = CONTEXT_NAME_LEVELS
    ),
    benchmark_regions = list(
      bench_version = BENCH_VERSION_LEVELS,
      ref = REF_LEVELS,
      bench_type = BENCH_TYPE_LEVELS,
      chrom = CHROM_LEVELS
    ),
    platinum_pedigree_regions = list(
      bench_version = c("PP"),
      ref = c("GRCh38"),
      bench_type = BENCH_TYPE_LEVELS,
      chrom = CHROM_LEVELS
    )
  )

  if (!dataset_name %in% names(factor_defs)) {
    stop(
      glue::glue(
        "Unknown dataset: '{dataset_name}'. ",
        "Available datasets: {paste(names(factor_defs), collapse = ', ')}"
      ),
      call. = FALSE
    )
  }

  factor_defs[[dataset_name]]
}


#' Get Validation Rules for a Dataset
#'
#' Returns a named list of predicate functions for dataset-specific validation.
#' Each predicate takes a column vector and returns TRUE if valid.
#'
#' @param dataset_name Character string identifying the dataset
#'
#' @return Named list of predicate functions keyed by column name
#'
#' @keywords internal
get_validation_rules <- function(dataset_name) {
  rules <- list(
    variant_table = list(
      pos = function(x) all(x >= 0, na.rm = TRUE),
      end = function(x) all(x >= 0, na.rm = TRUE),
      var_size = function(x) all(!is.na(x))
    ),
    diff_coverage = list(
      frac_cov = function(x) all(x >= 0 & x <= 1, na.rm = TRUE),
      ivl_len = function(x) all(x > 0, na.rm = TRUE),
      bases_cov = function(x) all(x >= 0, na.rm = TRUE)
    ),
    benchmark_regions = list(
      interval_size = function(x) all(x > 0, na.rm = TRUE),
      start = function(x) all(x >= 0, na.rm = TRUE)
    ),
    platinum_pedigree_regions = list(
      interval_size = function(x) all(x > 0, na.rm = TRUE),
      start = function(x) all(x >= 0, na.rm = TRUE)
    )
  )

  if (!dataset_name %in% names(rules)) {
    stop(
      glue::glue(
        "Unknown dataset: '{dataset_name}'. ",
        "Available datasets: {paste(names(rules), collapse = ', ')}"
      ),
      call. = FALSE
    )
  }

  rules[[dataset_name]]
}


#' Validate Data Against Schema
#'
#' Checks that required columns are present and runs dataset-specific
#' validation rules. Called at cache write time for fail-fast behavior.
#'
#' @param df Data frame to validate
#' @param dataset_name Character string identifying the dataset
#'
#' @return Invisible TRUE if validation passes; raises an error otherwise
#'
#' @examples
#' \dontrun{
#' validate_data(my_data, "diff_coverage")
#' }
#'
#' @export
validate_data <- function(df, dataset_name) {
  schema <- get_arrow_schema(dataset_name)
  rules <- get_validation_rules(dataset_name)

  # Check required columns
  schema_cols <- names(schema)
  df_cols <- names(df)
  missing_cols <- setdiff(schema_cols, df_cols)

  if (length(missing_cols) > 0) {
    stop(
      glue::glue(
        "Dataset '{dataset_name}' is missing required columns: ",
        "{paste(missing_cols, collapse = ', ')}"
      ),
      call. = FALSE
    )
  }

  # Run validation rules
  for (col_name in names(rules)) {
    if (col_name %in% df_cols) {
      rule_fn <- rules[[col_name]]
      if (!rule_fn(df[[col_name]])) {
        stop(
          glue::glue(
            "Validation failed for column '{col_name}' in dataset '{dataset_name}'"
          ),
          call. = FALSE
        )
      }
    }
  }

  invisible(TRUE)
}
