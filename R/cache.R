#' Parquet Caching Infrastructure
#'
#' Parquet-based caching for large datasets using Arrow. Cache invalidation
#' is based on source file modification times. Pipeline metadata is embedded
#' in Parquet file-level key-value metadata.
#'
#' Cache directory: analysis/cache/
#' File naming: {dataset_name}_{hash}.parquet
#' Compression: zstd, level 3

# Source schema registry
source(here::here("R/schemas.R"))

# ---- Constants --------------------------------------------------------------
CACHE_COMPRESSION <- "zstd"
CACHE_COMPRESSION_LEVEL <- 3L

#' Get Cache Directory Path
#'
#' Returns the cache directory, configurable via the `q100.cache_dir` option.
#' Defaults to `analysis/cache/` in the project root.
#'
#' @return Character string path to cache directory
#'
#' @keywords internal
.get_cache_dir <- function() {
  getOption("q100.cache_dir", default = here::here("analysis/cache"))
}


# ---- Internal helpers -------------------------------------------------------

#' Build Deterministic Cache Key
#'
#' Creates a hash from dataset name, source file modification times, and
#' parameters. When source files change (pipeline re-run), the mtime changes,
#' producing a new hash and a cache miss.
#'
#' @param dataset_name Character string identifying the dataset
#' @param source_files Character vector of source file paths
#' @param params Named list of parameters (filters, benchmark_id, etc.)
#'
#' @return Character string hash
#'
#' @keywords internal
.build_cache_key <- function(dataset_name, source_files, params = list()) {
  # Get modification times for all source files
  mtimes <- file.info(source_files)$mtime
  mtimes_str <- as.character(mtimes)

  # Build deterministic key components
  key_data <- list(
    dataset = dataset_name,
    mtimes = mtimes_str,
    params = params
  )

  rlang::hash(key_data)
}


#' Find Cache File for a Given Key
#'
#' @param dataset_name Character string identifying the dataset
#' @param cache_key Character string hash
#'
#' @return File path if cache exists, NULL otherwise
#'
#' @keywords internal
.find_cache_file <- function(dataset_name, cache_key) {
  cache_path <- fs::path(.get_cache_dir(), glue::glue("{dataset_name}_{cache_key}.parquet"))

  if (fs::file_exists(cache_path)) {
    return(cache_path)
  }

  NULL
}


#' Strip Factors to Character
#'
#' Converts factor columns to character before Parquet write, since Parquet
#' stores strings natively. Factor levels are restored on read via the
#' schema registry.
#'
#' @param df Data frame
#'
#' @return Data frame with factors converted to character
#'
#' @keywords internal
.strip_factors <- function(df) {
  factor_cols <- vapply(df, is.factor, logical(1))
  if (any(factor_cols)) {
    df[factor_cols] <- lapply(df[factor_cols], as.character)
  }
  df
}


#' Restore Factor Levels from Schema Registry
#'
#' Converts character columns back to factors using levels defined in
#' the schema registry.
#'
#' @param df Data frame
#' @param dataset_name Character string identifying the dataset
#'
#' @return Data frame with factor columns restored
#'
#' @keywords internal
.restore_factors <- function(df, dataset_name) {
  factor_levels <- get_factor_levels(dataset_name)

  for (col_name in names(factor_levels)) {
    if (col_name %in% names(df)) {
      df[[col_name]] <- factor(df[[col_name]], levels = factor_levels[[col_name]])
    }
  }

  df
}


# ---- Pipeline metadata ------------------------------------------------------

#' Collect Pipeline Metadata
#'
#' Gathers pipeline configuration, R version, and package versions for
#' embedding in cached Parquet files.
#'
#' @return Named list with cache_date, r_version, package_versions, and
#'   pipeline_config_summary
#'
#' @export
collect_pipeline_metadata <- function() {
  # Package versions for key dependencies
  pkg_names <- c("arrow", "dplyr", "vroom", "tidyverse", "here", "glue")
  pkg_versions <- vapply(pkg_names, function(pkg) {
    tryCatch(
      as.character(utils::packageVersion(pkg)),
      error = function(e) "not installed"
    )
  }, character(1))
  names(pkg_versions) <- pkg_names

  # Pipeline config summary (best-effort)
  config_summary <- tryCatch(
    {
      config <- parse_pipeline_config()
      list(
        benchmarks = config$benchmarks,
        references = config$references
      )
    },
    error = function(e) {
      list(note = "config.yaml not available")
    }
  )

  list(
    cache_date = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    r_version = R.version.string,
    package_versions = pkg_versions,
    pipeline_config_summary = config_summary
  )
}


# ---- Public cache functions --------------------------------------------------

#' Write Data to Parquet Cache
#'
#' Validates the data, strips factors, embeds pipeline metadata as JSON in
#' Parquet file-level key-value metadata, and writes to the cache directory.
#'
#' @param data Data frame to cache
#' @param dataset_name Character string identifying the dataset
#' @param source_files Character vector of source file paths used to generate
#'   the data (for cache invalidation)
#' @param params Named list of parameters (filters, benchmark_id, etc.)
#'
#' @return Invisible path to the cached file
#'
#' @export
write_cache <- function(data, dataset_name, source_files, params = list()) {
  # Validate data against schema

  validate_data(data, dataset_name)

  # Build cache key and path
  cache_key <- .build_cache_key(dataset_name, source_files, params)
  cache_dir <- .get_cache_dir()
  cache_path <- fs::path(cache_dir, glue::glue("{dataset_name}_{cache_key}.parquet"))

  # Ensure cache directory exists
  fs::dir_create(cache_dir, recurse = TRUE)

  # Strip factors for Parquet storage
  data_clean <- .strip_factors(data)

  # Get Arrow schema for type enforcement
  schema <- get_arrow_schema(dataset_name)

  # Select only columns in the schema (data may have extra columns)
  schema_cols <- names(schema)
  data_cols <- names(data_clean)

  # Use schema columns that exist in the data
  common_cols <- intersect(schema_cols, data_cols)
  extra_cols <- setdiff(data_cols, schema_cols)

  # Build Arrow table with schema for typed columns, keep extras as-is
  if (length(extra_cols) > 0) {
    # Build partial schema for columns that match
    partial_fields <- lapply(common_cols, function(col) schema$GetFieldByName(col))
    partial_schema <- arrow::schema(
      !!!stats::setNames(
        lapply(partial_fields, function(f) f$type),
        common_cols
      )
    )
    # Can't enforce full schema with extra columns; write as-is
    tbl <- arrow::arrow_table(data_clean)
  } else {
    tbl <- arrow::arrow_table(data_clean, schema = schema)
  }

  # Collect and embed pipeline metadata
  metadata <- collect_pipeline_metadata()
  metadata_json <- jsonlite::toJSON(metadata, auto_unbox = TRUE, pretty = FALSE)

  # Set file-level metadata
  existing_meta <- tbl$metadata
  if (is.null(existing_meta)) existing_meta <- list()
  existing_meta[["q100_pipeline_metadata"]] <- as.character(metadata_json)
  tbl$metadata <- existing_meta

  # Write Parquet file
  arrow::write_parquet(
    tbl,
    sink = cache_path,
    compression = CACHE_COMPRESSION,
    compression_level = CACHE_COMPRESSION_LEVEL
  )

  message(glue::glue("Cache written: {cache_path}"))
  invisible(cache_path)
}


#' Read Data from Parquet Cache
#'
#' Reads a cached Parquet file, restores factor levels from the schema
#' registry, and returns a tibble. Returns NULL if no valid cache exists.
#'
#' @param dataset_name Character string identifying the dataset
#' @param source_files Character vector of source file paths
#' @param params Named list of parameters
#' @param validate Logical; if TRUE, validate data after reading (default FALSE)
#'
#' @return A tibble if cache hit, NULL if cache miss
#'
#' @export
read_cache <- function(dataset_name, source_files, params = list(),
                       validate = FALSE) {
  cache_key <- .build_cache_key(dataset_name, source_files, params)
  cache_path <- .find_cache_file(dataset_name, cache_key)

  if (is.null(cache_path)) {
    return(NULL)
  }

  # Read Parquet file
  df <- arrow::read_parquet(cache_path, as_data_frame = TRUE)
  df <- tibble::as_tibble(df)

  # Restore factor levels
  df <- .restore_factors(df, dataset_name)

  # Optional validation on read

  if (validate) {
    validate_data(df, dataset_name)
  }

  message(glue::glue("Cache hit: {cache_path}"))
  df
}


#' Check if Valid Cache Exists
#'
#' @param dataset_name Character string identifying the dataset
#' @param source_files Character vector of source file paths
#' @param params Named list of parameters
#'
#' @return Logical
#'
#' @export
cache_is_valid <- function(dataset_name, source_files, params = list()) {
  cache_key <- .build_cache_key(dataset_name, source_files, params)
  !is.null(.find_cache_file(dataset_name, cache_key))
}


#' Invalidate Cache for a Specific Dataset
#'
#' Removes all cached Parquet files for the given dataset name.
#'
#' @param dataset_name Character string identifying the dataset
#'
#' @return Invisible number of files removed
#'
#' @export
invalidate_cache <- function(dataset_name) {
  cache_dir <- .get_cache_dir()
  if (!fs::dir_exists(cache_dir)) {
    return(invisible(0L))
  }

  pattern <- paste0("^", dataset_name, "_[a-f0-9]+\\.parquet$")
  all_files <- fs::dir_ls(cache_dir, glob = "*.parquet", fail = FALSE)
  files <- all_files[grepl(pattern, fs::path_file(all_files))]

  if (length(files) > 0) {
    fs::file_delete(files)
    message(glue::glue("Removed {length(files)} cache file(s) for '{dataset_name}'"))
  }

  invisible(length(files))
}


#' Clear All Cache Files
#'
#' Removes the entire cache directory and all cached Parquet files.
#'
#' @return Invisible TRUE
#'
#' @export
clear_cache <- function() {
  cache_dir <- .get_cache_dir()
  if (fs::dir_exists(cache_dir)) {
    fs::dir_delete(cache_dir)
    message(glue::glue("Cache directory removed: {cache_dir}"))
  } else {
    message("No cache directory to remove.")
  }

  invisible(TRUE)
}


#' Get Cache Information
#'
#' Lists all cached files with their size, age, and associated dataset name.
#'
#' @return Tibble with columns: file, dataset, size_mb, created, age_hours
#'
#' @export
cache_info <- function() {
  cache_dir <- .get_cache_dir()
  if (!fs::dir_exists(cache_dir)) {
    return(
      tibble::tibble(
        file = character(0),
        dataset = character(0),
        size_mb = numeric(0),
        created = as.POSIXct(character(0)),
        age_hours = numeric(0)
      )
    )
  }

  files <- fs::dir_ls(cache_dir, glob = "*.parquet", fail = FALSE)

  if (length(files) == 0) {
    return(
      tibble::tibble(
        file = character(0),
        dataset = character(0),
        size_mb = numeric(0),
        created = as.POSIXct(character(0)),
        age_hours = numeric(0)
      )
    )
  }

  file_info <- fs::file_info(files)

  tibble::tibble(
    file = fs::path_file(files),
    dataset = stringr::str_extract(fs::path_file(files), "^[^_]+(?=_[a-f0-9]+\\.parquet$)"),
    size_mb = round(as.numeric(file_info$size) / 1024^2, 2),
    created = file_info$modification_time,
    age_hours = round(
      as.numeric(difftime(Sys.time(), file_info$modification_time, units = "hours")),
      1
    )
  )
}
