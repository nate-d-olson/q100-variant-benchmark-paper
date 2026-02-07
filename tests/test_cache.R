library(testthat)
library(here)

# Source the modules under test
source(here::here("R/schemas.R"))
source(here::here("R/data_loading.R"))
source(here::here("R/cache.R"))


# ---- Schema tests -----------------------------------------------------------

test_that("get_arrow_schema returns valid Schema for each dataset", {
  for (ds in c("variant_table", "diff_coverage", "benchmark_regions")) {
    schema <- get_arrow_schema(ds)
    expect_s3_class(schema, "Schema")
    expect_gt(length(schema), 0)
  }
})

test_that("get_arrow_schema errors for unknown datasets", {
  expect_error(
    get_arrow_schema("nonexistent_dataset"),
    "Unknown dataset"
  )
})

test_that("get_factor_levels returns correct levels", {
  levels <- get_factor_levels("benchmark_regions")
  expect_type(levels, "list")
  expect_true("bench_version" %in% names(levels))
  expect_equal(levels$bench_version, c("v0.6", "v4.2.1", "v5.0q"))
  expect_equal(levels$ref, c("GRCh37", "GRCh38", "CHM13v2.0"))
  expect_equal(levels$chrom, paste0("chr", c(1:22, "X", "Y")))
})

test_that("get_factor_levels errors for unknown datasets", {
  expect_error(
    get_factor_levels("nonexistent_dataset"),
    "Unknown dataset"
  )
})

test_that("validate_data passes valid data", {
  valid_df <- tibble::tibble(
    context_name = c("HP", "TR"),
    chrom = c("chr1", "chr2"),
    start = c(100L, 200L),
    end = c(200L, 300L),
    n_overlap = c(1L, 2L),
    bases_cov = c(50L, 75L),
    ivl_len = c(100L, 100L),
    frac_cov = c(0.5, 0.75)
  )

  expect_true(validate_data(valid_df, "diff_coverage"))
})

test_that("validate_data fails on missing columns", {
  incomplete_df <- tibble::tibble(
    context_name = "HP",
    chrom = "chr1"
  )

  expect_error(
    validate_data(incomplete_df, "diff_coverage"),
    "missing required columns"
  )
})

test_that("validate_data fails on rule violations", {
  # frac_cov out of [0, 1] range
  bad_df <- tibble::tibble(
    context_name = c("HP"),
    chrom = c("chr1"),
    start = c(100L),
    end = c(200L),
    n_overlap = c(1L),
    bases_cov = c(50L),
    ivl_len = c(100L),
    frac_cov = c(1.5)
  )

  expect_error(
    validate_data(bad_df, "diff_coverage"),
    "Validation failed for column 'frac_cov'"
  )
})

test_that("validate_data fails on negative interval_size", {
  bad_df <- tibble::tibble(
    bench_version = "v5.0q",
    ref = "GRCh38",
    bench_type = "smvar",
    chrom = "chr1",
    start = 100L,
    end = 200L,
    interval_size = -10L
  )

  expect_error(
    validate_data(bad_df, "benchmark_regions"),
    "Validation failed for column 'interval_size'"
  )
})


# ---- Cache round-trip tests -------------------------------------------------

test_that("write_cache + read_cache preserves data and factor levels", {
  cache_dir <- withr::local_tempdir()
  withr::local_options(q100.cache_dir = cache_dir)

  # Create synthetic data with factors
  test_data <- tibble::tibble(
    context_name = factor(c("HP", "TR", "SD"), levels = CONTEXT_NAME_LEVELS),
    chrom = c("chr1", "chr2", "chr3"),
    start = c(100L, 200L, 300L),
    end = c(200L, 300L, 400L),
    n_overlap = c(1L, 2L, 3L),
    bases_cov = c(50L, 75L, 80L),
    ivl_len = c(100L, 100L, 100L),
    frac_cov = c(0.5, 0.75, 0.8)
  )

  # Create a temporary source file to hash against
  source_file <- fs::path(cache_dir, "test_source.bed")
  writeLines("test content", source_file)

  # Write cache
  cache_path <- write_cache(
    test_data,
    "diff_coverage",
    source_files = source_file,
    params = list(benchmark_id = "test")
  )

  expect_true(fs::file_exists(cache_path))

  # Read cache
  result <- read_cache(
    "diff_coverage",
    source_files = source_file,
    params = list(benchmark_id = "test")
  )

  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 3)
  expect_equal(as.character(result$context_name), c("HP", "TR", "SD"))
  expect_true(is.factor(result$context_name))
  expect_equal(levels(result$context_name), CONTEXT_NAME_LEVELS)
  expect_equal(result$frac_cov, c(0.5, 0.75, 0.8))
})

test_that("cache_is_valid returns FALSE when no cache exists", {
  cache_dir <- withr::local_tempdir()
  withr::local_options(q100.cache_dir = cache_dir)

  source_file <- fs::path(cache_dir, "test_source.bed")
  writeLines("test", source_file)

  expect_false(cache_is_valid("diff_coverage", source_files = source_file))
})

test_that("invalidate_cache removes files for a specific dataset", {
  cache_dir <- withr::local_tempdir()
  withr::local_options(q100.cache_dir = cache_dir)

  # Create test data
  test_data <- tibble::tibble(
    context_name = c("HP"),
    chrom = c("chr1"),
    start = 100L,
    end = 200L,
    n_overlap = 1L,
    bases_cov = 50L,
    ivl_len = 100L,
    frac_cov = 0.5
  )

  source_file <- fs::path(cache_dir, "test_source.bed")
  writeLines("test", source_file)

  # Write cache
  write_cache(test_data, "diff_coverage", source_files = source_file)

  # Verify file exists
  parquet_files <- fs::dir_ls(cache_dir, glob = "*.parquet")
  expect_equal(length(parquet_files), 1)

  # Invalidate
  invalidate_cache("diff_coverage")

  # Verify file removed
  parquet_files <- fs::dir_ls(cache_dir, glob = "*.parquet", fail = FALSE)
  expect_equal(length(parquet_files), 0)
})

test_that("clear_cache removes all files", {
  cache_dir <- withr::local_tempdir()
  withr::local_options(q100.cache_dir = cache_dir)
  fs::dir_create(cache_dir)

  # Create a dummy file
  writeLines("test", fs::path(cache_dir, "test.parquet"))

  clear_cache()

  expect_false(fs::dir_exists(cache_dir))
})

test_that("cache_info returns correct tibble", {
  cache_dir <- withr::local_tempdir()
  withr::local_options(q100.cache_dir = cache_dir)

  test_data <- tibble::tibble(
    context_name = c("HP"),
    chrom = c("chr1"),
    start = 100L,
    end = 200L,
    n_overlap = 1L,
    bases_cov = 50L,
    ivl_len = 100L,
    frac_cov = 0.5
  )

  source_file <- fs::path(cache_dir, "test_source.bed")
  writeLines("test", source_file)

  write_cache(test_data, "diff_coverage", source_files = source_file)

  info <- cache_info()

  expect_s3_class(info, "tbl_df")
  expect_equal(nrow(info), 1)
  expect_true("file" %in% names(info))
  expect_true("dataset" %in% names(info))
  expect_true("size_mb" %in% names(info))
  expect_true("created" %in% names(info))
  expect_true("age_hours" %in% names(info))
})


# ---- Metadata tests ---------------------------------------------------------

test_that("collect_pipeline_metadata returns expected structure", {
  meta <- collect_pipeline_metadata()

  expect_type(meta, "list")
  expect_true("cache_date" %in% names(meta))
  expect_true("r_version" %in% names(meta))
  expect_true("package_versions" %in% names(meta))
  expect_true("pipeline_config_summary" %in% names(meta))

  expect_type(meta$r_version, "character")
  expect_type(meta$package_versions, "character")
})

test_that("written Parquet files contain q100_pipeline_metadata", {
  cache_dir <- withr::local_tempdir()
  withr::local_options(q100.cache_dir = cache_dir)

  test_data <- tibble::tibble(
    context_name = c("HP"),
    chrom = c("chr1"),
    start = 100L,
    end = 200L,
    n_overlap = 1L,
    bases_cov = 50L,
    ivl_len = 100L,
    frac_cov = 0.5
  )

  source_file <- fs::path(cache_dir, "test_source.bed")
  writeLines("test", source_file)

  cache_path <- write_cache(
    test_data,
    "diff_coverage",
    source_files = source_file
  )

  # Read file metadata
  pq_file <- arrow::read_parquet(cache_path, as_data_frame = FALSE)
  file_meta <- pq_file$metadata

  expect_true("q100_pipeline_metadata" %in% names(file_meta))

  # Parse JSON and verify structure
  parsed_meta <- jsonlite::fromJSON(file_meta[["q100_pipeline_metadata"]])
  expect_true("cache_date" %in% names(parsed_meta))
  expect_true("r_version" %in% names(parsed_meta))
})
