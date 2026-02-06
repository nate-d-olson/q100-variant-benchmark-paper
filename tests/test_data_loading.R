library(testthat)
library(here)

# Source the data loading functions
source(here::here("R/data_loading.R"))

test_that("parse_benchmark_id extracts all components", {
  # Test with full path
  result <- parse_benchmark_id(
    "results/var_counts/v5.0q_GRCh38_smvar/stratification_combined_metrics.csv"
  )

  expect_equal(result$bench_version, "v5.0q")
  expect_equal(result$ref, "GRCh38")
  expect_equal(result$var_type, "smvar")
})

test_that("parse_benchmark_id works with benchmark ID alone", {
  result <- parse_benchmark_id("v5.0q_GRCh38_smvar")

  expect_equal(result$bench_version, "v5.0q")
  expect_equal(result$ref, "GRCh38")
  expect_equal(result$var_type, "smvar")
})

test_that("parse_benchmark_id handles different versions", {
  # Test v4.2.1
  result <- parse_benchmark_id("v4.2.1_GRCh38_smvar")
  expect_equal(result$bench_version, "v4.2.1")

  # Test v0.6
  result <- parse_benchmark_id("v0.6_GRCh37_stvar")
  expect_equal(result$bench_version, "v0.6")
})

test_that("parse_benchmark_id handles different references", {
  result <- parse_benchmark_id("v5.0q_CHM13v2.0_smvar")
  expect_equal(result$ref, "CHM13v2.0")
})

test_that("parse_benchmark_id raises error for invalid format", {
  expect_error(
    parse_benchmark_id("invalid_benchmark_id"),
    "Could not parse benchmark ID"
  )
})

test_that("load_stratification_metrics returns correct structure", {
  metrics <- load_stratification_metrics()

  # Check that it's a tibble
  expect_s3_class(metrics, "tbl_df")

  # Check required columns exist
  required_cols <- c(
    "bench_version", "ref", "var_type", "strat_name",
    "strat_bp", "intersect_bp", "pct_of_strat", "pct_of_bench",
    "total_variants", "snp_count", "indel_count",
    "variant_density_per_mb"
  )
  expect_true(all(required_cols %in% names(metrics)))

  # Check that we have multiple benchmarks
  expect_gt(nrow(metrics), 40)
})

test_that("load_stratification_metrics has expected stratifications", {
  metrics <- load_stratification_metrics()

  expected_strats <- c("HP", "MAP", "SD", "SD10kb", "TR", "TR10kb")
  actual_strats <- unique(metrics$strat_name)

  expect_true(all(expected_strats %in% actual_strats))
})

test_that("load_stratification_metrics respects benchmark filter", {
  metrics_filtered <- load_stratification_metrics(
    benchmark_filter = c("v5.0q_GRCh38_smvar")
  )

  # Should have exactly 6 rows (one per stratification)
  expect_equal(nrow(metrics_filtered), 6)
  expect_equal(unique(metrics_filtered$bench_version), "v5.0q")
  expect_equal(unique(metrics_filtered$ref), "GRCh38")
  expect_equal(unique(metrics_filtered$var_type), "smvar")
})

test_that("load_stratification_metrics variant counts are non-negative", {
  metrics <- load_stratification_metrics()

  expect_true(all(metrics$total_variants >= 0))
  expect_true(all(metrics$snp_count >= 0))
  expect_true(all(metrics$indel_count >= 0))
})

test_that("load_exclusion_metrics returns correct structure or empty", {
  exclusions <- load_exclusion_metrics()

  # Either returns empty tibble or valid data
  if (nrow(exclusions) > 0) {
    expect_s3_class(exclusions, "tbl_df")

    required_cols <- c(
      "bench_version", "ref", "var_type", "exclusions",
      "exclusion_bp", "intersect_bp", "pct_of_exclusion", "pct_of_dip"
    )
    expect_true(all(required_cols %in% names(exclusions)))
  } else {
    expect_s3_class(exclusions, "tbl_df")
    expect_equal(nrow(exclusions), 0)
  }
})

test_that("load_exclusion_metrics only returns v5.0q data", {
  exclusions <- load_exclusion_metrics()

  if (nrow(exclusions) > 0) {
    # All exclusion data should be v5.0q
    expect_true(all(stringr::str_detect(
      exclusions$bench_version,
      "^v5\\.0q"
    )))
  }
})

test_that("load_reference_sizes returns correct structure", {
  ref_sizes <- load_reference_sizes()

  expect_s3_class(ref_sizes, "tbl_df")

  required_cols <- c("ref", "chrom", "length", "ns", "asm_bp")
  expect_true(all(required_cols %in% names(ref_sizes)))
})

test_that("load_reference_sizes calculates asm_bp correctly", {
  ref_sizes <- load_reference_sizes()

  # asm_bp should equal length - ns
  ref_sizes <- ref_sizes %>%
    dplyr::mutate(
      expected_asm_bp = length - ns,
      asm_correct = asm_bp == expected_asm_bp
    )

  expect_true(all(ref_sizes$asm_correct))
})

test_that("load_reference_sizes chromosomes have chr prefix", {
  ref_sizes <- load_reference_sizes()

  # All chromosome names should start with "chr"
  expect_true(all(stringr::str_detect(ref_sizes$chrom, "^chr")))
})

test_that("load_reference_sizes has expected references", {
  ref_sizes <- load_reference_sizes()

  expected_refs <- c("GRCh37", "GRCh38", "CHM13v2.0")
  actual_refs <- unique(ref_sizes$ref)

  expect_true(all(expected_refs %in% actual_refs))
})

test_that("load_variant_table handles invalid benchmark ID", {
  expect_error(
    load_variant_table("invalid_benchmark_id"),
    "Could not parse benchmark ID"
  )
})

test_that("load_variant_table errors on missing file", {
  expect_error(
    load_variant_table("v5.0q_nonexistent_genome_smvar"),
    "Variant file not found"
  )
})

test_that("load_diff_coverage returns correct structure", {
  # Use a known benchmark
  coverage <- load_diff_coverage("v5.0q_GRCh38_smvar")

  expect_s3_class(coverage, "tbl_df")

  required_cols <- c(
    "strat_name", "chrom", "start", "end",
    "n_overlap", "bases_cov", "ivl_len", "frac_cov"
  )
  expect_true(all(required_cols %in% names(coverage)))
})

test_that("load_diff_coverage has expected stratifications", {
  coverage <- load_diff_coverage("v5.0q_GRCh38_smvar")

  expected_strats <- c("HP", "MAP", "SD", "SD10kb", "TR", "TR10kb")
  actual_strats <- unique(coverage$strat_name)

  expect_true(all(expected_strats %in% actual_strats))
})

test_that("load_diff_coverage respects stratification filter", {
  coverage_filtered <- load_diff_coverage(
    "v5.0q_GRCh38_smvar",
    strat_filter = c("HP", "TR")
  )

  expect_true(all(unique(coverage_filtered$strat_name) %in% c("HP", "TR")))
})

test_that("load_diff_coverage fraction covered is between 0 and 1", {
  coverage <- load_diff_coverage("v5.0q_GRCh38_smvar")

  expect_true(all(coverage$frac_cov >= 0))
  expect_true(all(coverage$frac_cov <= 1))
})

test_that("load_diff_coverage errors on invalid benchmark", {
  expect_error(
    load_diff_coverage("v5.0q_nonexistent_genome_smvar"),
    "Coverage directory not found"
  )
})

test_that("stratification metrics and coverage have matching benchmarks", {
  metrics <- load_stratification_metrics()
  unique_benchmarks_metrics <- metrics %>%
    dplyr::distinct(bench_version, ref, var_type) %>%
    nrow()

  # Test with first benchmark
  benchmark_id <- "v5.0q_GRCh38_smvar"
  coverage <- load_diff_coverage(benchmark_id)

  # They should have the same stratifications
  metrics_strats <- metrics %>%
    dplyr::filter(
      bench_version == "v5.0q",
      ref == "GRCh38",
      var_type == "smvar"
    ) %>%
    dplyr::pull(strat_name) %>%
    unique() %>%
    sort()

  coverage_strats <- unique(coverage$strat_name) %>% sort()

  expect_equal(metrics_strats, coverage_strats)
})
