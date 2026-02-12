library(testthat)
library(here)

source(here::here("R/data_loading.R"))

test_that("load_exclusion_metrics loads variant counts", {
  # Mocking results directory or using existing if available
  # We know results/ exists from previous check
  df <- load_exclusion_metrics()

  if (nrow(df) > 0) {
    expect_true("total_variants" %in% names(df), info = "Should contain total_variants")
    expect_true("snv_count" %in% names(df), info = "Should contain snv_count (renamed from snp_count)")
    expect_true(is.numeric(df$total_variants))
  } else {
    skip("No exclusion files found to test")
  }
})
