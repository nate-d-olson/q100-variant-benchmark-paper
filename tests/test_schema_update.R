library(testthat)
library(here)

source(here::here("R/schemas.R"))

test_that("variant_table schema contains variant_length field", {
  schema <- get_arrow_schema("variant_table")
  expect_true("variant_length" %in% names(schema), 
              info = "Schema must contain variant_length field as per track requirements")
  expect_true(schema$variant_length$type == arrow::int32(),
              info = "variant_length must be int32")
})
