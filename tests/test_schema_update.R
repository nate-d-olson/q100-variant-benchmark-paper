library(testthat)
library(here)

source(here::here("R/schemas.R"))

test_that("variant_table schema does not contain redundant variant_length field", {
  schema <- get_arrow_schema("variant_table")
  expect_false("variant_length" %in% names(schema),
    info = "Schema must not contain redundant variant_length field"
  )
  expect_true("var_size" %in% names(schema),
    info = "Schema must contain var_size field"
  )
  expect_true(schema$var_size$type == arrow::int32(),
    info = "var_size must be int32"
  )
})
