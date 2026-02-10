library(testthat)
library(dplyr)

# Helper function to assign size bins (simulating what will be in the QMD or a helper)
assign_size_bin <- function(var_size) {
  abs_size <- abs(var_size)
  case_when(
    abs_size < 15 ~ "<15bp",
    abs_size >= 15 & abs_size < 50 ~ "15-49bp",
    TRUE ~ ">=50bp"
  )
}

test_that("assign_size_bin correctly categorizes variants", {
  sizes <- c(0, 1, 14, 15, 49, 50, 100, -5, -15, -49, -50)
  bins <- assign_size_bin(sizes)
  
  expected <- c(
    "<15bp", "<15bp", "<15bp",  # 0, 1, 14
    "15-49bp", "15-49bp",       # 15, 49
    ">=50bp", ">=50bp",         # 50, 100
    "<15bp",                    # -5
    "15-49bp", "15-49bp",       # -15, -49
    ">=50bp"                    # -50
  )
  
  expect_equal(bins, expected)
})
