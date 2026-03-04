#!/usr/bin/env Rscript

# Test Script for Table Compatibility
# Run this to render all test notebooks and verify outputs

library(here)

cat("========================================\n")
cat("Table Compatibility Testing\n")
cat("========================================\n\n")

# Check if required packages are available
cat("Checking package availability...\n")

required_packages <- c("tidyverse", "flextable", "officer", "kableExtra", "huxtable", "gt", "quarto")
missing_packages <- c()

for (pkg in required_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  ✓ %s installed\n", pkg))
  } else {
    cat(sprintf("  ✗ %s NOT installed\n", pkg))
    missing_packages <- c(missing_packages, pkg)
  }
}

if (length(missing_packages) > 0) {
  cat("\n❌ Missing packages detected. Install with:\n")
  cat(sprintf("  renv::install(c(%s))\n", paste(shQuote(missing_packages), collapse = ", ")))
  cat("\nExiting...\n")
  quit(status = 1)
}

cat("\n✅ All required packages available\n\n")

# Render test notebooks
cat("Rendering test notebooks...\n\n")

notebooks <- c(
  "tests/table_compatibility_test.qmd",
  "tests/flextable_theme_implementation.qmd"
)

for (notebook in notebooks) {
  notebook_path <- here::here(notebook)

  if (!file.exists(notebook_path)) {
    cat(sprintf("  ⚠️  Skipping %s (not found)\n", notebook))
    next
  }

  cat(sprintf("  📄 Rendering %s...\n", notebook))

  tryCatch({
    quarto::quarto_render(notebook_path, quiet = TRUE)
    cat(sprintf("     ✓ Success\n"))
  }, error = function(e) {
    cat(sprintf("     ✗ Error: %s\n", e$message))
  })
}

cat("\n========================================\n")
cat("Testing Complete\n")
cat("========================================\n\n")

cat("Next steps:\n")
cat("1. Open HTML outputs in browser:\n")
cat("   - tests/table_compatibility_test.html\n")
cat("   - tests/flextable_theme_implementation.html\n\n")
cat("2. Open Word output:\n")
cat("   - tests/table_compatibility_test.docx\n\n")
cat("3. Test Google Docs compatibility:\n")
cat("   - Upload .docx to Google Drive\n")
cat("   - Open with Google Docs\n")
cat("   - Check formatting preservation\n\n")
cat("4. Read documentation:\n")
cat("   - docs/flextable-migration-summary.md\n\n")

cat("Done! 🎉\n")
