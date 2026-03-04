#!/bin/bash

# Pilot Migration Script
# Convert benchmark_unique_regions.qmd as pilot test

set -euo pipefail

echo "========================================="
echo "Pilot Migration: benchmark_unique_regions"
echo "========================================="
echo ""

# Check if setup was completed
if [ ! -f "tests/table_compatibility_test.html" ]; then
    echo "❌ Error: Test notebooks not rendered yet"
    echo "   Run setup_flextable_migration.sh first"
    exit 1
fi

echo "✓ Test notebooks found"
echo ""

# Source flextable theme in test
echo "Step 1: Testing flextable theme function..."
echo "--------------------------------------------"

Rscript -e '
source(here::here("R/plot_themes.R"))
library(tidyverse)
library(flextable)

# Quick test
iris %>%
  head(5) %>%
  flextable() %>%
  theme_flextable_manuscript()

cat("✓ Theme function works!\n")
'

echo ""
echo "Step 2: Creating backup of original notebook..."
echo "-----------------------------------------------"

cp analysis/benchmark_unique_regions.qmd analysis/benchmark_unique_regions.qmd.backup
echo "✓ Backed up to benchmark_unique_regions.qmd.backup"

echo ""
echo "Step 3: Applying flextable conversion..."
echo "----------------------------------------"

# The actual conversion will be done via replacement
echo "⚠️  Manual conversion required - see analysis/benchmark_unique_regions_flextable_pilot.qmd"
echo "   A pilot version has been created for you to review"

echo ""
echo "Step 4: Rendering pilot notebook..."
echo "-----------------------------------"

if [ -f "analysis/benchmark_unique_regions_flextable_pilot.qmd" ]; then
    quarto render analysis/benchmark_unique_regions_flextable_pilot.qmd --to html
    echo "✓ Rendered HTML"
    
    quarto render analysis/benchmark_unique_regions_flextable_pilot.qmd --to docx
    echo "✓ Rendered Word"
    
    # Open outputs
    open analysis/benchmark_unique_regions_flextable_pilot.html
    open analysis/benchmark_unique_regions_flextable_pilot.docx
else
    echo "⚠️  Pilot notebook not found, skipping render"
fi

echo ""
echo "========================================="
echo "Pilot Migration Complete!"
echo "========================================="
echo ""
echo "Next steps:"
echo "1. Review the pilot HTML and Word outputs"
echo "2. Test Word document in Google Docs"
echo "3. If satisfied, proceed with full migration:"
echo "   bash full_migration.sh"
echo ""
echo "To revert pilot changes:"
echo "   mv analysis/benchmark_unique_regions.qmd.backup analysis/benchmark_unique_regions.qmd"
echo ""
