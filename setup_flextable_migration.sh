#!/bin/bash

# Setup script for flextable migration
# Run this outside the sandbox environment

set -euo pipefail

echo "========================================="
echo "Flextable Migration Setup"
echo "========================================="
echo ""

# Check current branch
BRANCH=$(git branch --show-current)
echo "Current branch: $BRANCH"

if [ "$BRANCH" != "feat/migrate-tables-to-flextable" ]; then
    echo "⚠️  Warning: Not on feat/migrate-tables-to-flextable branch"
    read -p "Continue anyway? (y/n) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
fi

echo ""
echo "Step 1: Checking dependencies..."
echo "--------------------------------"

# Check for gettext (needed for gdtools compilation)
if command -v brew &> /dev/null; then
    if ! brew list gettext &> /dev/null; then
        echo "⚠️  gettext not found (required for gdtools compilation)"
        echo "   Installing via Homebrew..."
        brew install gettext
    else
        echo "✓ gettext found"
    fi
else
    echo "⚠️  Homebrew not found - may have issues compiling gdtools"
    echo "   If compilation fails, install Homebrew and run: brew install gettext"
fi

echo ""
echo "Step 2: Installing R packages..."
echo "--------------------------------"

Rscript -e '
cat("Installing flextable and officer...\n")
renv::install(c("flextable", "officer"))
cat("✓ Packages installed\n\n")

cat("Recording in renv.lock...\n")
renv::record(c("flextable", "officer"))
cat("✓ Recorded in renv.lock\n")
'

echo ""
echo "Step 3: Rendering test notebooks..."
echo "-----------------------------------"

quarto render tests/table_compatibility_test.qmd --to html
echo "✓ Rendered table_compatibility_test.html"

quarto render tests/table_compatibility_test.qmd --to docx
echo "✓ Rendered table_compatibility_test.docx"

quarto render tests/flextable_theme_implementation.qmd --to html
echo "✓ Rendered flextable_theme_implementation.html"

quarto render tests/flextable_theme_implementation.qmd --to docx
echo "✓ Rendered flextable_theme_implementation.docx"

echo ""
echo "Step 4: Opening outputs for review..."
echo "-------------------------------------"

# Open HTML in browser
open tests/table_compatibility_test.html
open tests/flextable_theme_implementation.html

# Open Word documents
open tests/table_compatibility_test.docx
open tests/flextable_theme_implementation.docx

echo ""
echo "========================================="
echo "Setup Complete!"
echo "========================================="
echo ""
echo "Next steps:"
echo "1. Review the HTML outputs in your browser"
echo "2. Review the Word documents"
echo "3. Upload Word documents to Google Drive and open in Google Docs"
echo "4. If formatting looks good, proceed with pilot migration"
echo ""
echo "To proceed with pilot migration:"
echo "  bash pilot_migration.sh"
echo ""
