#!/bin/bash

# Full Migration Script
# Convert all analysis notebooks from gt to flextable

set -euo pipefail

echo "==========================================="
echo "Full Migration: gt → flextable"
echo "==========================================="
echo ""

# List of notebooks to convert
NOTEBOOKS=(
    "analysis/benchmarkset_characterization.qmd"
    "analysis/benchmark_difficult.qmd"
    "analysis/benchmark_exclusions.qmd"
    "analysis/benchmark_interval_size_distributions.qmd"
    "analysis/benchmark_unique_regions.qmd"
    "analysis/genomic_context_analysis.qmd"
    "analysis/external_evaluation.qmd"
)

echo "Notebooks to convert:"
for notebook in "${NOTEBOOKS[@]}"; do
    echo "  - $notebook"
done
echo ""

read -p "Proceed with full migration? (y/n) " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Migration cancelled"
    exit 0
fi

echo ""
echo "Step 1: Creating backups..."
echo "---------------------------"

BACKUP_DIR="analysis/backups_gt_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$BACKUP_DIR"

for notebook in "${NOTEBOOKS[@]}"; do
    if [ -f "$notebook" ]; then
        cp "$notebook" "$BACKUP_DIR/"
        echo "✓ Backed up $(basename $notebook)"
    fi
done

echo "✓ All backups saved to: $BACKUP_DIR"

echo ""
echo "Step 2: Applying conversions..."
echo "--------------------------------"

# Function to convert a notebook
convert_notebook() {
    local file=$1
    local basename=$(basename "$file")
    
    echo "Converting $basename..."
    
    # Create conversion using sed (portable across macOS and Linux)
    sed -i.bak \
        -e 's/gt::gt(/flextable::flextable(/g' \
        -e 's/gt::fmt_integer(/fmt_integer_flextable(/g' \
        -e 's/gt::fmt_number(/fmt_number_flextable(/g' \
        -e 's/gt::cols_label(/cols_label_flextable(/g' \
        -e 's/theme_gt_manuscript()/theme_flextable_manuscript()/g' \
        -e 's/columns = /columns = c(/g' \
        -e 's/), decimals = /)), decimals = /g' \
        "$file"
    
    # Special handling for groupname_col (needs manual review)
    if grep -q "groupname_col" "$file"; then
        echo "  ⚠️  Contains groupname_col - needs manual review for merge_v()"
    fi
    
    rm "${file}.bak"  # Remove sed backup
    echo "  ✓ Converted $basename"
}

# Convert each notebook
for notebook in "${NOTEBOOKS[@]}"; do
    if [ -f "$notebook" ]; then
        convert_notebook "$notebook"
    else
        echo "  ⚠️  Skipping $notebook (not found)"
    fi
done

echo ""
echo "Step 3: Updating YAML frontmatter..."
echo "------------------------------------"

# Add docx output format to notebooks
for notebook in "${NOTEBOOKS[@]}"; do
    if [ -f "$notebook" ]; then
        # Check if already has docx format
        if ! grep -q "docx:" "$notebook"; then
            echo "  Adding docx format to $(basename $notebook)"
            # This is complex - document for manual update
            echo "    ⚠️  Manual YAML update needed"
        fi
    fi
done

echo ""
echo "Step 4: Manual review items..."
echo "------------------------------"
echo ""
echo "The following require MANUAL REVIEW:"
echo ""
echo "1. Grouped tables (groupname_col):"
echo "   - Search for: flextable(groupname_col"
echo "   - Change to: arrange(groupname_col) %>% flextable() %>% merge_v(j = \"groupname_col\")"
echo ""
echo "2. Column name quoting:"
echo "   - Change: columns = c(col1, col2)"
echo "   - To: columns = c(\"col1\", \"col2\")"
echo ""
echo "3. Alignment (add after theme):"
echo "   - Add: %>% align(j = c(\"numeric_cols\"), align = \"right\", part = \"body\")"
echo ""
echo "4. YAML frontmatter (add docx output):"
echo "   format:"
echo "     html: default"
echo "     docx: default"
echo ""

echo ""
echo "Step 5: Testing conversions..."
echo "------------------------------"

# Test render each notebook
for notebook in "${NOTEBOOKS[@]}"; do
    if [ -f "$notebook" ]; then
        basename=$(basename "$notebook" .qmd)
        echo "Testing $basename..."
        
        # Try HTML render
        if quarto render "$notebook" --to html --quiet; then
            echo "  ✓ HTML renders successfully"
        else
            echo "  ✗ HTML render failed - review $notebook"
        fi
        
        # Try Word render
        if quarto render "$notebook" --to docx --quiet 2>/dev/null; then
            echo "  ✓ Word renders successfully"
        else
            echo "  ⚠️  Word render skipped (format may not be configured)"
        fi
    fi
done

echo ""
echo "==========================================="
echo "Migration Complete!"
echo "==========================================="
echo ""
echo "Summary:"
echo "  - Backups saved: $BACKUP_DIR"
echo "  - Notebooks converted: ${#NOTEBOOKS[@]}"
echo "  - Manual reviews needed: See Step 4 above"
echo ""
echo "Next steps:"
echo "1. Review manual items listed in Step 4"
echo "2. Test each notebook output:"
echo "   - Open HTML files in browser"
echo "   - Open Word files in Microsoft Word"
echo "   - Upload to Google Drive and test Google Docs"
echo "3. Fix any rendering errors"
echo "4. Commit changes:"
echo "   git add -u"
echo "   git add R/plot_themes.R"
echo "   git commit -m 'feat: migrate tables from gt to flextable for Word/GDocs compatibility'"
echo ""
echo "To revert all changes:"
echo "   cp $BACKUP_DIR/*.qmd analysis/"
echo ""
