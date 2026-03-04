# Table Migration Execution Guide

**Branch:** `feat/migrate-tables-to-flextable`  
**Created:** 2026-03-03  
**Purpose:** Step-by-step guide for migrating from gt to flextable packages

---

## Overview

This guide walks through the complete migration process from gt (HTML-optimized tables) to flextable (Word/Google Docs-compatible tables). The migration is organized into three phases:

1. **Setup & Testing** - Install packages and verify compatibility
2. **Pilot Migration** - Convert one notebook as test
3. **Full Migration** - Convert all remaining notebooks

---

## Prerequisites

- [x] Feature branch created: `feat/migrate-tables-to-flextable`
- [x] Theme functions added to `R/plot_themes.R`
- [x] Test notebooks created in `tests/`
- [x] Migration scripts created (setup, pilot, full)

---

## Phase 1: Setup & Testing (15-30 minutes)

### Step 1.1: Run Setup Script

```bash
# From project root
bash setup_flextable_migration.sh
```

**What it does:**
- Installs `flextable` and `officer` packages via renv
- Records packages in `renv.lock`
- Renders test notebooks to HTML and Word
- Opens outputs for review

**Expected outputs:**
- `tests/table_compatibility_test.html`
- `tests/table_compatibility_test.docx`
- `tests/flextable_theme_implementation.html`
- `tests/flextable_theme_implementation.docx`

### Step 1.2: Review Test Outputs

1. **HTML Review (in browser)**
   - Check table formatting (borders, fonts, alignment)
   - Verify striped rows
   - Confirm color scheme matches manuscript theme

2. **Word Review (in Microsoft Word)**
   - ✅ Tables should be native Word tables (editable)
   - ✅ Formatting should be preserved (borders, fonts, colors)
   - ✅ Number formatting should work (commas, decimals)
   - ✅ Headers should be bold with gray background

3. **Google Docs Review**
   - Upload `.docx` files to Google Drive
   - Open with Google Docs
   - ✅ Formatting should mostly be preserved
   - ✅ Tables should be editable
   - ⚠️ Minor color differences are acceptable

### Step 1.3: Decision Point

**If test outputs look good:**
→ Proceed to Phase 2 (Pilot Migration)

**If issues found:**
→ Adjust theme function in `R/plot_themes.R`
→ Re-run setup script
→ Review again

---

## Phase 2: Pilot Migration (30-60 minutes)

### Step 2.1: Review Pilot Script

```bash
# Review what the pilot will do
cat pilot_migration.sh
```

### Step 2.2: Run Pilot Migration

```bash
bash pilot_migration.sh
```

**What it does:**
- Tests theme function with real data
- Creates backup of `benchmark_unique_regions.qmd`
- Notifies about manual conversion steps

### Step 2.3: Manual Conversion (ONE TABLE)

Pick one simple table from `analysis/benchmark_unique_regions.qmd` and convert it manually:

**Before (gt):**

```r
unique_summary_tbl %>%
  select(comparison_label, old_benchmark, new_benchmark, old_only_bp) %>%
  gt::gt(groupname_col = "comparison_label") %>%
  gt::fmt_integer(columns = old_only_bp) %>%
  gt::cols_label(old_only_bp = "Previous-Only Bases") %>%
  theme_gt_manuscript()
```

**After (flextable):**

```r
unique_summary_tbl %>%
  select(comparison_label, old_benchmark, new_benchmark, old_only_bp) %>%
  arrange(comparison_label) %>%
  flextable() %>%
  fmt_integer_flextable(columns = c("old_only_bp")) %>%
  cols_label_flextable(old_only_bp = "Previous-Only Bases") %>%
  theme_flextable_manuscript() %>%
  merge_v(j = "comparison_label") %>%
  fix_border_issues() %>%
  align(j = c("old_only_bp"), align = "right", part = "body")
```

**Key changes:**
1. `gt::gt(groupname_col = "x")` → `arrange(x) %>% flextable() %>% merge_v(j = "x")`
2. `gt::fmt_integer` → `fmt_integer_flextable`
3. Column names must be quoted: `old_only_bp` → `"old_only_bp"`
4. Add alignment explicitly for numeric columns
5. Add `fix_border_issues()` after `merge_v()`

### Step 2.4: Test Pilot Render

```bash
# Render the pilot notebook
quarto render analysis/benchmark_unique_regions.qmd --to html
quarto render analysis/benchmark_unique_regions.qmd --to docx

# Open outputs
open analysis/benchmark_unique_regions.html
open analysis/benchmark_unique_regions.docx
```

### Step 2.5: Review Pilot Outputs

**HTML Checklist:**
- [ ] Table renders without errors
- [ ] Formatting looks professional
- [ ] Numbers are formatted correctly
- [ ] Headers are styled properly

**Word Checklist:**
- [ ] Table is native Word table
- [ ] Formatting is preserved
- [ ] Can edit table content
- [ ] Borders look correct

**Google Docs Checklist:**
- [ ] Upload .docx to Google Drive
- [ ] Open with Google Docs
- [ ] Table is editable
- [ ] Formatting mostly preserved

### Step 2.6: Decision Point

**If pilot successful:**
→ Proceed to Phase 3 (Full Migration)

**If issues found:**
→ Debug and fix
→ Update theme or conversion approach
→ Re-test pilot

**To revert pilot:**
```bash
mv analysis/benchmark_unique_regions.qmd.backup analysis/benchmark_unique_regions.qmd
```

---

## Phase 3: Full Migration (2-3 hours)

### Step 3.1: Review Migration Plan

```bash
# See which notebooks will be converted
cat full_migration.sh
```

**Notebooks to convert:**
- `analysis/benchmarkset_characterization.qmd`
- `analysis/benchmark_difficult.qmd`
- `analysis/benchmark_exclusions.qmd`
- `analysis/benchmark_interval_size_distributions.qmd`
- `analysis/benchmark_unique_regions.qmd`
- `analysis/genomic_context_analysis.qmd`
- `analysis/external_evaluation.qmd`

### Step 3.2: Run Full Migration

```bash
bash full_migration.sh
```

**What it does:**
- Creates timestamped backup directory
- Backs up all notebooks
- Applies automatic conversions (sed replacements)
- Identifies items needing manual review
- Tests renders

**Expected output:**
```
Backups saved to: analysis/backups_gt_20260303_143022/
Notebooks converted: 7
Manual reviews needed: See Step 4
```

### Step 3.3: Manual Review

The script will identify these items needing manual attention:

#### 1. Grouped Tables (groupname_col)

**Find:** `flextable(groupname_col`

**Fix:**
```r
# Before
gt::gt(groupname_col = "comparison_label")

# After
arrange(comparison_label) %>%
flextable() %>%
merge_v(j = "comparison_label") %>%
fix_border_issues()
```

#### 2. Column Name Quoting

**Find:** `columns = c([^"]`

**Fix:**
```r
# Before
columns = c(old_only_bp, new_only_bp)

# After
columns = c("old_only_bp", "new_only_bp")
```

#### 3. Add Alignment

**Required for:** All numeric columns

**Add after theme:**
```r
theme_flextable_manuscript() %>%
align(j = c("numeric_col1", "numeric_col2"), align = "right", part = "body")
```

#### 4. YAML Frontmatter

**Find:** Notebooks without `docx:` output

**Add:**
```yaml
---
title: "Analysis"
format:
  html:
    toc: true
    code-fold: true
  docx: default  # Add this line
---
```

### Step 3.4: Test All Renders

```bash
# Render each notebook
for notebook in analysis/*.qmd; do
    echo "Rendering $(basename $notebook)..."
    quarto render "$notebook" --to html
    quarto render "$notebook" --to docx
done
```

### Step 3.5: Review All Outputs

**Create review checklist:**

| Notebook | HTML OK | Word OK | GDocs OK | Notes |
|----------|---------|---------|----------|-------|
| benchmarkset_characterization | [ ] | [ ] | [ ] | |
| benchmark_difficult | [ ] | [ ] | [ ] | |
| benchmark_exclusions | [ ] | [ ] | [ ] | |
| benchmark_interval_size_distributions | [ ] | [ ] | [ ] | |
| benchmark_unique_regions | [ ] | [ ] | [ ] | |
| genomic_context_analysis | [ ] | [ ] | [ ] | |
| external_evaluation | [ ] | [ ] | [ ] | |

**For each notebook:**
1. Open HTML in browser - verify appearance
2. Open Word doc - verify formatting and editability
3. Upload to Google Drive - test Google Docs compatibility
4. Document any issues in Notes column

### Step 3.6: Fix Issues

**Common issues and fixes:**

| Issue | Solution |
|-------|----------|
| Column not found | Check column name quoting |
| Alignment wrong | Add explicit `align()` call |
| Borders missing after merge | Add `fix_border_issues()` |
| Error: object not found | Check function name (use flextable:: prefix if needed) |
| Numbers not formatted | Check `fmt_integer_flextable()` syntax |

---

## Phase 4: Finalization (30 minutes)

### Step 4.1: Update Documentation

1. **Update `R/plot_themes.R` documentation:**
   - Add examples of flextable functions
   - Document when to use gt vs flextable

2. **Update project README:**
   ```markdown
   ## Table Packages
   
   - **gt**: For HTML-only outputs (interactive reports, website)
   - **flextable**: For Word/Google Docs exports (manuscripts, reports)
   ```

3. **Update `docs/` if applicable**

### Step 4.2: Run Tests

```bash
# Run R tests
Rscript -e 'testthat::test_file("tests/test_data_loading.R")'
Rscript -e 'testthat::test_file("tests/test_cache.R")'

# Lint R code
Rscript -e 'lintr::lint_dir("R/")'
Rscript -e 'lintr::lint_dir("analysis/")'

# Format R code
Rscript -e 'styler::style_dir("R/")'
```

### Step 4.3: Update renv.lock

```r
# From R console
renv::snapshot()
```

### Step 4.4: Commit Changes

```bash
# Stage changes
git add R/plot_themes.R
git add analysis/*.qmd
git add tests/
git add *.sh
git add renv.lock

# Check status
git status

# Commit
git commit -m "feat: migrate tables from gt to flextable for Word/GDocs compatibility

- Add theme_flextable_manuscript() and helper functions to R/plot_themes.R
- Convert all analysis notebooks from gt to flextable
- Update YAML frontmatter to include docx output format
- Add test notebooks and migration scripts
- Update renv.lock with flextable and officer packages

Resolves formatting issues when exporting tables to Word or Google Docs.
Tables now render as native, editable Word tables with preserved formatting."

# Push
git push -u origin feat/migrate-tables-to-flextable
```

### Step 4.5: Create Pull Request

```bash
# Using GitHub CLI
gh pr create \
  --title "feat: Migrate tables from gt to flextable for Word/Google Docs compatibility" \
  --body "## Summary

Migrates table generation from gt package to flextable package to improve Word and Google Docs compatibility.

## Changes

- **R/plot_themes.R**: Added `theme_flextable_manuscript()` and helper functions
- **analysis/*.qmd**: Converted all gt tables to flextable
- **tests/**: Added test notebooks and comparison documentation
- **Scripts**: Added migration scripts for reproducibility

## Testing

- ✅ All notebooks render to HTML successfully
- ✅ All notebooks render to Word (.docx) successfully
- ✅ Word documents tested in Microsoft Word (native tables, editable)
- ✅ Word documents uploaded to Google Drive and tested in Google Docs
- ✅ Formatting preserved across all output formats

## Benefits

- Tables now compatible with Word and Google Docs
- Better collaboration with reviewers who use Word
- Native editable Word tables instead of images
- Consistent formatting across output formats

## Breaking Changes

None - gt package still available for HTML-only use cases.

## Documentation

See:
- `tests/TABLE_RECOMMENDATIONS_SUMMARY.md` - Decision rationale
- `tests/gt_to_flextable_quick_reference.md` - Migration guide
- `tests/table_packages_research.md` - Package comparison
" \
  --base main

# Or create PR via GitHub web interface
```

---

## Troubleshooting

### Issue: Packages not installing

```bash
# Try from R console
R
> renv::restore()
> renv::install(c("flextable", "officer"))
```

### Issue: Quarto rendering fails

```bash
# Check Quarto version
quarto --version  # Should be >= 1.3

# Update Quarto if needed
# Download from: https://quarto.org/docs/get-started/
```

### Issue: R temporary directory error

This occurs in sandboxed environments. Run scripts outside VS Code terminal:

```bash
# From native Terminal.app or iTerm
cd /Users/nolson/projects/q100-variant-benchmark-paper
bash setup_flextable_migration.sh
```

### Issue: Tables look different in Google Docs

**Expected:** Minor color differences are normal (Google Docs doesn't support all colors)

**Fix if major issues:**
- Simplify colors (use basic colors: #000000, #CCCCCC)
- Reduce border complexity
- Test with simpler theme

### Issue: Grouped tables have double borders

**Fix:** Add `fix_border_issues()` after `merge_v()`

```r
flextable() %>%
merge_v(j = "group_col") %>%
fix_border_issues() %>%  # Add this
theme_flextable_manuscript()
```

---

## Rollback Plan

### Rollback pilot only:

```bash
mv analysis/benchmark_unique_regions.qmd.backup analysis/benchmark_unique_regions.qmd
```

### Rollback full migration:

```bash
# Find backup directory
ls -la analysis/backups_gt_*

# Restore all files
cp analysis/backups_gt_YYYYMMDD_HHMMSS/*.qmd analysis/

# Verify
git diff analysis/
```

### Rollback branch:

```bash
git checkout main
git branch -D feat/migrate-tables-to-flextable
```

---

## Success Criteria

✅ **Phase 1 Complete:**
- [ ] Test notebooks render successfully
- [ ] Word outputs look professional
- [ ] Google Docs test shows acceptable formatting

✅ **Phase 2 Complete:**
- [ ] Pilot notebook renders without errors
- [ ] Pilot Word output verified
- [ ] Pilot Google Docs test successful

✅ **Phase 3 Complete:**
- [ ] All 7 notebooks converted
- [ ] All renders succeed (HTML + Word)
- [ ] Manual review items addressed
- [ ] All outputs tested

✅ **Phase 4 Complete:**
- [ ] Documentation updated
- [ ] Tests passing
- [ ] Changes committed
- [ ] Pull request created

---

## Timeline Estimate

| Phase | Task | Time | Total |
|-------|------|------|-------|
| 1 | Setup & Testing | 15-30 min | 30 min |
| 1 | Review outputs | 15 min | |
| 2 | Pilot migration | 30 min | 60 min |
| 2 | Test & review | 30 min | |
| 3 | Full migration | 1 hour | 3 hours |
| 3 | Manual reviews | 1 hour | |
| 3 | Testing | 1 hour | |
| 4 | Finalization | 30 min | 30 min |
| **Total** | | | **4-5 hours** |

---

## Next Actions

**RIGHT NOW:**

```bash
# 1. Run setup
bash setup_flextable_migration.sh

# 2. Review test outputs
open tests/table_compatibility_test.html
open tests/table_compatibility_test.docx

# 3. Upload to Google Drive for testing
# (Manual step - use Google Drive web interface)
```

**AFTER SETUP SUCCESSFUL:**

```bash
# 4. Run pilot migration
bash pilot_migration.sh
```

**AFTER PILOT SUCCESSFUL:**

```bash
# 5. Run full migration
bash full_migration.sh
```

---

## Contact & Support

- **Documentation**: See `tests/` directory for detailed guides
- **Issues**: Review `docs/troubleshooting.md`
- **Questions**: Create GitHub issue or discussion

---

**Good luck with the migration! 🚀**
