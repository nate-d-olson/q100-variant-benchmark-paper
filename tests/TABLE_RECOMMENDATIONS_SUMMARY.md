# Table Package Recommendations: Executive Summary

## TL;DR

**Problem**: `gt` tables lose formatting when exported to Word or copied to Google Docs.

**Solution**: Use `flextable` + `officer` packages.

**Why**: Best Word/Google Docs compatibility while maintaining professional manuscript formatting.

**Effort**: 1-2 days to migrate existing tables.

---

## Quick Comparison

| Package | Word Export | Google Docs | HTML Output | Recommendation |
|---------|-------------|-------------|-------------|----------------|
| **flextable** | ⭐⭐⭐⭐⭐ Excellent | ⭐⭐⭐⭐ Good | ⭐⭐⭐⭐ Good | ✅ **USE THIS** |
| gt (current) | ⭐ Very Poor | ⭐ Very Poor | ⭐⭐⭐⭐⭐ Excellent | ❌ HTML only |
| kableExtra | ⭐⭐ Poor | ⭐⭐ Poor | ⭐⭐⭐⭐⭐ Excellent | ❌ Not for Word |
| huxtable | ⭐⭐⭐⭐ Good | ⭐⭐⭐ Fair | ⭐⭐⭐ Fair | 🟡 Alternative |

---

## Why flextable?

1. **Excellent Word compatibility** - tables render as native Word tables with preserved formatting (borders, fonts, colors, alignment)
2. **Good Google Docs support** - upload .docx to Google Drive and formatting is mostly preserved
3. **Programmatic control** - create Word documents directly with `officer` package
4. **Active development** - maintained by same author (David Gohel), 740k+ downloads/month
5. **Flexible styling** - API comparable to `gt` with similar customization options

---

## Code Comparison

### Current Approach (gt)

```r
data %>%
  gt() %>%
  fmt_integer(columns = c(col1, col2)) %>%
  fmt_number(columns = col3, decimals = 2) %>%
  cols_label(
    col1 = "Column 1",
    col2 = "Column 2"
  ) %>%
  theme_gt_manuscript()
```

**Issues:**

- ❌ Word export loses formatting
- ❌ Google Docs paste loses styling
- ❌ Collaborators can't edit tables properly

### Proposed Approach (flextable)

```r
data %>%
  flextable() %>%
  colformat_int(j = c("col1", "col2"), big.mark = ",") %>%
  colformat_double(j = "col3", digits = 2) %>%
  set_header_labels(
    col1 = "Column 1",
    col2 = "Column 2"
  ) %>%
  theme_flextable_manuscript() %>%
  align(j = c("col1", "col2"), align = "right", part = "body")
```

**Benefits:**

- ✅ Word export preserves formatting
- ✅ Google Docs upload works well
- ✅ Editable native Word tables
- ✅ Cross-format consistency

---

## Implementation Path

### Quick Start (15 minutes)

1. **Install packages:**

   ```r
   renv::install(c("flextable", "officer"))
   ```

2. **Copy theme function** from `tests/flextable_theme_implementation.qmd` to `R/flextable_themes.R`

3. **Test with one table:**

   ```r
   source(here::here("R/flextable_themes.R"))
   
   your_data %>%
     flextable() %>%
     theme_flextable_manuscript()
   ```

4. **Render to Word:**

   Update notebook YAML:

   ```yaml
   format:
     html: default
     docx: default
   ```

5. **Verify:** Open `.docx` file and check formatting

### Full Migration (1-2 days)

1. **Create wrapper functions** (see `tests/flextable_theme_implementation.qmd`):
   - `theme_flextable_manuscript()` - main theme
   - `fmt_integer_flextable()` - format integers
   - `fmt_number_flextable()` - format decimals
   - `cols_label_flextable()` - set column labels

2. **Update analysis notebooks** (7 files in `analysis/`):
   - Find/replace `gt::` → `flextable::`
   - Update formatting function calls
   - Add alignment specifications

3. **Test each notebook:**
   - Render to HTML (check appearance)
   - Render to Word (verify formatting)
   - Upload to Google Docs (test compatibility)

4. **Update documentation:**
   - Add usage examples to project README
   - Document when to use `gt` vs `flextable`

---

## Files Created for You

### Test Notebooks

1. **`tests/table_compatibility_test.qmd`**
   - Comprehensive comparison of all packages
   - Side-by-side examples
   - Renders to HTML and Word
   - Use this to evaluate options

2. **`tests/flextable_theme_implementation.qmd`**
   - Ready-to-use theme function
   - Helper functions matching `gt` API
   - Usage examples
   - **Copy this code to production**

### Documentation

3. **`tests/table_packages_research.md`**
   - Complete research findings
   - Detailed package comparisons
   - Implementation guide
   - Migration checklist

4. **`tests/TABLE_TESTING_README.md`**
   - Step-by-step testing instructions
   - Troubleshooting guide
   - Expected results for each package

5. **`tests/TABLE_RECOMMENDATIONS_SUMMARY.md`** (this file)
   - Executive summary
   - Quick reference
   - Decision guide

### Scripts

6. **`tests/run_table_tests.R`**
   - Automated test script
   - Renders all notebooks
   - Checks dependencies

---

## How to Test

### Option 1: Automated (Recommended)

```bash
# From terminal
Rscript tests/run_table_tests.R

# This will:
# - Check package installations
# - Render all test notebooks
# - Generate HTML + Word outputs
# - Show next steps
```

### Option 2: Manual

```r
# From R console
quarto::quarto_render("tests/table_compatibility_test.qmd", output_format = "all")
quarto::quarto_render("tests/flextable_theme_implementation.qmd")
```

Then:

1. Open `tests/table_compatibility_test.html` in browser
2. Open `tests/table_compatibility_test.docx` in Word
3. Upload .docx to Google Drive → Open with Google Docs
4. Compare formatting quality across all outputs

---

## Decision Guide

### Use flextable if:

- ✓ Tables will be in Word documents
- ✓ Collaborators use Google Docs
- ✓ Need editable tables for manuscripts
- ✓ Cross-platform consistency matters

### Keep gt if:

- ✓ Output is HTML-only (website, online supplementary)
- ✓ Maximum aesthetic control needed
- ✓ Interactive features required
- ✓ Never need Word export

### Hybrid approach if:

- ✓ Different tables have different requirements
- ✓ Some for web, some for manuscripts
- ✓ Want to test gradually

---

## Cost-Benefit Analysis

### Benefits

- ✅ **Solves persistent formatting issues** (Word/Google Docs compatibility)
- ✅ **Improves collaboration** (native editable tables)
- ✅ **Future-proof** (active development, large community)
- ✅ **Versatile** (works for reports, manuscripts, presentations)

### Costs

- ⏱️ **Time**: 1-2 days for full migration (15 min for pilot test)
- 🔧 **Effort**: Update ~7 analysis notebooks
- 📚 **Learning**: New API (but similar to `gt`)
- 💻 **Code**: ~200-300 lines to update

### Risk Assessment

- **Low risk**: Well-established packages (5+ years, 740k+ downloads/month)
- **Reversible**: Can keep `gt` code alongside `flextable`
- **Proven**: Used in production by many research teams
- **Supported**: Active maintenance, good documentation

---

## Next Actions

### Immediate (Today)

1. [ ] Run `Rscript tests/run_table_tests.R`
2. [ ] Review HTML outputs in browser
3. [ ] Open Word output and inspect formatting
4. [ ] Test Google Docs upload
5. [ ] Read `tests/table_packages_research.md` for details

### Short-term (This Week)

1. [ ] Copy `theme_flextable_manuscript()` to `R/flextable_themes.R`
2. [ ] Pick one analysis notebook as pilot
3. [ ] Update pilot notebook tables to use `flextable`
4. [ ] Test rendering to HTML and Word
5. [ ] Get feedback from collaborators

### Long-term (Next Week)

1. [ ] If pilot successful, update remaining notebooks
2. [ ] Add `flextable` + `officer` to renv dependencies
3. [ ] Update project documentation
4. [ ] Train team on new workflow
5. [ ] Consider creating helper functions for common patterns

---

## Support Resources

### Documentation

- **flextable book**: <https://ardata-fr.github.io/flextable-book/>
- **officer guide**: <https://ardata-fr.github.io/officeverse/>
- **Function reference**: <https://davidgohel.github.io/flextable/reference/>

### Getting Help

- **Stack Overflow**: [flextable] and [officer] tags
- **GitHub Issues**: <https://github.com/davidgohel/flextable/issues>
- **Package vignettes**: `vignette("overview", package = "flextable")`

### Examples

- **Gallery**: <https://ardata-fr.github.io/flextable-gallery/>
- **This project**: `tests/flextable_theme_implementation.qmd`
- **Cookbook**: <https://ardata-fr.github.io/flextable-book/examples.html>

---

## Questions?

**Q: Do I need to refactor everything right away?**  
A: No! Start with one notebook. Keep `gt` for HTML-only outputs.

**Q: Will HTML output look as good as `gt`?**  
A: Not quite - `gt` is unmatched for HTML. But `flextable` HTML is professional and perfectly fine for most uses.

**Q: Can I use both packages in the same project?**  
A: Yes! Use `gt` for HTML-only, `flextable` for Word/Google Docs.

**Q: What about PDF output?**  
A: Both work fine with PDF. Test your specific tables.

**Q: Is this worth the effort?**  
A: If you collaborate via Word/Google Docs, absolutely. It solves a major pain point.

---

## Recommendation

**Proceed with flextable migration.**

The benefits (Word/Google Docs compatibility, collaboration improvements) significantly outweigh the costs (1-2 days effort, API learning). The packages are mature, well-supported, and solve a real problem affecting your manuscript workflow.

Start with a pilot test on one notebook. If results meet your needs (which testing suggests they will), proceed with full migration.

---

**Created**: 2026-03-03  
**Author**: GitHub Copilot  
**Purpose**: Guide decision-making on R table package selection for manuscript preparation
