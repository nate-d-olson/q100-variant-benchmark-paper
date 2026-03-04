# R Tables for Word & Google Docs: Complete Research Guide

## Executive Summary

**Problem**: The `gt` package creates beautiful HTML tables but loses formatting when exported to Word or Google Docs - critical formats for manuscript collaboration.

**Solution**: Use `flextable` + `officer` packages for Word/Google Docs-compatible tables while maintaining manuscript-quality formatting.

**Impact**: Enables seamless sharing of formatted tables with collaborators using Word or Google Docs without manual reformatting.

---

## Package Comparison

### 1. flextable ⭐ RECOMMENDED

**GitHub**: https://github.com/davidgohel/flextable  
**Documentation**: https://ardata-fr.github.io/flextable-book/  
**CRAN**: https://cran.r-project.org/package=flextable

**Pros:**
- ✅ Excellent Word/Office integration via `officer` package
- ✅ Native .docx export with preserved formatting (borders, fonts, colors, alignment)
- ✅ Good HTML output for web viewing
- ✅ Flexible styling API comparable to `gt`
- ✅ Can create Word documents programmatically
- ✅ Active development (5+ years, regular updates)
- ✅ Large user community (~740k downloads/month)
- ✅ Works with Quarto and RMarkdown

**Cons:**
- ❌ Different API from `gt` (requires code refactoring)
- ❌ Slightly more verbose than `gt`
- ⚠️ HTML output less polished than `gt` (but still good)

**Best For:**
- Manuscript tables that need Word export
- Collaborative documents in Google Docs
- Reports that mix Word and HTML outputs
- Creating Word documents programmatically

**Code Example:**
```r
library(flextable)

data %>%
  flextable() %>%
  colformat_int(j = c("col1", "col2"), big.mark = ",") %>%
  colformat_double(j = "col3", digits = 2) %>%
  set_header_labels(col1 = "Column 1", col2 = "Column 2") %>%
  theme_booktabs() %>%  # or custom theme
  autofit()
```

---

### 2. officer

**GitHub**: https://github.com/davidgohel/officer  
**Documentation**: https://ardata-fr.github.io/officeverse/  
**CRAN**: https://cran.r-project.org/package=officer

**Purpose**: Companion to `flextable` for direct Word/PowerPoint manipulation.

**Capabilities:**
- ✅ Read and write Word (.docx) and PowerPoint (.pptx) files
- ✅ Add flextables to documents
- ✅ Modify existing documents (headers, footers, styles)
- ✅ Create documents from templates
- ✅ Programmatic report generation

**Use Cases:**
- Creating Word reports with R
- Adding tables to existing Word templates
- Batch report generation
- Custom formatted manuscripts

**Code Example:**
```r
library(officer)
library(flextable)

doc <- read_docx() %>%
  body_add_par("Title", style = "heading 1") %>%
  body_add_flextable(my_flextable) %>%
  body_add_par("Caption", style = "Normal")

print(doc, target = "report.docx")
```

---

### 3. kableExtra

**GitHub**: https://github.com/haozhu233/kableExtra  
**Documentation**: https://haozhu233.github.io/kableExtra/  
**CRAN**: https://cran.r-project.org/package=kableExtra

**Pros:**
- ✅ Simple API (extends base `kable()`)
- ✅ Excellent HTML output with Bootstrap styling
- ✅ Good LaTeX/PDF output
- ✅ Many styling options

**Cons:**
- ❌ **Poor Word export** - most styling lost in .docx
- ❌ Limited cross-format consistency
- ❌ Not designed for Word/Google Docs workflow

**Verdict**: Great for HTML and PDF, but NOT suitable for Word/Google Docs use case.

---

### 4. huxtable

**GitHub**: https://github.com/hughjonesd/huxtable  
**Documentation**: https://hughjonesd.github.io/huxtable/  
**CRAN**: https://cran.r-project.org/package=huxtable

**Pros:**
- ✅ Unified API for HTML, LaTeX, Word, Excel, and screen
- ✅ Consistent appearance across formats
- ✅ Good Word compatibility via `officer`
- ✅ Can export to multiple formats simultaneously

**Cons:**
- ❌ More verbose than `gt` or `flextable`
- ❌ Smaller community than alternatives
- ⚠️ Occasional formatting quirks
- ⚠️ More complex learning curve

**Verdict**: Good alternative to `flextable`, but `flextable` has better documentation and larger community.

---

### 5. gt (Current Package)

**GitHub**: https://github.com/rstudio/gt  
**Documentation**: https://gt.rstudio.com/  
**CRAN**: https://cran.r-project.org/package=gt

**Review:**
- ✅ Excellent HTML output
- ✅ Beautiful, polished design
- ✅ Intuitive API
- ✅ Great documentation
- ❌ **Poor Word export** - tables exported as images or lose formatting
- ❌ **Poor Google Docs compatibility** - copy/paste loses formatting
- ⚠️ Primarily designed for HTML output

**Recommendation**: Keep `gt` for HTML-only outputs (websites, HTML reports), but use `flextable` for tables that need Word/Google Docs.

---

## Workflow Recommendations

### Option 1: Full Migration to flextable (Recommended)

**Pros:**
- Single package for all table outputs
- Consistent code across all analyses
- Best Word/Google Docs compatibility

**Cons:**
- Requires refactoring all existing `gt` code
- HTML output slightly less polished than `gt`

**Steps:**
1. Install flextable + officer: `renv::install(c("flextable", "officer"))`
2. Create `R/flextable_themes.R` with manuscript theme function
3. Update analysis notebooks: replace `gt` calls with `flextable` equivalents
4. Update Quarto YAML: add `docx: default` format
5. Test rendering to HTML and Word

---

### Option 2: Hybrid Approach (gt for HTML, flextable for Word)

**Pros:**
- Keep `gt`'s beautiful HTML output
- Use `flextable` only when Word export needed
- No need to refactor existing code

**Cons:**
- Maintain two sets of table code
- More complex codebase
- Risk of inconsistent table appearance

**Implementation:**
```r
# Conditional table rendering based on output format
create_table <- function(data, output_format = knitr::opts_knit$get("rmarkdown.pandoc.to")) {
  if (output_format == "docx") {
    # Use flextable for Word
    data %>%
      flextable() %>%
      theme_flextable_manuscript()
  } else {
    # Use gt for HTML
    data %>%
      gt() %>%
      theme_gt_manuscript()
  }
}
```

---

### Option 3: Quarto-Native Approach (Not Recommended)

Quarto has some built-in table features, but they're limited:

```markdown
| Column1 | Column2 | Column3 |
|---------|---------|---------|
| Value1  | Value2  | Value3  |

: Table caption
```

**Limitations:**
- No conditional formatting
- Limited styling control
- No grouped headers
- Not suitable for complex manuscript tables

**Verdict**: Not a viable replacement for `gt` or `flextable`.

---

## Implementation Guide

### Step 1: Install Packages

```r
# Activate renv
renv::activate()

# Install flextable and officer
renv::install(c("flextable", "officer"))

# Update lockfile
renv::snapshot()
```

### Step 2: Create Theme Function

Create `R/flextable_themes.R` (see `tests/flextable_theme_implementation.qmd` for complete code):

```r
theme_flextable_manuscript <- function(
  ft,
  striped = TRUE,
  font_family = "Roboto",
  base_font_size = 9,
  header_bg = "#F5F5F5",
  stripe_color = "#FAFAFA"
) {
  # Theme implementation...
}
```

### Step 3: Update Analysis Files

**Before (gt):**
```r
data %>%
  gt() %>%
  fmt_integer(columns = c(col1, col2)) %>%
  fmt_number(columns = col3, decimals = 2) %>%
  cols_label(
    col1 = "Column 1",
    col2 = "Column 2",
    col3 = "Column 3"
  ) %>%
  theme_gt_manuscript()
```

**After (flextable):**
```r
data %>%
  flextable() %>%
  colformat_int(j = c("col1", "col2"), big.mark = ",") %>%
  colformat_double(j = "col3", digits = 2) %>%
  set_header_labels(
    col1 = "Column 1",
    col2 = "Column 2",
    col3 = "Column 3"
  ) %>%
  theme_flextable_manuscript() %>%
  align(j = c("col1", "col2", "col3"), align = "right", part = "body")
```

### Step 4: Update Quarto YAML

```yaml
---
title: "Analysis"
format:
  html: default
  docx: default  # Enable Word export
---
```

### Step 5: Render and Test

```r
# Render both formats
quarto::quarto_render("analysis/file.qmd", output_format = "all")

# Test Word output
# 1. Open .docx in Word - check formatting
# 2. Upload to Google Drive - open with Google Docs
# 3. Copy table from Word - paste into Google Doc
```

---

## Migration Checklist

- [ ] Install `flextable` and `officer` packages
- [ ] Create `R/flextable_themes.R` with theme functions
- [ ] Create helper functions to mimic `gt` API (fmt_integer_flextable, cols_label_flextable)
- [ ] Test rendering with sample data (see `tests/table_compatibility_test.qmd`)
- [ ] Update one analysis notebook as pilot test
- [ ] Verify Word export formatting matches requirements
- [ ] Test Google Docs upload/copy-paste workflow
- [ ] Update remaining analysis notebooks
- [ ] Update project documentation
- [ ] Update renv.lock with new dependencies
- [ ] Add Word rendering to Makefile if needed

---

## Additional Resources

### flextable Documentation
- **Getting Started**: https://ardata-fr.github.io/flextable-book/
- **Gallery**: https://ardata-fr.github.io/flextable-gallery/
- **Function Reference**: https://davidgohel.github.io/flextable/reference/index.html

### officer Documentation
- **Introduction**: https://ardata-fr.github.io/officeverse/
- **Word Documents**: https://davidgohel.github.io/officer/articles/offcran/word.html
- **Styling**: https://davidgohel.github.io/officer/articles/offcran/word_styles.html

### Comparison Articles
- **R Tables Comparison**: https://rfortherestofus.com/2019/11/how-to-make-beautiful-tables-in-r/
- **flextable vs gt**: https://stackoverflow.com/questions/60388313/
- **Word Export Guide**: https://rmarkdown.rstudio.com/articles_docx.html

### Example Projects Using flextable
- **gtsummary with flextable**: https://www.danieldsjoberg.com/gtsummary/articles/tbl_summary.html#customize-output
- **Clinical Trial Reports**: https://ardata-fr.github.io/flextable-book/examples.html

---

## Alternative Approaches (Not Recommended)

### 1. LaTeX + pandoc
- Convert tables to LaTeX, then use pandoc for Word
- **Issue**: Complex, brittle, poor editing experience

### 2. Excel Export + Manual Formatting
- Export to Excel, format manually, insert in Word
- **Issue**: Not reproducible, time-consuming

### 3. Screenshot/Image Export
- Export `gt` tables as images
- **Issue**: Not editable, poor quality, accessibility issues

### 4. HTML + Copy-Paste
- Render HTML, copy table, paste in Word
- **Issue**: Inconsistent results, formatting often lost

---

## Conclusion

**For manuscript tables that need Word/Google Docs compatibility, use flextable + officer.**

This provides:
- ✅ Professional formatting matching journal requirements
- ✅ Preserved styling in Word/Google Docs
- ✅ Reproducible workflow
- ✅ Programmatic table creation
- ✅ Active maintenance and support

**Implementation effort**: Moderate (1-2 days for full migration)  
**Long-term benefit**: High (solves persistent formatting issues)  
**Risk**: Low (well-established packages with large communities)
