# gt → flextable Migration Quick Reference

## Installation

```r
# Install packages
renv::install(c("flextable", "officer"))
renv::snapshot()

# Load in notebooks
library(flextable)
library(officer)
```

## Function Mapping

| gt Function | flextable Equivalent | Notes |
|-------------|----------------------|-------|
| `gt()` | `flextable()` | Create table |
| `fmt_integer()` | `colformat_int(j = , big.mark = ",")` | Use `j =` not `columns =` |
| `fmt_number(decimals = )` | `colformat_double(j = , digits = )` | Use `digits` not `decimals` |
| `fmt_percent()` | `colformat_double(j = , suffix = "%")` | Manual percent formatting |
| `cols_label()` | `set_header_labels()` | Named list syntax |
| `tab_header()` | `set_caption()` | Caption/title |
| `tab_stubhead()` | `set_header_labels(stub = "")` | Stub column label |
| `cols_align()` | `align(j = , align = , part = "body")` | Must specify `part` |
| `tab_options()` | Various `set_*()` functions | Different approach |
| `theme_gt_manuscript()` | `theme_flextable_manuscript()` | Custom function (see below) |

## Common Patterns

### Basic Table

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

### Grouped Table

**Before (gt):**

```r
data %>%
  gt(groupname_col = "group_var") %>%
  fmt_integer(columns = c(col1, col2)) %>%
  theme_gt_manuscript()
```

**After (flextable):**

```r
data %>%
  arrange(group_var) %>%
  flextable() %>%
  colformat_int(j = c("col1", "col2"), big.mark = ",") %>%
  merge_v(j = "group_var") %>%
  theme_flextable_manuscript() %>%
  fix_border_issues()
```

### Caption/Title

**Before (gt):**

```r
data %>%
  gt() %>%
  tab_header(
    title = "Table 1: Summary Statistics",
    subtitle = "By reference genome"
  )
```

**After (flextable):**

```r
data %>%
  flextable() %>%
  set_caption("Table 1: Summary Statistics")
# Note: No subtitle in flextable - add as separate paragraph
```

### Column Width

**Before (gt):**

```r
data %>%
  gt() %>%
  cols_width(
    col1 ~ px(100),
    col2 ~ px(150)
  )
```

**After (flextable):**

```r
data %>%
  flextable() %>%
  width(j = "col1", width = 1.0) %>%  # inches
  width(j = "col2", width = 1.5) %>%
  autofit()  # or use autofit() for all columns
```

## Helper Functions

Add these to `R/flextable_themes.R`:

```r
# Wrapper for integer formatting
fmt_integer_flextable <- function(ft, columns, big_mark = ",") {
  ft %>% colformat_int(j = columns, big.mark = big_mark)
}

# Wrapper for number formatting
fmt_number_flextable <- function(ft, columns, decimals = 2) {
  ft %>% colformat_double(j = columns, digits = decimals)
}

# Wrapper for column labels
cols_label_flextable <- function(ft, ...) {
  labels <- list(...)
  ft %>% set_header_labels(values = labels)
}

# Then use like:
data %>%
  flextable() %>%
  fmt_integer_flextable(columns = c("col1", "col2")) %>%
  fmt_number_flextable(columns = "col3", decimals = 4) %>%
  cols_label_flextable(
    col1 = "Column 1",
    col2 = "Column 2"
  )
```

## Theme Function

```r
theme_flextable_manuscript <- function(
  ft,
  striped = TRUE,
  font_family = "Roboto",
  base_font_size = 9,
  header_bg = "#F5F5F5",
  stripe_color = "#FAFAFA"
) {
  n_rows <- nrow(ft$body$dataset)
  
  ft_styled <- ft %>%
    # Fonts
    font(fontname = font_family, part = "all") %>%
    fontsize(size = base_font_size, part = "body") %>%
    fontsize(size = base_font_size + 1, part = "header") %>%
    bold(part = "header") %>%
    
    # Borders
    border_remove() %>%
    hline_top(border = fp_border(width = 2, color = "black"), part = "header") %>%
    hline_bottom(border = fp_border(width = 2, color = "black"), part = "body") %>%
    hline_bottom(border = fp_border(width = 1, color = "black"), part = "header") %>%
    
    # Colors
    bg(bg = header_bg, part = "header") %>%
    
    # Alignment
    align(align = "center", part = "header") %>%
    valign(valign = "middle", part = "all") %>%
    
    # Padding
    padding(padding = 8, part = "all") %>%
    
    # Width
    set_table_properties(layout = "autofit", width = 1.0)
  
  # Striped rows
  if (striped && n_rows > 1) {
    ft_styled <- ft_styled %>%
      bg(i = seq(2, n_rows, by = 2), bg = stripe_color, part = "body")
  }
  
  return(ft_styled)
}
```

## Quarto YAML

Update your notebook YAML to enable Word export:

```yaml
---
title: "Analysis"
format:
  html: default
  docx: default  # Add this
---
```

Or for more control:

```yaml
---
title: "Analysis"
format:
  html:
    toc: true
    code-fold: true
  docx:
    toc: false
    number-sections: true
---
```

## Common Gotchas

### 1. Column Names Must Be Quoted

**Wrong:**

```r
flextable() %>%
  colformat_int(j = c(col1, col2))  # Will fail
```

**Right:**

```r
flextable() %>%
  colformat_int(j = c("col1", "col2"))  # Works
```

### 2. Must Specify `part` for Body Styling

**Wrong:**

```r
flextable() %>%
  align(j = "col1", align = "right")  # Only aligns header
```

**Right:**

```r
flextable() %>%
  align(j = "col1", align = "right", part = "body")  # Aligns body cells
```

### 3. Grouped Tables Need `merge_v()`

**Wrong:**

```r
data %>%
  flextable(groupname_col = "group")  # groupname_col doesn't exist
```

**Right:**

```r
data %>%
  arrange(group) %>%
  flextable() %>%
  merge_v(j = "group")
```

### 4. Fix Borders After Merging

```r
flextable() %>%
  merge_v(j = "group") %>%
  fix_border_issues()  # Important!
```

## Alignment Reference

```r
# Center all headers (typical)
align(align = "center", part = "header")

# Right-align numeric columns
align(j = c("col1", "col2"), align = "right", part = "body")

# Left-align text columns
align(j = c("name", "category"), align = "left", part = "body")

# Vertical alignment
valign(valign = "middle", part = "all")
```

## Border Reference

```r
# Remove all borders
border_remove()

# Top border (2px black)
hline_top(border = fp_border(width = 2, color = "black"), part = "header")

# Bottom border (2px black)
hline_bottom(border = fp_border(width = 2, color = "black"), part = "body")

# Inner borders (1px gray)
hline(border = fp_border(width = 1, color = "#CCCCCC"))

# Outer box
border_outer(border = fp_border(width = 2, color = "black"))
```

## Testing Checklist

After converting a table:

- [ ] Render to HTML - check appearance
- [ ] Render to Word - verify formatting
- [ ] Check borders (top/bottom 2px black?)
- [ ] Check header styling (bold, gray background?)
- [ ] Check number formatting (commas, decimals?)
- [ ] Check alignment (headers centered, numbers right?)
- [ ] Check striping (alternating row colors?)
- [ ] Open in Word - manually edit to test editability
- [ ] Upload to Google Drive - open with Google Docs

## Example from Project

From `analysis/benchmark_unique_regions.qmd`:

**Original gt code:**

```r
old_only_summary_tbl %>%
  select(
    comparison_label,
    old_benchmark,
    new_benchmark,
    old_only_bp,
    new_only_bp,
    old_only_variants,
    new_only_variants,
    pct_old_only_of_old,
    pct_new_only_of_v5
  ) %>%
  gt::gt(groupname_col = "comparison_label") %>%
  gt::fmt_integer(columns = c(old_only_bp, new_only_bp, old_only_variants, new_only_variants)) %>%
  gt::fmt_number(columns = c(pct_old_only_of_old, pct_new_only_of_v5), decimals = 2) %>%
  gt::cols_label(
    old_benchmark = "Previous Benchmark",
    new_benchmark = "v5 Benchmark",
    old_only_bp = "Previous-Only Bases",
    new_only_bp = "v5-Only Bases",
    old_only_variants = "Previous-Only Variants",
    new_only_variants = "v5-Only Variants",
    pct_old_only_of_old = "% Previous Benchmark Bases",
    pct_new_only_of_v5 = "% v5 Benchmark Bases"
  ) %>%
  theme_gt_manuscript()
```

**Converted flextable code:**

```r
old_only_summary_tbl %>%
  select(
    comparison_label,
    old_benchmark,
    new_benchmark,
    old_only_bp,
    new_only_bp,
    old_only_variants,
    new_only_variants,
    pct_old_only_of_old,
    pct_new_only_of_v5
  ) %>%
  arrange(comparison_label) %>%
  flextable() %>%
  fmt_integer_flextable(
    columns = c("old_only_bp", "new_only_bp", "old_only_variants", "new_only_variants")
  ) %>%
  fmt_number_flextable(
    columns = c("pct_old_only_of_old", "pct_new_only_of_v5"),
    decimals = 2
  ) %>%
  cols_label_flextable(
    old_benchmark = "Previous Benchmark",
    new_benchmark = "v5 Benchmark",
    old_only_bp = "Previous-Only Bases",
    new_only_bp = "v5-Only Bases",
    old_only_variants = "Previous-Only Variants",
    new_only_variants = "v5-Only Variants",
    pct_old_only_of_old = "% Previous Benchmark Bases",
    pct_new_only_of_v5 = "% v5 Benchmark Bases"
  ) %>%
  theme_flextable_manuscript() %>%
  merge_v(j = "comparison_label") %>%
  fix_border_issues() %>%
  align(
    j = c("old_only_bp", "new_only_bp", "old_only_variants", 
          "new_only_variants", "pct_old_only_of_old", "pct_new_only_of_v5"),
    align = "right",
    part = "body"
  )
```

## Quick Test

Try this in your R console:

```r
# Load packages
library(tidyverse)
library(flextable)

# Test data
iris %>%
  head(10) %>%
  flextable() %>%
  theme_flextable_manuscript()
```

If it renders without errors, you're ready to migrate tables!
