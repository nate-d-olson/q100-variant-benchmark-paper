# Plot Themes and Styling Guide

Consistent ggplot2 and table styling for manuscript-ready figures and tables.
Designed for Cell Genomics submission requirements.

All exports live in `R/plot_themes.R` and are loaded by `analysis/_notebook_setup.R`.

## Quick Start

```r
source(here::here("R/plot_themes.R"))
```

## ggplot2 Figures

### Basic Usage

```r
ggplot(data, aes(x = context_name, y = variant_count, fill = bench_version)) +
  geom_col(position = "dodge") +
  scale_benchmark_version(aesthetic = "fill") +
  scale_genomic_context(aesthetic = "x") +
  theme_manuscript() +
  labs(x = "Genomic Context", y = "Number of Variants")
```

### Theme Function

`theme_manuscript(...)` — base manuscript theme.

- 8–11 pt fonts; 85–180 mm figure widths (Cell Genomics)
- Bottom horizontal legend, gray facet panels, clean minimal grid
- Accepts `...` to override any theme element:
  ```r
  theme_manuscript(axis.title = element_text(size = 12, face = "italic"))
  ```

### Color Scales

All scales accept `aesthetic`, `name`, `guide`, and `...` (passed to underlying ggplot2 scale).

| Function | Maps |
|---|---|
| `scale_benchmark_version()` | `v0.6`, `v4.2.1`, `v5.0q`, `PP` |
| `scale_bench_type()` | `smvar` → "Small Variants", `stvar` → "Structural Variants" |
| `scale_genomic_context()` | `HP`, `MAP`, `SD`, `SD10kb`, `TR`, `TR10kb` (auto-applies readable labels) |

Example with overrides:

```r
scale_benchmark_version(fill = "color", limits = c("v5.0q"))
```

### Color Palettes

Direct access via `get_color_palettes()`:

| Palette | Keys |
|---|---|
| `bench_version` | v0.6, v4.2.1, v5.0q, PP |
| `ref` | GRCh37, GRCh38, CHM13v2.0 |
| `bench_type` | smvar, stvar |
| `context_name` | HP, MAP, SD, SD10kb, TR, TR10kb |
| `chrom_type` | autosomes, sex chromosomes |
| `binary` | TRUE/FALSE, Yes/No |

All palettes are colorblind-friendly (deuteranopia/protanopia/tritanopia tested) and print-friendly.

### Figure Export

```r
# Single column (85 mm)
ggsave("figure.pdf", width = 3.35, height = 3.5, dpi = 300)
# 1.5 column (120 mm)
ggsave("figure.pdf", width = 4.7,  height = 3.5, dpi = 300)
# Full page (180 mm)
ggsave("figure.pdf", width = 7.0,  height = 4.0, dpi = 300)
```

`get_figure_params()` returns standardized export params:

```r
params <- get_figure_params(figure_name = "variant_counts", figure_number = 2,
                            format = "pdf", width = 5, height = 3.5)
ggsave(params$filename, width = params$width, height = params$height,
       dpi = params$dpi, units = params$units, bg = params$bg)
```

## Tables — flextable (primary)

`flextable` is the table library used across all analysis notebooks; it produces
Word/Google Docs–compatible output suitable for the manuscript.

### Basic Usage

```r
library(flextable)

data %>%
  flextable() %>%
  theme_flextable_manuscript() %>%
  fmt_integer_flextable(columns = c("total_variants", "snv_count"))
```

### Theme Function

`theme_flextable_manuscript(ft, striped = TRUE, font_family = "Roboto", base_font_size = 9, header_bg = "#F5F5F5", stripe_color = "#FAFAFA", ...)`

- Roboto 9 pt body / 10 pt bold header
- 2 px black top/bottom borders, 1 px header separator
- Optional alternating row stripes
- Autofit table layout

### Helper Functions

| Helper | Wraps | Purpose |
|---|---|---|
| `fmt_integer_flextable(ft, columns, big_mark = ",")` | `colformat_int()` | Thousands-separated integers |
| `fmt_number_flextable(ft, columns, decimals = 2)` | `colformat_double()` | Fixed-decimal numbers |
| `fmt_percent_flextable(ft, columns, decimals = 1, scale_values = FALSE)` | `set_formatter()` | Percent strings (set `scale_values = TRUE` for fractions like 0.85) |
| `cols_label_flextable(ft, ...)` | `set_header_labels()` | Rename header columns |
| `as_grouped_flextable(df, groupname_col)` | `as_grouped_data()` + `as_flextable()` | Grouped row headers (replaces `gt::gt(groupname_col = ...)`) |

### Example: Grouped Table with Formatting

```r
metrics %>%
  filter(ref == "GRCh38", bench_type == "smvar") %>%
  select(bench_version, context_name, total_variants, snv_count, indel_count) %>%
  arrange(bench_version, context_name) %>%
  as_grouped_flextable(groupname_col = "bench_version") %>%
  theme_flextable_manuscript(striped = TRUE) %>%
  fmt_integer_flextable(columns = c("total_variants", "snv_count", "indel_count")) %>%
  cols_label_flextable(
    context_name = "Genomic Context",
    total_variants = "Total Variants",
    snv_count = "SNVs",
    indel_count = "Indels"
  ) %>%
  flextable::set_caption("Variant counts by genomic context (GRCh38, small variants).")
```

## Tables — gt (legacy support)

`theme_gt_manuscript(gt_object, striped = TRUE, ...)` is still exported for
backward compatibility, but new tables should use flextable. Two notebooks
(`analysis/external_evaluation.qmd`) still contain `gt()` calls that have not
been migrated.

```r
data %>% gt() %>% theme_gt_manuscript()
```

## Label Helpers

```r
get_context_labels()                  # all context labels (named vector)
get_context_labels(c("HP", "TR"))     # subset
get_bench_type_labels()               # smvar = "Small Variants", stvar = "Structural Variants"
```

## Cell Genomics Compliance

- Figures: 85–180 mm width, 9 pt body / 8 pt labels, 300 dpi, vector (PDF/EPS) preferred
- Tables: sans-serif 9 pt body / 10 pt header, top/bottom borders, optional striping
- Colors: colorblind-friendly, render in grayscale

## Troubleshooting

**Wrong legend position** — pass `legend.position` via `theme_manuscript(legend.position = "right")` or `+ guides(color = "none")` to hide.

**Roboto font missing** — falls back to default sans; install Roboto or pass `font_family = "Arial"` to `theme_flextable_manuscript()`.

**Word export issues with `flextable::save_as_docx()`** — confirm `officer` is installed (already required by `flextable`).

**Stripes look wrong on small tables** — `theme_flextable_manuscript()` skips striping when `nrow <= 1`; for two-row tables disable with `striped = FALSE`.
