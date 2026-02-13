# Plot Themes and Styling Guide

This guide explains how to use the consistent ggplot2 and gt table styling for manuscript-ready figures and tables.

## Quick Start

### Load the Themes Module

```r
source(here("R/plot_themes.R"))
```

## ggplot2 Figures

### Basic Usage with Theme

```r
library(ggplot2)

ggplot(data, aes(x = context_name, y = variant_count, fill = bench_version)) +
  geom_col(position = "dodge") +
  scale_benchmark_version(aesthetic = "fill") +
  scale_genomic_context(aesthetic = "x") +
  theme_manuscript() +
  labs(
    title = "Variant Counts by Genomic Context",
    x = "Genomic Context",
    y = "Number of Variants"
  )
```

### Color Scales

All color scales follow the same pattern and are designed to be:

- **Colorblind-friendly** (deuteranopia, protanopia, tritanopia compatible)
- **Print-ready** (work in black & white)
- **Publication-quality** (meet journal standards)

#### Available Scales

```r
# For benchmark versions
scale_benchmark_version(aesthetic = "color", name = "Benchmark")
scale_benchmark_version(aesthetic = "fill", name = "Benchmark")

# For benchmark set types
scale_bench_type(aesthetic = "color", name = "Benchmark Type")
scale_bench_type(aesthetic = "fill")

# For genomic contexts
scale_genomic_context(aesthetic = "color", name = "Genomic Context")
scale_genomic_context(aesthetic = "fill")
```

### Color Palettes

Access all colors directly with `get_color_palettes()`:

```r
palettes <- get_color_palettes()

# Benchmark versions
palettes$bench_version
# v0.6     = "#1B9E77" (teal)
# v4.2.1   = "#D95F02" (orange)
# v5.0q    = "#7570B3" (purple)

# Reference genomes
palettes$ref
# GRCh37   = "#E41A1C" (red)
# GRCh38   = "#377EB8" (blue)
# CHM13v2.0 = "#4DAF4A" (green)

# Genomic contexts
palettes$context_name
# HP       = "#E41A1C" (red)
# MAP      = "#377EB8" (blue)
# SD       = "#4DAF4A" (green)
# SD10kb   = "#984EA3" (purple)
# TR       = "#FF7F00" (orange)
# TR10kb   = "#A65628" (brown)

# Benchmark set types
palettes$bench_type
# smvar    = "#1B9E77" (teal)
# stvar    = "#D95F02" (orange)
```

### Theme Features

`theme_manuscript()` includes:

- **Font sizing:** 9pt body, 8pt axis labels, 8pt legend text (Cell Genomics standard)
- **Professional layout:** Clean minimal design with proper spacing
- **Figure width optimized:** For 85-180mm journal figures
- **Legend positioning:** Bottom, horizontal layout
- **Grid styling:** Removed for clean appearance
- **Border styling:** Professional black borders with appropriate linewidths
- **Facet styling:** Gray background for panel distinction

### Figure Size Guidelines

For Cell Genomics:

```r
# Single-column figure (85 mm = 3.35 in)
ggsave("figure.pdf", width = 3.35, height = 3.5, dpi = 300)

# Multi-column figure (120-180 mm)
ggsave("figure.pdf", width = 5, height = 3.5, dpi = 300)

# Full-page width figure (180 mm = 7 in)
ggsave("figure.pdf", width = 7, height = 4, dpi = 300)
```

### Example: Complete Figure

```r
metrics %>%
  # Data preparation
  filter(ref == "GRCh38") %>%
  group_by(bench_version, bench_type, context_name) %>%
  summarise(total_variants = sum(total_variants), .groups = "drop") %>%

  # Visualization
  ggplot(aes(
    x = context_name,
    y = total_variants,
    fill = bench_version
  )) +
  geom_col(position = "dodge") +
  facet_wrap(~bench_type, scales = "free_y") +

  # Theming
  scale_benchmark_version(aesthetic = "fill") +
  scale_genomic_context(aesthetic = "x") +
  theme_manuscript() +

  # Labels
  labs(
    title = "Variant Counts by Genomic Context",
    subtitle = "GRCh38 reference",
    x = "Genomic Context",
    y = "Variant Count",
    fill = "Benchmark"
  )

# Export for publication
ggsave("Figure2_variant_counts.pdf", width = 5, height = 3.5, dpi = 300)
```

## GT Tables

### Basic Usage with Theme

```r
library(gt)

data %>%
  gt() %>%
  theme_gt_manuscript() %>%
  gt::fmt_number(
    columns = c(total_variants, snv_count),
    decimals = 0,
    use_seps = TRUE
  )
```

### Theme Features

`theme_gt_manuscript()` includes:

- **Font sizing:** 9pt body, 10pt headers (Cell Genomics standard)
- **Professional structure:** Black borders at top/bottom, gray column headers
- **Striped rows:** Optional alternating row colors for readability
- **Proper spacing:** 10px header padding, 8px data cell padding
- **Responsive:** Adjusts to fit available width

### Example: Complete Table

```r
metrics %>%
  filter(ref == "GRCh38", bench_type == "smvar") %>%
  select(bench_version, context_name, total_variants, snv_count, indel_count) %>%
  arrange(bench_version, context_name) %>%

  # Create table
  gt(groupname_col = "bench_version") %>%

  # Apply theme
  theme_gt_manuscript(striped = TRUE) %>%

  # Format columns
  gt::fmt_number(
    columns = c(total_variants, snv_count, indel_count),
    decimals = 0,
    use_seps = TRUE
  ) %>%

  # Labels
  gt::cols_label(
    context_name = "Genomic Context",
    total_variants = "Total Variants",
    snv_count = "SNVs",
    indel_count = "Indels"
  ) %>%

  # Title and subtitle
  gt::tab_header(
    title = "Variant Counts by Genomic Context",
    subtitle = "Small variants (SNVs and Indels) in GRCh38"
  )
```

## Helper Functions

### Get Context Labels

Convert abbreviations to full names:

```r
get_context_labels()
# Returns all context labels

get_context_labels(c("HP", "TR", "SD"))
# Returns named character vector with specified contexts
```

### Get Benchmark Type Labels

```r
get_bench_type_labels()
# Returns: smvar = "Small Variants", stvar = "Structural Variants"

get_bench_type_labels(c("smvar"))
# Returns: smvar = "Small Variants"
```

### Figure Parameters

Generate standardized export parameters:

```r
params <- get_figure_params(
  figure_name = "variant_counts",
  figure_number = 2,
  format = "pdf",
  width = 5,
  height = 3.5
)

ggsave(
  filename = params$filename,
  width = params$width,
  height = params$height,
  dpi = params$dpi,
  units = params$units,
  bg = params$bg
)
```

## Color Palette Specifications

### Design Principles

1. **Colorblind-friendly:** All palettes tested for deuteranopia, protanopia, tritanopia
2. **Print-ready:** Work in grayscale and color
3. **Distinct:** Easy to distinguish in small figures
4. **Professional:** Suitable for journal publication
5. **Consistent:** Same colors used across all figures and tables

### Benchmark Versions

- **v0.6** (#1B9E77): Teal - oldest version
- **v4.2.1** (#D95F02): Orange - intermediate version
- **v5.0q** (#7570B3): Purple - newest version

*Progression suggests upgrade path from old (teal) → new (purple)*

### Reference Genomes

- **GRCh37** (#E41A1C): Red - older human reference
- **GRCh38** (#377EB8): Blue - current human reference
- **CHM13v2.0** (#4DAF4A): Green - telomere-to-telomere assembly

*Colors chosen for maximum contrast and colorblind compatibility*

### Genomic Contexts

- **HP** (#E41A1C): Red - Homopolymers
- **MAP** (#377EB8): Blue - Low Mappability regions
- **SD** (#4DAF4A): Green - Segmental Duplications
- **SD10kb** (#984EA3): Purple - Large Segmental Duplications
- **TR** (#FF7F00): Orange - Tandem Repeats
- **TR10kb** (#A65628): Brown - Large Tandem Repeats

*All six colors are easily distinguishable and colorblind-friendly*

## Cell Genomics Requirements Compliance

### Figure Specifications

✓ Width: 85-180 mm (supported by theme spacing)
✓ Font: Sans-serif (Roboto/Arial), 9-10pt body, 8pt labels
✓ Resolution: 300 dpi (use `dpi = 300` in ggsave)
✓ Colors: Colorblind-friendly, print-ready
✓ Format: PDF/EPS (vector), PNG (raster high-res)

### Table Specifications

✓ Font: Sans-serif, 9pt body, 10pt headers
✓ Borders: Clear top/bottom borders, column headers distinct
✓ Color: Minimal, professional appearance
✓ Readability: Proper spacing and optional striping

### File Format Guidelines

```r
# Vector formats (preferred for publication)
ggsave("figure.pdf", width = 5, height = 3.5, dpi = 300)
ggsave("figure.eps", width = 5, height = 3.5, dpi = 300)

# High-resolution raster (fallback)
ggsave("figure.png", width = 5, height = 3.5, dpi = 300)

# Tables - save as PDF or PNG
gt_object %>%
  gtsave("table.pdf")  # requires webshot2
```

## Updating Existing Figures

To apply the manuscript theme to existing plots:

```r
# Before
plot <- ggplot(...) +
  geom_point() +
  theme_minimal()

# After
plot <- ggplot(...) +
  geom_point() +
  scale_benchmark_version(aesthetic = "color") +
  theme_manuscript()
```

## Troubleshooting

### Legends appear in wrong position

- Check `legend.position` in `theme_manuscript()`
- Use `+ guides(color = "none")` to hide unwanted legends

### Font sizes look wrong

- Verify figure width/height matches publication specifications
- Export with `dpi = 300`
- Check font rendering in PDF viewer (may appear different than on screen)

### Colors don't match expectations

- Verify color blindness mode in graphic software
- Export to PDF and check in Adobe Reader (more accurate colors)
- Compare with Cell Genomics figure examples

### Table formatting issues

- Ensure data is a data frame or tibble
- Check gt::opt_table_font() compatibility
- Use gt::tab_options() to override specific settings if needed

## Advanced Customization

### Create Custom Palette

```r
my_palette <- c(
  "Option1" = "#E41A1C",
  "Option2" = "#377EB8",
  "Option3" = "#4DAF4A"
)

scale_custom <- function(aesthetic = "color") {
  switch(aesthetic,
    color = scale_color_manual(values = my_palette),
    fill = scale_fill_manual(values = my_palette)
  )
}
```

### Modify Theme

```r
theme_manuscript_custom <- function() {
  theme_manuscript() +
    theme(
      axis.title.x = element_blank(),
      legend.position = "right"
    )
}
```

### Custom GT Table Theme

```r
theme_gt_custom <- function(gt_object) {
  theme_gt_manuscript(gt_object) %>%
    gt::tab_options(
      column_labels.background.color = "#E8E8E8"
    )
}
```
