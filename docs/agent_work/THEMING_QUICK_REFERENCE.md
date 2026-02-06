# Plot Theming Quick Reference

## Overview

Publication-ready ggplot2 and gt table styling system for Cell Genomics journal submissions with consistent colors, typography, and layout.

## Quick Start

### Load themes
```r
source(here("R/plot_themes.R"))
```

### Create a themed figure
```r
ggplot(data, aes(x = context_name, y = value, fill = bench_version)) +
  geom_col(position = "dodge") +
  scale_benchmark_version(aesthetic = "fill") +
  scale_genomic_context(aesthetic = "x") +
  theme_manuscript() +
  labs(title = "My Figure", x = "Context", y = "Value")
```

### Create a themed table
```r
data %>%
  gt() %>%
  theme_gt_manuscript() %>%
  gt::fmt_number(columns = c(col1, col2), decimals = 0, use_seps = TRUE)
```

## Color Palettes

### Benchmark Versions
```r
scale_benchmark_version(aesthetic = "fill")  # or "color"
# v0.6   → #1B9E77 (teal)
# v4.2.1 → #D95F02 (orange)
# v5.0q  → #7570B3 (purple)
```

### Reference Genomes
```r
scale_reference_genome(aesthetic = "fill")
# GRCh37    → #E41A1C (red)
# GRCh38    → #377EB8 (blue)
# CHM13v2.0 → #4DAF4A (green)
```

### Genomic Contexts
```r
scale_genomic_context(aesthetic = "fill")
# HP    → #E41A1C (red)
# MAP   → #377EB8 (blue)
# SD    → #4DAF4A (green)
# SD10kb → #984EA3 (purple)
# TR    → #FF7F00 (orange)
# TR10kb → #A65628 (brown)
```

### Variant Types
```r
scale_variant_type(aesthetic = "fill")
# smvar → #1B9E77 (teal)
# stvar → #D95F02 (orange)
```

## Theme Features

### theme_manuscript()
- ✓ 9pt body font, 8pt labels
- ✓ 85-180mm figure width optimized
- ✓ Professional minimal design
- ✓ Bottom legend positioning
- ✓ Black borders with proper linewidth
- ✓ Facet styling with gray background

### theme_gt_manuscript()
- ✓ 9pt body, 10pt headers
- ✓ Black top/bottom borders
- ✓ Gray column headers
- ✓ Optional striped rows
- ✓ Professional spacing (10px headers, 8px data)
- ✓ Responsive to figure width

## Export for Publication

```r
# Single column (85mm = 3.35in)
ggsave("figure.pdf", width = 3.35, height = 3.5, dpi = 300)

# Multi-column (120-180mm)
ggsave("figure.pdf", width = 5.5, height = 3.5, dpi = 300)

# Full width (180mm = 7in)
ggsave("figure.pdf", width = 7, height = 4, dpi = 300)

# Tables
gt_object %>% gtsave("table.pdf")
```

## Helper Functions

### Get context labels
```r
get_context_labels()                    # All labels
get_context_labels(c("HP", "TR"))       # Specific labels
```

### Get variant type labels
```r
get_variant_type_labels()               # All labels
get_variant_type_labels(c("smvar"))     # Specific labels
```

### Get all colors
```r
palettes <- get_color_palettes()
palettes$bench_version
palettes$context_name
```

### Figure export parameters
```r
params <- get_figure_params("my_figure", figure_number = 1)
# Returns: filename, width, height, dpi, units, bg, format
```

## Common Patterns

### Multi-series comparison
```r
ggplot(data, aes(x = ref, y = value, fill = bench_version)) +
  geom_col(position = "dodge") +
  scale_benchmark_version(aesthetic = "fill") +
  theme_manuscript()
```

### Faceted by genomic context
```r
ggplot(data, aes(x = bench_version, y = count, fill = var_type)) +
  geom_col() +
  facet_wrap(~context_name) +
  scale_variant_type(aesthetic = "fill") +
  theme_manuscript()
```

### Grouped table
```r
data %>%
  gt(groupname_col = "bench_version") %>%
  theme_gt_manuscript(striped = TRUE) %>%
  gt::fmt_number(columns = c(col1, col2), decimals = 0)
```

### Multi-panel figure
```r
library(patchwork)
p1 + p2 + p3 +
  plot_layout(guides = "collect") &
  theme_manuscript()
```

## Cell Genomics Requirements

| Requirement | Specification | Status |
|-------------|---------------|--------|
| Figure width | 85-180 mm | ✓ Supported |
| Font (body) | Sans-serif, 9pt | ✓ theme_manuscript() |
| Font (labels) | Sans-serif, 8pt | ✓ theme_manuscript() |
| Font (headers) | Sans-serif, 10pt | ✓ theme_gt_manuscript() |
| Resolution | 300 dpi | ✓ Use dpi = 300 |
| Colorblind | Yes | ✓ All palettes tested |
| Grayscale | Yes | ✓ Print-ready |
| Borders | Clear | ✓ Black borders |
| Color | Distinct | ✓ 6+ unique colors |

## Customization

### Override specific theme elements
```r
plot + theme_manuscript() +
  theme(
    legend.position = "right",
    axis.title.x = element_blank()
  )
```

### Create custom palette
```r
my_colors <- c("Option1" = "#E41A1C", "Option2" = "#377EB8")
scale_custom <- function(aes = "color") {
  switch(aes,
    color = scale_color_manual(values = my_colors),
    fill = scale_fill_manual(values = my_colors)
  )
}
```

## Examples

See full examples in `R/example_themed_analysis.R`:
- Variant counts by genomic context
- Variant density by context
- Summary statistics table
- Context coverage breakdown
- Multi-panel comprehensive figure
- Batch report generation

## Documentation

Complete guide: `docs/plot_themes_guide.md`
- Detailed usage instructions
- Color palette rationale
- Troubleshooting guide
- Advanced customization

## Design Principles

✓ **Colorblind-friendly:** Deuteranopia, protanopia, tritanopia compatible
✓ **Print-ready:** Works in black & white and color
✓ **Professional:** Suitable for peer-reviewed publication
✓ **Consistent:** Same colors across all figures and tables
✓ **Accessible:** Large fonts, high contrast
✓ **Simple:** Minimal design removes distractions

## File Formats

- **Vector (preferred):** PDF, EPS
- **Raster (fallback):** PNG at 300 dpi
- **Tables:** PDF, PNG

All formats render consistently with the manuscript theme.
