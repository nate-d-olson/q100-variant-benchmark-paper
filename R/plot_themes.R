#' Plot Themes and Styling for Manuscript
#'
#' Provides consistent ggplot2 and gt table styling for all figures and tables.
#' Designed for Cell Genomics journal submission requirements.
#'
#' @keywords internal
#'
#' @details
#' ## Color Palettes
#'
#' **Benchmark Versions:** Distinct colors for v0.6, v4.2.1, v5.0q
#' **Reference Genomes:** Colors for GRCh37, GRCh38, CHM13v2.0
#' **Variant Types:** Colors for smvar (small variants) and stvar (structural variants)
#' **Genomic Contexts:** Distinct colors for each difficult region type
#'
#' ## Theme Requirements (Cell Genomics)
#'
#' - Figure width: 85-180 mm (3.35-7 inches)
#' - Minimum font size: 7 pt (for legends), 8 pt (for axis labels)
#' - Resolution: 300 dpi for publication
#' - Color: Use colorblind-friendly palettes where possible
#' - File formats: PDF, EPS (scalable), or high-res PNG

library(tidyverse)
library(gt)

#' Color Palettes for Variables
#'
#' @return List of named color vectors for use in ggplot scales
#'
#' @keywords internal
#'
#' @details
#' All palettes are designed to be:
#' - Colorblind-friendly (deuteranopia, protanopia, tritanopia)
#' - Print-friendly (work in grayscale)
#' - Distinct under journal reproduction standards
#'
#' @export
get_color_palettes <- function() {
  list(
    # Benchmark versions - distinct, progressive colors
    bench_version = c(
      "v0.6" = "#1B9E77",      # Teal (oldest)
      "v4.2.1" = "#D95F02",    # Orange (intermediate)
      "v5.0q" = "#7570B3"      # Purple (newest)
    ),

    # Reference genomes - clearly distinguishable
    ref = c(
      "GRCh37" = "#E41A1C",    # Red
      "GRCh38" = "#377EB8",    # Blue
      "CHM13v2.0" = "#4DAF4A"  # Green
    ),

    # Variant types - simple contrast
    var_type = c(
      "smvar" = "#1B9E77",     # Teal (small variants)
      "stvar" = "#D95F02"      # Orange (structural variants)
    ),

    # Genomic contexts - distinctive palette for 6 categories
    context_name = c(
      "HP" = "#E41A1C",        # Red (Homopolymers)
      "MAP" = "#377EB8",       # Blue (Low Mappability)
      "SD" = "#4DAF4A",        # Green (Segmental Duplications)
      "SD10kb" = "#984EA3",    # Purple (Large SD)
      "TR" = "#FF7F00",        # Orange (Tandem Repeats)
      "TR10kb" = "#A65628"     # Brown (Large TR)
    ),

    # Chromosomes - grayscale for clarity
    chrom_type = c(
      "autosomes" = "#333333", # Dark gray
      "sex_chrom" = "#888888"  # Light gray
    ),

    # Boolean/categorical - binary contrast
    binary = c(
      "TRUE" = "#1B9E77",
      "FALSE" = "#D95F02",
      "Yes" = "#1B9E77",
      "No" = "#D95F02"
    )
  )
}

#' ggplot2 Theme for Manuscript Figures
#'
#' @return A ggplot2 theme object for use in figures
#'
#' @details
#' Theme specifications:
#' - Clean, minimal design suitable for journal publication
#' - Proper spacing for legends and labels
#' - Font sizes appropriate for 85-180mm figure width
#' - Removes gridlines for cleaner appearance
#' - Uses grayscale-friendly color specifications
#'
#' @examples
#' \dontrun{
#' ggplot(data) +
#'   geom_point() +
#'   theme_manuscript()
#' }
#'
#' @export
theme_manuscript <- function() {
  theme_minimal() +
    theme(
      # Font sizing for journal (assumes ~100mm figure width)
      text = element_text(
        family = "sans",
        size = 9,
        color = "black"
      ),
      plot.title = element_text(
        size = 11,
        face = "bold",
        hjust = 0,
        margin = margin(b = 6)
      ),
      plot.subtitle = element_text(
        size = 9,
        hjust = 0,
        margin = margin(b = 6)
      ),
      axis.title = element_text(
        size = 9,
        face = "bold",
        color = "black"
      ),
      axis.text = element_text(
        size = 8,
        color = "black"
      ),
      legend.title = element_text(
        size = 9,
        face = "bold"
      ),
      legend.text = element_text(
        size = 8
      ),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.margin = margin(t = 8),
      legend.key.width = unit(0.8, "cm"),
      legend.key.height = unit(0.3, "cm"),

      # Panel and axis styling
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.3
      ),
      axis.line = element_blank(),
      axis.ticks = element_line(
        color = "black",
        linewidth = 0.25
      ),
      axis.ticks.length = unit(0.2, "cm"),

      # Facet styling
      strip.text = element_text(
        size = 8,
        face = "bold",
        color = "black",
        margin = margin(5, 5, 5, 5)
      ),
      strip.background = element_rect(
        color = "gray90",
        fill = "gray95",
        linewidth = 0.3
      ),

      # Plot spacing
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "mm")
    )
}

#' Scale for Benchmark Versions
#'
#' @param aesthetic The ggplot aesthetic to apply (color, fill, etc.)
#' @param name Display name for the scale (default: "Benchmark")
#' @param guide Guide type ("legend" or "none")
#'
#' @return A ggplot2 scale object
#'
#' @keywords internal
#'
#' @export
scale_benchmark_version <- function(
    aesthetic = "color",
    name = "Benchmark",
    guide = "legend") {
  palettes <- get_color_palettes()

  switch(aesthetic,
    color = scale_color_manual(
      name = name,
      values = palettes$bench_version,
      guide = guide
    ),
    fill = scale_fill_manual(
      name = name,
      values = palettes$bench_version,
      guide = guide
    ),
    stop("Unknown aesthetic: ", aesthetic)
  )
}

#' Scale for Reference Genomes
#'
#' @param aesthetic The ggplot aesthetic to apply (color, fill, etc.)
#' @param name Display name for the scale (default: "Reference")
#' @param guide Guide type ("legend" or "none")
#'
#' @return A ggplot2 scale object
#'
#' @keywords internal
#'
#' @export
scale_reference_genome <- function(
    aesthetic = "color",
    name = "Reference",
    guide = "legend") {
  palettes <- get_color_palettes()

  switch(aesthetic,
    color = scale_color_manual(
      name = name,
      values = palettes$ref,
      guide = guide
    ),
    fill = scale_fill_manual(
      name = name,
      values = palettes$ref,
      guide = guide
    ),
    stop("Unknown aesthetic: ", aesthetic)
  )
}

#' Scale for Variant Types
#'
#' @param aesthetic The ggplot aesthetic to apply (color, fill, etc.)
#' @param name Display name for the scale (default: "Variant Type")
#' @param guide Guide type ("legend" or "none")
#'
#' @return A ggplot2 scale object
#'
#' @keywords internal
#'
#' @export
scale_variant_type <- function(
    aesthetic = "color",
    name = "Variant Type",
    guide = "legend") {
  palettes <- get_color_palettes()

  switch(aesthetic,
    color = scale_color_manual(
      name = name,
      values = palettes$var_type,
      labels = c(
        "smvar" = "Small Variants",
        "stvar" = "Structural Variants"
      ),
      guide = guide
    ),
    fill = scale_fill_manual(
      name = name,
      values = palettes$var_type,
      labels = c(
        "smvar" = "Small Variants",
        "stvar" = "Structural Variants"
      ),
      guide = guide
    ),
    stop("Unknown aesthetic: ", aesthetic)
  )
}

#' Scale for Genomic Contexts
#'
#' @param aesthetic The ggplot aesthetic to apply (color, fill, etc.)
#' @param name Display name for the scale (default: "Genomic Context")
#' @param guide Guide type ("legend" or "none")
#'
#' @return A ggplot2 scale object
#'
#' @keywords internal
#'
#' @export
scale_genomic_context <- function(
    aesthetic = "color",
    name = "Genomic Context",
    guide = "legend") {
  palettes <- get_color_palettes()

  # Create readable labels for genomic contexts
  context_labels <- c(
    "HP" = "Homopolymers",
    "MAP" = "Low Mappability",
    "SD" = "Segmental Duplications",
    "SD10kb" = "Large SDs (>10kb)",
    "TR" = "Tandem Repeats",
    "TR10kb" = "Large TRs (>10kb)"
  )

  switch(aesthetic,
    color = scale_color_manual(
      name = name,
      values = palettes$context_name,
      labels = context_labels,
      guide = guide
    ),
    fill = scale_fill_manual(
      name = name,
      values = palettes$context_name,
      labels = context_labels,
      guide = guide
    ),
    stop("Unknown aesthetic: ", aesthetic)
  )
}

#' GT Table Theme for Manuscript
#'
#' @param gt_object A gt table object
#' @param striped Logical: add striped rows for readability (default: TRUE)
#'
#' @return A gt table object with styling applied
#'
#' @details
#' Table styling specifications:
#' - Clean, professional appearance
#' - Appropriate font sizing for readability
#' - Proper spacing and borders for publication
#' - Alternating row colors for readability (optional)
#' - Color palette consistent with figures
#'
#' @examples
#' \dontrun{
#' data %>%
#'   gt() %>%
#'   theme_gt_manuscript()
#' }
#'
#' @export
theme_gt_manuscript <- function(gt_object, striped = TRUE) {
  # Base styling
  gt_object <- gt_object %>%
    gt::opt_table_font(
      font = list(
        gt::google_font("Roboto"),
        "Arial",
        "sans-serif"
      )
    ) %>%
    gt::tab_options(
      # Font sizing
      table.font.size = "9pt",
      table.heading.font.size = "10pt",
      stub.font.size = "9pt",
      summary_row.text_transform = "uppercase",

      # Table structure
      table.border.top.style = "solid",
      table.border.top.width = px(2),
      table.border.top.color = "black",
      table.border.bottom.style = "solid",
      table.border.bottom.width = px(2),
      table.border.bottom.color = "black",
      heading.border.bottom.style = "solid",
      heading.border.bottom.width = px(1),
      heading.border.bottom.color = "black",

      # Column labels
      column_labels.background.color = "#F5F5F5",
      column_labels.text_transform = "capitalize",
      column_labels.padding = px(10),
      column_labels.border.top.style = "solid",
      column_labels.border.top.width = px(1),
      column_labels.border.top.color = "black",
      column_labels.border.bottom.style = "solid",
      column_labels.border.bottom.width = px(1),
      column_labels.border.bottom.color = "black",

      # Data cells
      data_row.padding = px(8),
      table.width = pct(100),

      # Borders and spacing
      table.margin.left = "auto",
      table.margin.right = "auto"
    )

  # Add striped rows if requested
  if (striped) {
    gt_object <- gt_object %>%
      gt::opt_row_striping(row_striping = TRUE)
  }

  return(gt_object)
}

#' Format Numeric Columns for Publication
#'
#' @param gt_object A gt table object
#' @param columns Column names to format
#' @param decimals Number of decimal places to show
#' @param sep_mark Thousands separator ("," or ".")
#'
#' @return A gt table object with formatted columns
#'
#' @keywords internal
#'
#' @export
fmt_publication <- function(
    gt_object,
    columns,
    decimals = 2,
    sep_mark = ",") {
  gt_object %>%
    gt::fmt_number(
      columns = all_of(columns),
      decimals = decimals,
      sep_mark = sep_mark,
      use_seps = TRUE
    )
}

#' Helper: Add context labels to legend
#'
#' Creates standardized context labels for use in legends and tables
#'
#' @param context_names Vector of context abbreviations (HP, MAP, SD, etc.)
#'
#' @return Named character vector mapping abbreviations to full names
#'
#' @keywords internal
#'
#' @export
get_context_labels <- function(context_names = NULL) {
  all_labels <- c(
    "HP" = "Homopolymers",
    "MAP" = "Low Mappability",
    "SD" = "Segmental Duplications",
    "SD10kb" = "Large SDs (>10kb)",
    "TR" = "Tandem Repeats",
    "TR10kb" = "Large TRs (>10kb)"
  )

  if (is.null(context_names)) {
    return(all_labels)
  }

  all_labels[context_names]
}

#' Helper: Get benchmark set type labels
#'
#' Creates standardized benchmark set type labels
#'
#' @param bench_types Vector of benchmark set type codes ("smvar", "stvar")
#'
#' @return Named character vector mapping codes to full names
#'
#' @keywords internal
#'
#' @export
get_bench_type_labels <- function(bench_types = NULL) {
  all_labels <- c(
    "smvar" = "Small Variants",
    "stvar" = "Structural Variants"
  )

  if (is.null(bench_types)) {
    return(all_labels)
  }

  all_labels[bench_types]
}

#' Helper: Create publication-ready figure filename
#'
#' Generates standardized filenames for figure exports
#'
#' @param figure_name Name of the figure (e.g., "variant_counts", "coverage_map")
#' @param figure_number Figure number for manuscript
#' @param format File format ("pdf", "png", "eps")
#' @param width Figure width in inches
#' @param height Figure height in inches
#'
#' @return List with filename and export parameters
#'
#' @keywords internal
#'
#' @export
get_figure_params <- function(
    figure_name,
    figure_number,
    format = "pdf",
    width = 4.5,
    height = 3.5) {
  list(
    filename = sprintf("Figure%d_%s.%s", figure_number, figure_name, format),
    width = width,
    height = height,
    dpi = 300,
    units = "in",
    bg = "white",
    format = format
  )
}
