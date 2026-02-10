#' Example: Using Plot Themes in Analysis
#'
#' This file demonstrates how to apply the manuscript themes to create
#' publication-ready figures and tables using the data loading functions
#' and genomic context terminology.
#'
#' @keywords internal

# Setup
library(tidyverse)
library(here)
library(gt)

# Load functions
source(here("R/data_loading.R"))
source(here("R/plot_themes.R"))

# ============================================================================
# EXAMPLE 1: Variant Counts by Genomic Context (Figure)
# ============================================================================

example_figure_variant_counts <- function() {
  # Load data
  metrics <- load_genomic_context_metrics()

  # Create figure
  plot <- metrics %>%
    # Filter to specific reference and context
    filter(ref == "GRCh38") %>%
    group_by(bench_version, var_type, context_name) %>%
    summarise(
      snv_count = sum(snv_count),
      indel_count = sum(indel_count),
      total_variants = snv_count + indel_count,
      .groups = "drop"
    ) %>%
    # Visualization
    ggplot(aes(
      x = context_name,
      y = total_variants,
      fill = bench_version
    )) +
    geom_col(position = "dodge", color = "black", linewidth = 0.3) +
    facet_wrap(~var_type, scales = "free_y") +

    # Apply consistent theming
    scale_benchmark_version(aesthetic = "fill", name = "Benchmark") +
    theme_manuscript() +

    # Labels and annotations
    labs(
      title = "Variant Counts by Genomic Context",
      subtitle = "GRCh38 reference genome",
      x = "Genomic Context",
      y = "Total Variant Count",
      fill = "Benchmark Version"
    ) +

    # Fine-tune spacing
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    )

  return(plot)
}

# Export figure
# p <- example_figure_variant_counts()
# ggsave(
#   filename = "Figure1_variant_counts_by_context.pdf",
#   plot = p,
#   width = 5.5,
#   height = 3.5,
#   dpi = 300,
#   units = "in"
# )


# ============================================================================
# EXAMPLE 2: Variant Density by Context (Faceted Figure)
# ============================================================================

example_figure_variant_density <- function() {
  metrics <- load_genomic_context_metrics()

  plot <- metrics %>%
    filter(ref == "GRCh38", var_type == "smvar") %>%
    ggplot(aes(
      x = bench_version,
      y = variant_density_per_mb,
      color = context_name,
      group = context_name
    )) +
    geom_line(linewidth = 0.7, alpha = 0.8) +
    geom_point(size = 2.5) +
    facet_wrap(~context_name, ncol = 3) +

    # Apply themes
    scale_genomic_context(aesthetic = "color") +
    theme_manuscript() +
    labs(
      title = "Variant Density by Genomic Context",
      x = "Benchmark Version",
      y = "Variants per Megabase",
      color = "Genomic Context"
    ) +
    theme(
      legend.position = "none" # Remove legend since facets show context
    )

  return(plot)
}


# ============================================================================
# EXAMPLE 3: Summary Statistics Table
# ============================================================================

example_table_summary <- function() {
  metrics <- load_genomic_context_metrics()

  table <- metrics %>%
    filter(ref == "GRCh38") %>%
    group_by(bench_version, var_type) %>%
    summarise(
      `Total Variants` = sum(total_variants),
      `SNVs` = sum(snv_count),
      `Indels` = sum(indel_count),
      `Avg Density` = mean(variant_density_per_mb),
      .groups = "drop"
    ) %>%
    arrange(bench_version, var_type) %>%
    # Create gt table
    gt(groupname_col = "bench_version") %>%
    # Apply theme
    theme_gt_manuscript(striped = TRUE) %>%
    # Format numeric columns
    gt::fmt_number(
      columns = c(`Total Variants`, `SNVs`, `Indels`),
      decimals = 0,
      use_seps = TRUE
    ) %>%
    gt::fmt_number(
      columns = `Avg Density`,
      decimals = 1
    ) %>%
    # Labels
    gt::cols_label(
      var_type = "Variant Type",
      `Total Variants` = "Total Variants",
      `SNVs` = "Single Nucleotide Variants",
      `Indels` = "Insertions/Deletions",
      `Avg Density` = "Avg Density (per Mb)"
    ) %>%
    # Title
    gt::tab_header(
      title = "Variant Summary by Benchmark and Type",
      subtitle = "GRCh38 reference genome"
    )

  return(table)
}

# Export table
# table <- example_table_summary()
# gt::gtsave(table, filename = "Table1_variant_summary.pdf")


# ============================================================================
# EXAMPLE 4: Context Coverage Breakdown Table
# ============================================================================

example_table_context_coverage <- function() {
  metrics <- load_genomic_context_metrics()

  table <- metrics %>%
    filter(ref == "GRCh38", bench_version == "v5.0q") %>%
    select(var_type, context_name, context_bp, intersect_bp, pct_of_context) %>%
    rename(
      `Variant Type` = var_type,
      `Genomic Context` = context_name,
      `Context Size (bp)` = context_bp,
      `Overlapping (bp)` = intersect_bp,
      `Coverage %` = pct_of_context
    ) %>%
    arrange(`Variant Type`, `Genomic Context`) %>%
    # Create table
    gt(groupname_col = "`Variant Type`") %>%
    # Apply theme
    theme_gt_manuscript(striped = FALSE) %>%
    # Format columns
    gt::fmt_number(
      columns = c(`Context Size (bp)`, `Overlapping (bp)`),
      decimals = 0,
      use_seps = TRUE
    ) %>%
    gt::fmt_number(
      columns = `Coverage %`,
      decimals = 1
    ) %>%
    # Color code coverage
    gt::tab_style(
      style = list(
        gt::cell_fill(color = "#E8F5E9"),
        gt::cell_text(weight = "bold")
      ),
      locations = gt::cells_body(
        columns = `Coverage %`,
        rows = `Coverage %` > 80
      )
    ) %>%
    gt::tab_style(
      style = list(
        gt::cell_fill(color = "#FFF3E0")
      ),
      locations = gt::cells_body(
        columns = `Coverage %`,
        rows = `Coverage %` > 50 & `Coverage %` <= 80
      )
    ) %>%
    # Title
    gt::tab_header(
      title = "Genomic Context Coverage",
      subtitle = "v5.0q benchmark on GRCh38"
    )

  return(table)
}


# ============================================================================
# EXAMPLE 5: Benchmark Comparison (Multiple Series)
# ============================================================================

example_figure_benchmark_comparison <- function() {
  metrics <- load_genomic_context_metrics()

  plot <- metrics %>%
    filter(var_type == "smvar") %>%
    group_by(bench_version, ref) %>%
    summarise(
      avg_density = mean(variant_density_per_mb),
      total_variants = sum(total_variants),
      .groups = "drop"
    ) %>%
    ggplot(aes(
      x = ref,
      y = avg_density,
      fill = bench_version,
      color = bench_version
    )) +
    geom_col(position = "dodge", alpha = 0.8, linewidth = 0.4) +

    # Apply themes
    scale_benchmark_version(aesthetic = "fill") +
    scale_benchmark_version(aesthetic = "color", guide = "none") +
    theme_manuscript() +
    labs(
      title = "Mean Variant Density Across Benchmarks",
      x = "Reference Genome",
      y = "Mean Variants per Megabase",
      fill = "Benchmark",
      subtitle = "Small variants only"
    )

  return(plot)
}


# ============================================================================
# EXAMPLE 6: Comprehensive Multi-Panel Figure
# ============================================================================

example_figure_comprehensive <- function() {
  library(patchwork)

  metrics <- load_genomic_context_metrics()

  # Panel A: Variant counts
  p_a <- metrics %>%
    filter(ref == "GRCh38", var_type == "smvar") %>%
    group_by(bench_version, context_name) %>%
    summarise(total = sum(total_variants), .groups = "drop") %>%
    ggplot(aes(x = context_name, y = total, fill = bench_version)) +
    geom_col(position = "dodge") +
    scale_benchmark_version(aesthetic = "fill") +
    theme_manuscript() +
    labs(
      title = "A. Variant Counts",
      x = "Genomic Context",
      y = "Count",
      fill = "Benchmark"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # Panel B: Coverage percentage
  p_b <- metrics %>%
    filter(ref == "GRCh38") %>%
    group_by(bench_version, context_name) %>%
    summarise(pct = mean(pct_of_context), .groups = "drop") %>%
    ggplot(aes(x = context_name, y = pct, fill = bench_version)) +
    geom_col(position = "dodge") +
    scale_benchmark_version(aesthetic = "fill") +
    theme_manuscript() +
    labs(
      title = "B. Coverage %",
      x = "Genomic Context",
      y = "Percent (%)",
      fill = "Benchmark"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # Panel C: Variant density
  p_c <- metrics %>%
    filter(ref == "GRCh38") %>%
    group_by(bench_version, var_type) %>%
    summarise(density = mean(variant_density_per_mb), .groups = "drop") %>%
    ggplot(aes(x = var_type, y = density, fill = bench_version)) +
    geom_col(position = "dodge") +
    scale_benchmark_version(aesthetic = "fill") +
    theme_manuscript() +
    labs(
      title = "C. Mean Variant Density",
      x = "Variant Type",
      y = "Variants/Mb",
      fill = "Benchmark"
    )

  # Combine panels
  figure <- (p_a | p_b) / p_c +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  return(figure)
}


# ============================================================================
# Helper Function: Generate Standard Report
# ============================================================================

#' Generate Manuscript-Ready Report
#'
#' Creates a standardized report with consistent figures and tables
#'
#' @param output_dir Directory to save figures and tables
#' @param format File format(s) to export ("pdf", "png", or both)
#'
#' @keywords internal
#'
#' @export
generate_themed_report <- function(output_dir = "figures", format = "pdf") {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Generate and export figures
  message("Generating figures...")

  p1 <- example_figure_variant_counts()
  ggsave(
    file.path(output_dir, "Figure1_variant_counts.pdf"),
    p1,
    width = 5.5, height = 3.5, dpi = 300
  )
  message("✓ Figure 1 saved")

  p2 <- example_figure_variant_density()
  ggsave(
    file.path(output_dir, "Figure2_variant_density.pdf"),
    p2,
    width = 7, height = 4, dpi = 300
  )
  message("✓ Figure 2 saved")

  p3 <- example_figure_comprehensive()
  ggsave(
    file.path(output_dir, "Figure3_comprehensive.pdf"),
    p3,
    width = 7, height = 5, dpi = 300
  )
  message("✓ Figure 3 saved")

  # Generate and export tables
  message("\nGenerating tables...")

  t1 <- example_table_summary()
  gt::gtsave(t1, file.path(output_dir, "Table1_summary.pdf"))
  message("✓ Table 1 saved")

  t2 <- example_table_context_coverage()
  gt::gtsave(t2, file.path(output_dir, "Table2_coverage.pdf"))
  message("✓ Table 2 saved")

  message("\n✓ Report generation complete!")
  message(sprintf("Files saved to: %s", normalizePath(output_dir)))
}

# Example usage:
# generate_themed_report(output_dir = "manuscript/figures")
