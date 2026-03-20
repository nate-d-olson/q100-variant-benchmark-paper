#!/usr/bin/env Rscript
# Export all manuscript tables to a Word document for copy-paste into the draft.
# AI Disclosure: Developed with assistance from Claude (Anthropic).
#
# Usage: Rscript scripts/export_tables_docx.R
# Output: manuscript/tables.docx

library(tidyverse)
library(here)
library(officer)
library(flextable)

source(here("R/data_loading.R"))
source(here("R/plot_themes.R"))

output_path <- here("manuscript/tables.docx")

# --- Data Loading -----------------------------------------------------------

message("Loading primary analysis data...")
analysis_data <- load_primary_analysis_data(
  include_variants = TRUE,
  include_reference_sizes = TRUE,
  include_hg002q100_size = TRUE,
  include_benchmark_regions = TRUE
)
list2env(analysis_data, envir = environment())

variants_df <- variants_df %>%
  filter(region_ids == "BMKREGIONS", is_pass)

message("Loading exclusion data...")
exclusion_metrics_tbl <- load_exclusion_metrics()
exclusion_interaction_tbl <- load_exclusion_interactions()

# Exclusion name mapping (from benchmark_exclusions.qmd)
exclusion_name_map <- c(
  "segdups" = "Segmental Duplications",
  "tandem-repeats" = "Tandem Repeats",
  "satellites" = "Satellites",
  "pav-inversions" = "PAV-inversions",
  "flanks" = "Flanks",
  "VDJ" = "VDJ",
  "HG002Q100-errors" = "HG002Q100v1.1 errors",
  "TSPY2-segdups" = "TSPY2-segdups",
  "consecutive-svs" = "Consecutive SVs",
  "gaps" = "Gaps",
  "dipcall-bugs-T2TACE" = "Dipcall-bugs T2TACE",
  "svs-and-simple-repeats" = "SV regions",
  "self-discrep" = "Self-Discrepancies",
  "HG002-mosaic" = "HG002-mosaic",
  "complex-var" = "Complex Variants"
)

rename_exclusion <- function(x) {
  dplyr::if_else(x %in% names(exclusion_name_map), exclusion_name_map[x], x)
}

rename_exclusion_combination <- function(x) {
  vapply(x, function(combo) {
    parts <- strsplit(combo, "\\|")[[1]]
    paste(rename_exclusion(parts), collapse = "|")
  }, character(1), USE.NAMES = FALSE)
}

exclusion_metrics_tbl <- exclusion_metrics_tbl %>%
  mutate(exclusions = rename_exclusion(exclusions))

exclusion_interaction_tbl <- exclusion_interaction_tbl %>%
  mutate(exclusion_combination = rename_exclusion_combination(exclusion_combination))

bench_type_labels <- get_bench_type_labels()

# --- Table Building ----------------------------------------------------------

## Benchmark summary data (from benchmarkset_characterization.qmd)
bench_region_stats_df <- bench_regions_df %>%
  filter(
    bench_version %in% c("v5.0q", "v4.2.1", "v0.6"),
    chrom %in% CHROM_LEVELS
  ) %>%
  group_by(bench_type, ref, bench_version) %>%
  summarise(total_bp = sum(interval_size), .groups = "drop") %>%
  left_join(total_ref_asm_bp, by = "ref") %>%
  mutate(
    pct_ref_cvg = 100 * (total_bp / total_asm_bp),
    pct_asm_cvg = 100 * (total_bp / hg002q100_size)
  ) %>%
  select(-total_asm_bp)

bench_var_cnt_df <- variants_df %>%
  count(bench_type, ref, bench_version, var_type)

bench_summary_df <- bench_var_cnt_df %>%
  pivot_wider(names_from = var_type, values_from = n, values_fill = 0L) %>%
  mutate(Total = rowSums(across(any_of(c("SNV", "INDEL", "INS", "DEL"))))) %>%
  full_join(bench_region_stats_df, by = c("bench_type", "ref", "bench_version")) %>%
  mutate(
    ref = std_references(ref),
    bench_version = std_bench_versions(bench_version),
    bench_mb = round(total_bp / 1e6, 0)
  ) %>%
  arrange(bench_type, ref, bench_version)

format_bench_table <- function(df, var_cols) {
  df %>%
    select(ref, bench_version, all_of(var_cols), Total, bench_mb, pct_ref_cvg, pct_asm_cvg) %>%
    as_flextable(show_coltype = FALSE) %>%
    fmt_integer_flextable(columns = c(var_cols, "Total", "bench_mb")) %>%
    fmt_number_flextable(columns = c("pct_ref_cvg", "pct_asm_cvg"), decimals = 1) %>%
    cols_label_flextable(
      ref = "Reference",
      bench_version = "Version",
      Total = "Total Variants",
      bench_mb = "Bench (Mb)",
      pct_ref_cvg = "% Ref",
      pct_asm_cvg = "% Asm"
    ) %>%
    flextable::align(j = "ref", align = "right", part = "body") %>%
    theme_flextable_manuscript()
}

## Table 1a: Small Variants (GRCh38)
table1_smvar <- bench_summary_df %>%
  filter(bench_type == "smvar", ref == "GRCh38") %>%
  format_bench_table(var_cols = c("SNV", "INDEL"))

## Table 1b: Structural Variants (GRCh37)
table1_stvar <- bench_summary_df %>%
  filter(bench_type == "stvar", ref == "GRCh37") %>%
  format_bench_table(var_cols = c("DEL", "INS"))

## Table S: Small Variants (All refs)
table_s_smvar <- bench_summary_df %>%
  filter(bench_type == "smvar") %>%
  format_bench_table(var_cols = c("SNV", "INDEL"))

## Table S: Structural Variants (All refs)
table_s_stvar <- bench_summary_df %>%
  filter(bench_type == "stvar") %>%
  format_bench_table(var_cols = c("DEL", "INS"))

## Exclusion table: GRCh38
grch38_exclusion_tbl <- exclusion_metrics_tbl %>%
  filter(ref == "GRCh38") %>%
  select(exclusions, bench_type, intersect_bp, pct_of_dip) %>%
  group_by(exclusions) %>%
  mutate(is_smvar_only = all(bench_type == "smvar")) %>%
  ungroup() %>%
  distinct(exclusions, .keep_all = TRUE) %>%
  arrange(desc(intersect_bp)) %>%
  select(exclusions, intersect_bp, pct_of_dip, is_smvar_only)

smvar_only_rows <- which(grch38_exclusion_tbl$is_smvar_only)

table_excl_grch38 <- grch38_exclusion_tbl %>%
  select(-is_smvar_only) %>%
  flextable() %>%
  fmt_integer_flextable(columns = "intersect_bp") %>%
  fmt_number_flextable(columns = "pct_of_dip", decimals = 3) %>%
  flextable::align(j = c("intersect_bp", "pct_of_dip"), align = "right", part = "all") %>%
  flextable::italic(i = smvar_only_rows, j = "exclusions", part = "body") %>%
  flextable::color(i = smvar_only_rows, j = "exclusions", color = "#555555", part = "body") %>%
  cols_label_flextable(
    exclusions = "Exclusion",
    intersect_bp = "Bases Excluded",
    pct_of_dip = "% of Dip.bed"
  ) %>%
  theme_flextable_manuscript()

## Exclusion table: All References (Supplemental)
supp_exclusion_tbl <- exclusion_metrics_tbl %>%
  select(ref, exclusions, bench_type, intersect_bp, pct_of_dip) %>%
  group_by(exclusions) %>%
  mutate(is_smvar_only = all(bench_type == "smvar")) %>%
  ungroup() %>%
  distinct(exclusions, ref, .keep_all = TRUE) %>%
  select(-bench_type) %>%
  pivot_wider(names_from = ref, values_from = c(intersect_bp, pct_of_dip)) %>%
  arrange(desc(intersect_bp_GRCh38))

supp_smvar_only_rows <- which(supp_exclusion_tbl$is_smvar_only)

bp_cols <- c("intersect_bp_GRCh37", "intersect_bp_GRCh38", "intersect_bp_CHM13v2.0")
pct_cols <- c("pct_of_dip_GRCh37", "pct_of_dip_GRCh38", "pct_of_dip_CHM13v2.0")

table_excl_all_refs <- supp_exclusion_tbl %>%
  select(-is_smvar_only) %>%
  flextable() %>%
  fmt_integer_flextable(columns = bp_cols) %>%
  fmt_number_flextable(columns = pct_cols, decimals = 3) %>%
  flextable::align(j = c(bp_cols, pct_cols), align = "right", part = "all") %>%
  flextable::italic(i = supp_smvar_only_rows, j = "exclusions", part = "body") %>%
  flextable::color(i = supp_smvar_only_rows, j = "exclusions", color = "#555555", part = "body") %>%
  cols_label_flextable(
    exclusions = "Exclusion",
    intersect_bp_GRCh37 = "Bases",
    intersect_bp_GRCh38 = "Bases",
    intersect_bp_CHM13v2.0 = "Bases",
    pct_of_dip_GRCh37 = "% Dip",
    pct_of_dip_GRCh38 = "% Dip",
    pct_of_dip_CHM13v2.0 = "% Dip"
  ) %>%
  flextable::add_header_row(
    values = c("", "GRCh37", "GRCh38", "CHM13v2.0", "GRCh37", "GRCh38", "CHM13v2.0"),
    top = TRUE
  ) %>%
  flextable::merge_h(part = "header") %>%
  theme_flextable_manuscript()

## Exclusion interaction summary
interaction_summary_tbl <- exclusion_interaction_tbl %>%
  group_by(ref, bench_type) %>%
  mutate(
    total_bases = sum(bases_bp),
    total_variants = sum(variant_count),
    pct_bases = 100 * bases_bp / total_bases,
    pct_variants = 100 * variant_count / total_variants
  ) %>%
  ungroup()

## Excluded variant counts by exclusion (GRCh38)
table_excl_variants <- interaction_summary_tbl %>%
  filter(n_exclusions == 1, ref == "GRCh38") %>%
  select(exclusion_combination, bench_type, variant_count, bases_bp) %>%
  mutate(
    bench_type = bench_type_labels[bench_type],
    variants_per_mb = round(variant_count / (bases_bp / 1e6), 1)
  ) %>%
  arrange(desc(variant_count)) %>%
  flextable() %>%
  fmt_integer_flextable(columns = c("variant_count", "bases_bp")) %>%
  fmt_number_flextable(columns = "variants_per_mb", decimals = 1) %>%
  flextable::align(
    j = c("variant_count", "bases_bp", "variants_per_mb"),
    align = "right", part = "all"
  ) %>%
  cols_label_flextable(
    exclusion_combination = "Exclusion",
    bench_type = "Type",
    variant_count = "Variants",
    bases_bp = "Bases",
    variants_per_mb = "Variants/Mb"
  ) %>%
  theme_flextable_manuscript()

## Excluded variant counts - all refs (Supplemental)
table_excl_variants_all <- interaction_summary_tbl %>%
  filter(n_exclusions == 1) %>%
  select(ref, exclusion_combination, bench_type, variant_count) %>%
  mutate(bench_type = bench_type_labels[bench_type]) %>%
  pivot_wider(names_from = ref, values_from = variant_count, values_fill = 0L) %>%
  arrange(desc(GRCh38)) %>%
  flextable() %>%
  fmt_integer_flextable(columns = c("GRCh37", "GRCh38", "CHM13v2.0")) %>%
  flextable::align(
    j = c("GRCh37", "GRCh38", "CHM13v2.0"), align = "right", part = "all"
  ) %>%
  cols_label_flextable(
    exclusion_combination = "Exclusion",
    bench_type = "Type"
  ) %>%
  theme_flextable_manuscript()

## Exclusion summary statistics
table_excl_summary <- interaction_summary_tbl %>%
  group_by(ref, bench_type) %>%
  summarise(
    total_combinations = n(),
    single_exclusions = sum(n_exclusions == 1),
    multi_exclusions = sum(n_exclusions >= 2),
    max_overlap = max(n_exclusions),
    total_bases = sum(bases_bp),
    total_variants = sum(variant_count),
    .groups = "drop"
  ) %>%
  mutate(bench_type = bench_type_labels[bench_type]) %>%
  flextable() %>%
  fmt_integer_flextable(
    columns = c(
      "total_combinations", "single_exclusions", "multi_exclusions",
      "total_bases", "total_variants"
    )
  ) %>%
  flextable::align(
    j = c(
      "total_combinations", "single_exclusions", "multi_exclusions",
      "max_overlap", "total_bases", "total_variants"
    ),
    align = "right", part = "all"
  ) %>%
  cols_label_flextable(
    ref = "Reference",
    bench_type = "Type",
    total_combinations = "Total Combos",
    single_exclusions = "Single Only",
    multi_exclusions = "Multi-Overlap",
    max_overlap = "Max Overlap",
    total_bases = "Total Bases",
    total_variants = "Total Variants"
  ) %>%
  theme_flextable_manuscript()

# --- Build Word Document -----------------------------------------------------

message("Building Word document...")

# --- Table Legends ------------------------------------------------------------

legend_table1a <- paste(
  "Table 1a. Summary of the v5.0q and v4.2.1 small variant benchmark sets for HG002 on GRCh38.",
  "Variant counts include single-nucleotide variants (SNVs) and insertions/deletions (indels) less than 50 bp",
  "within benchmark regions, restricted to autosomes and sex chromosomes.",
  "Bench (Mb) indicates the total size of the benchmark regions in megabases.",
  "% Ref, percentage of the assembled reference genome covered by benchmark regions;",
  "% Asm, percentage of the HG002 Q100 diploid assembly covered by benchmark regions.",
  "See also Table S1."
)

legend_table1b <- paste(
  "Table 1b. Summary of the v5.0q and v0.6 structural variant (SV) benchmark sets for HG002 on GRCh37.",
  "Variant counts include deletions (DEL) and insertions (INS) of 50 bp or larger within benchmark regions.",
  "The v0.6 SV benchmark was only generated for GRCh37.",
  "Bench (Mb), % Ref, and % Asm are defined as in Table 1a.",
  "See also Table S2."
)

legend_table2 <- paste(
  "Table 2. Bases removed by each exclusion category from the GRCh38 v5.0q benchmark regions.",
  "Exclusion regions were subtracted from the dipcall assembly alignment regions (dip.bed) to define the final benchmark regions.",
  "Bases Excluded, number of base pairs in the intersection of each exclusion with dip.bed;",
  "% of Dip.bed, percentage of dip.bed removed by each exclusion.",
  "Exclusion names in italics apply only to the small variant benchmark.",
  "Individual exclusion regions may overlap; the total bases excluded is less than the sum of individual exclusions.",
  "See also Table S3."
)

legend_table3 <- paste(
  "Table 3. Number of variants in regions removed by each individual exclusion category for GRCh38.",
  "Variants were counted within single-exclusion regions (i.e., regions attributed to exactly one exclusion category)",
  "for both small variants (<50 bp) and structural variants (>=50 bp).",
  "Variants/Mb, variant density per megabase of excluded region.",
  "See also Table S4."
)

legend_table_s1 <- paste(
  "Table S1. Summary of small variant benchmark sets for HG002 across three reference genomes (GRCh37, GRCh38, CHM13v2.0).",
  "Variant counts include single-nucleotide variants (SNVs) and insertions/deletions (indels) less than 50 bp",
  "within benchmark regions.",
  "Bench (Mb), total size of benchmark regions in megabases;",
  "% Ref, percentage of the assembled reference genome covered;",
  "% Asm, percentage of the HG002 Q100 diploid assembly covered.",
  "The v4.2.1 benchmark was generated for GRCh37 and GRCh38 only."
)

legend_table_s2 <- paste(
  "Table S2. Summary of structural variant (SV) benchmark sets for HG002 across three reference genomes.",
  "Variant counts include deletions (DEL) and insertions (INS) of 50 bp or larger within benchmark regions.",
  "The v0.6 SV benchmark was only generated for GRCh37.",
  "Bench (Mb), % Ref, and % Asm are defined as in Table S1."
)

legend_table_s3 <- paste(
  "Table S3. Bases removed by each exclusion category from v5.0q benchmark regions across all three reference genomes.",
  "For each exclusion, the table shows the number of base pairs in the intersection with dip.bed (Bases)",
  "and the percentage of dip.bed removed (% Dip).",
  "Exclusion names in italics apply only to the small variant benchmark.",
  "Not all exclusion categories are applicable to every reference genome."
)

legend_table_s4 <- paste(
  "Table S4. Number of variants in single-exclusion regions across all three reference genomes.",
  "Variant counts are shown separately for small variants and structural variants.",
  "Only regions attributed to exactly one exclusion category are included."
)

legend_table_s5 <- paste(
  "Table S5. Summary of exclusion region interactions across reference genomes and variant types.",
  "Total Combos, total number of distinct exclusion combinations observed;",
  "Single Only, number of exclusion categories that appear without overlap;",
  "Multi-Overlap, number of combinations involving two or more exclusion categories;",
  "Max Overlap, maximum number of exclusion categories overlapping at any single region;",
  "Total Bases, total base pairs across all exclusion combinations;",
  "Total Variants, total variants across all exclusion combinations."
)

# --- Build Word Document -----------------------------------------------------

doc <- read_docx() %>%
  # Table 1a
  body_add_par("Table 1a: Small Variant Benchmark Summary (GRCh38)", style = "heading 2") %>%
  body_add_par(legend_table1a, style = "Normal") %>%
  body_add_flextable(table1_smvar) %>%
  body_add_break() %>%
  # Table 1b
  body_add_par("Table 1b: Structural Variant Benchmark Summary (GRCh37)", style = "heading 2") %>%
  body_add_par(legend_table1b, style = "Normal") %>%
  body_add_flextable(table1_stvar) %>%
  body_add_break() %>%
  # Exclusion table (GRCh38)
  body_add_par("Table 2: Bases Removed by Exclusion (GRCh38)", style = "heading 2") %>%
  body_add_par(legend_table2, style = "Normal") %>%
  body_add_flextable(table_excl_grch38) %>%
  body_add_break() %>%
  # Excluded variant counts (GRCh38)
  body_add_par("Table 3: Excluded Variant Counts by Exclusion (GRCh38)", style = "heading 2") %>%
  body_add_par(legend_table3, style = "Normal") %>%
  body_add_flextable(table_excl_variants) %>%
  body_add_break() %>%
  # --- Supplemental Tables ---
  body_add_par("Supplemental Tables", style = "heading 1") %>%
  # Table S1
  body_add_par("Table S1: Small Variant Benchmark Summary (All References)", style = "heading 2") %>%
  body_add_par(legend_table_s1, style = "Normal") %>%
  body_add_flextable(table_s_smvar) %>%
  body_add_break() %>%
  # Table S2
  body_add_par("Table S2: Structural Variant Benchmark Summary (All References)", style = "heading 2") %>%
  body_add_par(legend_table_s2, style = "Normal") %>%
  body_add_flextable(table_s_stvar) %>%
  body_add_break() %>%
  # Table S3
  body_add_par("Table S3: Bases Removed by Exclusion (All References)", style = "heading 2") %>%
  body_add_par(legend_table_s3, style = "Normal") %>%
  body_add_flextable(table_excl_all_refs) %>%
  body_add_break() %>%
  # Table S4
  body_add_par("Table S4: Excluded Variant Counts (All References)", style = "heading 2") %>%
  body_add_par(legend_table_s4, style = "Normal") %>%
  body_add_flextable(table_excl_variants_all) %>%
  body_add_break() %>%
  # Table S5
  body_add_par("Table S5: Exclusion Interaction Summary", style = "heading 2") %>%
  body_add_par(legend_table_s5, style = "Normal") %>%
  body_add_flextable(table_excl_summary)

print(doc, target = output_path)
message("Tables written to: ", output_path)
