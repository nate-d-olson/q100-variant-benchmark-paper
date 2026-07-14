#!/usr/bin/env Rscript
# Generate genome-wide karyotype ideogram figure for the Q100 variant benchmark manuscript.
# AI Disclosure: Developed with assistance from Claude (Anthropic).
#
# Produces a multi-track ideogram showing benchmark region coverage across chromosomes for:
#   - Difficult regions (HP+TR+SD+MAP union)
#   - v0.6 stvar (lifted from GRCh37 to GRCh38)
#   - v4.2.1 smvar
#   - v5.0q smvar-only regions
#   - v5.0q stvar-only regions
#   - v5.0q regions covered by both smvar and stvar
#
# Usage: Rscript scripts/make_ideogram.R
# Output: manuscript/figs/ideogram.pdf, manuscript/figs/ideogram.png

suppressPackageStartupMessages({
  library(karyoploteR)
  library(rtracklayer)
  library(GenomicRanges)
  library(here)
})

source(here::here("R/plot_themes.R"))

# --- Colors -------------------------------------------------------------------
# Version track colors match the bench_version palette from plot_themes.R.
# The v5.0q sub-band uses a purple monochromatic sequence so all three
# sub-tracks read as "v5.0q" without colliding with the version track colors.
col_v06 <- "#1B9E77" # teal  (bench_version palette)
col_v421 <- "#D95F02" # orange (bench_version palette)
col_smvar_only <- "#9E9AC8" # light purple  (v5.0q smvar-only)
col_stvar_only <- "#756BB1" # medium purple (v5.0q stvar-only)
col_both <- "#54278F" # dark purple   (v5.0q both)
col_difficult <- "#888888"

# Variant density window size (1 Mb is appropriate for single-chromosome view)
variant_window_size <- 1e6

# Track r0/r1 coordinates.
# Each benchmark band is split into a thin region-coverage strip (top) and a
# taller density subplot (bottom). The three v5.0q sub-tracks share one region
# strip; their two VCFs each get a separate density sub-band.
tracks <- list(
  # Difficult regions — region strip only (no variant VCF)
  difficult     = list(r0 = 0.01, r1 = 0.06, col = col_difficult,  label = "Difficult"),
  # v4.2.1 smvar — region strip + density
  v421          = list(r0 = 0.08, r1 = 0.13, col = col_v421,       label = "v4.2.1 smvar"),
  v421_dens     = list(r0 = 0.13, r1 = 0.30, col = col_v421),
  # v5.0q — three overlapping region sub-tracks, then per-VCF density bands
  smvar_only    = list(r0 = 0.32, r1 = 0.37, col = col_smvar_only, label = "v5.0q smvar only"),
  stvar_only    = list(r0 = 0.32, r1 = 0.37, col = col_stvar_only, label = "v5.0q stvar only"),
  both          = list(r0 = 0.32, r1 = 0.37, col = col_both,       label = "v5.0q both"),
  v5_smvar_dens = list(r0 = 0.37, r1 = 0.55, col = col_smvar_only),
  v5_stvar_dens = list(r0 = 0.55, r1 = 0.73, col = col_stvar_only),
  # v0.6 stvar — region strip + density
  v06           = list(r0 = 0.75, r1 = 0.80, col = col_v06,        label = "v0.6 stvar"),
  v06_dens      = list(r0 = 0.80, r1 = 1.00, col = col_v06)
)

# --- Paths --------------------------------------------------------------------
res_dir <- here::here("resources")
bmk_dir <- file.path(res_dir, "benchmarksets")
strat_dir <- file.path(res_dir, "stratifications")
chain_file <- file.path(res_dir, "hg19ToHg38.over.chain.gz")

path_v06 <- file.path(bmk_dir, "v0.6_GRCh37_stvar_benchmark.bed")
path_v421 <- file.path(bmk_dir, "v4.2.1_GRCh38_smvar_benchmark.bed")
path_smvar <- file.path(bmk_dir, "v5.0q_GRCh38_smvar_benchmark.bed")
path_stvar <- file.path(bmk_dir, "v5.0q_GRCh38_stvar_benchmark.bed")

path_v06_vcf      <- file.path(bmk_dir, "v0.6_GRCh37_stvar_benchmark.vcf.gz")
path_v421_vcf     <- file.path(bmk_dir, "v4.2.1_GRCh38_smvar_benchmark.vcf.gz")
path_v5_smvar_vcf <- file.path(bmk_dir, "v5.0q_GRCh38_smvar_benchmark.vcf.gz")
path_v5_stvar_vcf <- file.path(bmk_dir, "v5.0q_GRCh38_stvar_benchmark.vcf.gz")

bedtools <- "/opt/homebrew/bin/bedtools"

# --- Data preparation ---------------------------------------------------------
message("Preparing temporary BED files...")
tmp_dir <- tempdir()
smvar_only_bed <- file.path(tmp_dir, "smvar_only.bed")
stvar_only_bed <- file.path(tmp_dir, "stvar_only.bed")
v5_both_bed <- file.path(tmp_dir, "v5_both.bed")
difficult_bed <- file.path(tmp_dir, "difficult.bed")

# smvar-only: regions in smvar but not in stvar
cmd_smvar_only <- sprintf(
  "%s subtract -a %s -b %s > %s",
  bedtools,
  path_smvar,
  path_stvar,
  smvar_only_bed
)

# stvar-only: regions in stvar but not in smvar
cmd_stvar_only <- sprintf(
  "%s subtract -a %s -b %s > %s",
  bedtools,
  path_stvar,
  path_smvar,
  stvar_only_bed
)

# both: intersection of smvar and stvar
cmd_both <- sprintf(
  "%s intersect -a %s -b %s > %s",
  bedtools,
  path_smvar,
  path_stvar,
  v5_both_bed
)
# difficult regions: merge HP, TR, SD, MAP
# Use gunzip -c for cross-platform compatibility (macOS zcat appends .Z)
cmd_difficult <- sprintf(
  "{ gunzip -c %s; gunzip -c %s; gunzip -c %s; gunzip -c %s; } | sort -k1,1 -k2,2n | %s merge -i - > %s",
  # "{ gunzip -c %s; gunzip -c %s; } | sort -k1,1 -k2,2n | %s merge -i - > %s",
  file.path(strat_dir, "GRCh38_HP.bed.gz"),
  file.path(strat_dir, "GRCh38_TR.bed.gz"),
  file.path(strat_dir, "GRCh38_SD.bed.gz"),
  file.path(strat_dir, "GRCh38_MAP.bed.gz"),
  bedtools,
  difficult_bed
)

for (cmd in list(cmd_smvar_only, cmd_stvar_only, cmd_both, cmd_difficult)) {
  ret <- system(cmd)
  if (ret != 0) stop("bedtools command failed:\n", cmd)
}
message("Temporary BED files ready.")

# --- liftOver v0.6 GRCh37 -> GRCh38 -----------------------------------------
message("Running liftOver for v0.6 stvar (GRCh37 -> GRCh38)...")
# import.chain() does not support .gz files directly; decompress to temp file first
chain_unzipped <- file.path(tmp_dir, "hg19ToHg38.over.chain")
if (!file.exists(chain_unzipped)) {
  ret <- system(sprintf("gunzip -c %s > %s", chain_file, chain_unzipped))
  if (ret != 0) stop("Failed to decompress chain file: ", chain_file)
}
chain <- rtracklayer::import.chain(chain_unzipped)
# v0.6 BED has no 'chr' prefix (GRCh37 uses '1', '2', etc.) — add it
v06_raw <- read.table(
  path_v06,
  header = FALSE,
  sep = "\t",
  col.names = c("chrom", "start", "end"),
  colClasses = c("character", "integer", "integer")
)
# Ensure chr prefix for liftOver (hg19 uses chr-prefixed names in the chain)
if (!grepl("^chr", v06_raw$chrom[1])) {
  v06_raw$chrom <- paste0("chr", v06_raw$chrom)
}
v06_gr <- GRanges(
  seqnames = v06_raw$chrom,
  ranges = IRanges(start = v06_raw$start + 1L, end = v06_raw$end)
)
v06_lifted_list <- rtracklayer::liftOver(v06_gr, chain)
v06_hg38 <- unlist(v06_lifted_list)
message(sprintf("v0.6 stvar: %d intervals -> %d after liftOver", length(v06_gr), length(v06_hg38)))
# Merge overlapping intervals introduced by liftOver
v06_hg38 <- GenomicRanges::reduce(sort(v06_hg38))
message(sprintf("v0.6 stvar: %d intervals after merging overlaps", length(v06_hg38)))

# --- Load GRanges -------------------------------------------------------------
message("Loading GRanges for all tracks...")

load_bed_gr <- function(path) {
  df <- read.table(
    path,
    header = FALSE,
    sep = "\t",
    col.names = c("chrom", "start", "end"),
    colClasses = c("character", "integer", "integer")
  )
  GRanges(
    seqnames = df$chrom,
    ranges = IRanges(start = df$start + 1L, end = df$end)
  )
}

gr_v421 <- load_bed_gr(path_v421)
gr_smvar_only <- load_bed_gr(smvar_only_bed)
gr_stvar_only <- load_bed_gr(stvar_only_bed)
gr_both <- load_bed_gr(v5_both_bed)
gr_difficult <- load_bed_gr(difficult_bed)

# Filter to standard autosomes + chrX (hg38 names)
keep_chroms <- paste0("chr", c(1:22, "X", "Y"))

filter_chroms <- function(gr, chroms = keep_chroms) {
  gr[seqnames(gr) %in% chroms]
}

gr_v421 <- filter_chroms(gr_v421)
gr_smvar_only <- filter_chroms(gr_smvar_only)
gr_stvar_only <- filter_chroms(gr_stvar_only)
gr_both <- filter_chroms(gr_both)
gr_difficult <- filter_chroms(gr_difficult)
v06_hg38 <- filter_chroms(v06_hg38)

message("All GRanges loaded.")

# --- Load variant positions from VCF files ------------------------------------
message("Loading variant positions from VCF files...")

# Efficiently extract CHROM/POS for target chromosomes at the shell level,
# avoiding loading the full VCF into R memory.
load_vcf_positions <- function(path, chroms) {
  chrom_filter <- paste(sprintf("$1==\"%s\"", chroms), collapse = " || ")
  df <- read.table(
    pipe(sprintf(
      "gunzip -c %s | grep -v '^#' | awk '(%s)' | cut -f1,2",
      path, chrom_filter
    )),
    col.names    = c("chrom", "pos"),
    colClasses   = c("character", "integer"),
    sep          = "\t"
  )
  if (nrow(df) == 0) return(GRanges())
  GRanges(seqnames = df$chrom, ranges = IRanges(start = df$pos, width = 1))
}

gr_vars_v421     <- load_vcf_positions(path_v421_vcf,     keep_chroms)
gr_vars_v5_smvar <- load_vcf_positions(path_v5_smvar_vcf, keep_chroms)
gr_vars_v5_stvar <- load_vcf_positions(path_v5_stvar_vcf, keep_chroms)

message(sprintf("Variant counts on %s — v4.2.1: %d  v5.0q smvar: %d  v5.0q stvar: %d",
  paste(keep_chroms, collapse = ","),
  length(gr_vars_v421), length(gr_vars_v5_smvar), length(gr_vars_v5_stvar)))

# v0.6 VCF is GRCh37 (no chr prefix) — liftOver to GRCh38
v06_chr37_names <- sub("^chr", "", keep_chroms)
chrom_filter_37 <- paste(sprintf("$1==\"%s\"", v06_chr37_names), collapse = " || ")
v06_vcf_raw <- read.table(
  pipe(sprintf(
    "gunzip -c %s | grep -v '^#' | awk '(%s)' | cut -f1,2",
    path_v06_vcf, chrom_filter_37
  )),
  col.names    = c("chrom", "pos"),
  colClasses   = c("character", "integer"),
  sep          = "\t"
)
v06_vcf_raw$chrom <- paste0("chr", v06_vcf_raw$chrom)
v06_vars_gr37 <- GRanges(
  seqnames = v06_vcf_raw$chrom,
  ranges   = IRanges(start = v06_vcf_raw$pos, width = 1)
)
gr_vars_v06 <- unlist(rtracklayer::liftOver(v06_vars_gr37, chain))
message(sprintf("v0.6 stvar variants on %s: %d -> %d after liftOver",
  paste(keep_chroms, collapse = ","), length(v06_vars_gr37), length(gr_vars_v06)))

message("Variant positions loaded.")

# --- Plot helpers -------------------------------------------------------------

# Build plot params; topmargin/bottommargin are per-chromosome row margins
# (default 120/100) — reducing them compacts the inter-chromosome spacing.
make_plot_params <- function(leftmargin = 0.14) {
  pp <- getDefaultPlotParams(plot.type = 1)
  pp$data1height    <- 100
  pp$ideogramheight <- 7
  pp$leftmargin     <- leftmargin
  pp$data1inmargin  <- 0   # no gap between ideogram and data
  pp$data1outmargin <- 5   # reduce from default 20
  pp$topmargin      <- 15  # reduce from default 120
  pp$bottommargin   <- 5   # reduce from default 100
  pp
}

# Add all data tracks to a karyoplot object.
# wsize controls the density window size (use a smaller value for zoomed views).
add_tracks <- function(kp, wsize = variant_window_size) {
  kpDataBackground(kp, r0 = tracks$difficult$r0, r1 = tracks$difficult$r1, color = "white")
  kpPlotRegions(kp, gr_difficult, col = tracks$difficult$col, border = NA,
    r0 = tracks$difficult$r0, r1 = tracks$difficult$r1)

  kpDataBackground(kp, r0 = tracks$v421$r0, r1 = tracks$v421_dens$r1, color = "white")
  kpPlotRegions(kp, gr_v421, col = tracks$v421$col, border = NA,
    r0 = tracks$v421$r0, r1 = tracks$v421$r1)
  kpPlotDensity(kp, gr_vars_v421, window.size = wsize,
    col = tracks$v421$col, border = NA,
    r0 = tracks$v421_dens$r0, r1 = tracks$v421_dens$r1)

  kpDataBackground(kp, r0 = tracks$both$r0, r1 = tracks$v5_stvar_dens$r1, color = "white")
  kpPlotRegions(kp, gr_smvar_only, col = tracks$smvar_only$col, border = NA,
    r0 = tracks$smvar_only$r0, r1 = tracks$smvar_only$r1)
  kpPlotRegions(kp, gr_stvar_only, col = tracks$stvar_only$col, border = NA,
    r0 = tracks$stvar_only$r0, r1 = tracks$stvar_only$r1)
  kpPlotRegions(kp, gr_both, col = tracks$both$col, border = NA,
    r0 = tracks$both$r0, r1 = tracks$both$r1)
  kpPlotDensity(kp, gr_vars_v5_smvar, window.size = wsize,
    col = tracks$v5_smvar_dens$col, border = NA,
    r0 = tracks$v5_smvar_dens$r0, r1 = tracks$v5_smvar_dens$r1)
  kpPlotDensity(kp, gr_vars_v5_stvar, window.size = wsize,
    col = tracks$v5_stvar_dens$col, border = NA,
    r0 = tracks$v5_stvar_dens$r0, r1 = tracks$v5_stvar_dens$r1)

  kpDataBackground(kp, r0 = tracks$v06$r0, r1 = tracks$v06_dens$r1, color = "white")
  kpPlotRegions(kp, v06_hg38, col = tracks$v06$col, border = NA,
    r0 = tracks$v06$r0, r1 = tracks$v06$r1)
  kpPlotDensity(kp, gr_vars_v06, window.size = wsize,
    col = tracks$v06$col, border = NA,
    r0 = tracks$v06_dens$r0, r1 = tracks$v06_dens$r1)
}

# Add left-margin track labels to a karyoplot object.
add_labels <- function(kp, cex_main = 0.38, cex_sub = 0.30) {
  kpAddLabels(kp, labels = "Difficult",
    r0 = tracks$difficult$r0, r1 = tracks$difficult$r1,
    cex = cex_main, col = tracks$difficult$col, label.margin = 0.02)
  kpAddLabels(kp, labels = "v4.2.1 smvar",
    r0 = tracks$v421$r0, r1 = tracks$v421_dens$r1,
    cex = cex_main, col = tracks$v421$col, label.margin = 0.02)
  kpAddLabels(kp, labels = "v5.0q",
    r0 = tracks$both$r0, r1 = tracks$v5_stvar_dens$r1,
    cex = cex_main, col = col_both, label.margin = 0.02)
  kpAddLabels(kp, labels = "smvar",
    r0 = tracks$v5_smvar_dens$r0, r1 = tracks$v5_smvar_dens$r1,
    cex = cex_sub, col = col_smvar_only, label.margin = 0.02)
  kpAddLabels(kp, labels = "stvar",
    r0 = tracks$v5_stvar_dens$r0, r1 = tracks$v5_stvar_dens$r1,
    cex = cex_sub, col = col_stvar_only, label.margin = 0.02)
  kpAddLabels(kp, labels = "v0.6 stvar",
    r0 = tracks$v06$r0, r1 = tracks$v06_dens$r1,
    cex = cex_main, col = tracks$v06$col, label.margin = 0.02)
}

# --- Plot functions -----------------------------------------------------------

# Figure 1: all chromosomes, no track labels (labels belong on the inset only)
plot_main <- function() {
  kp <- plotKaryotype(
    genome = "hg38", chromosomes = keep_chroms, plot.type = 1,
    plot.params = make_plot_params(leftmargin = 0.10),
    cex = 0.5, cex.main = 0.8,
    main = "HG002 Q100 Variant Benchmark"
  )
  add_tracks(kp)
}

# Figure 2: chr8 full chromosome with track labels
plot_chr8_full <- function() {
  kp <- plotKaryotype(
    genome = "hg38", chromosomes = "chr8", plot.type = 1,
    plot.params = make_plot_params(leftmargin = 0.22),
    cex = 0.8
  )
  add_tracks(kp)
  add_labels(kp, cex_main = 0.7, cex_sub = 0.55)
}

# Figure 3: chr8 zoomed to 6–14 Mb with track labels
plot_chr8_zoom <- function() {
  zoom_region <- GRanges("chr8", IRanges(6e6, 14e6))
  kp <- plotKaryotype(
    genome = "hg38", chromosomes = "chr8", plot.type = 1,
    plot.params = make_plot_params(leftmargin = 0.22),
    zoom = zoom_region, cex = 0.8
  )
  add_tracks(kp, wsize = 5e5)  # 500 kb windows for the 8 Mb zoom region
  add_labels(kp, cex_main = 0.7, cex_sub = 0.55)
}

# --- Output -------------------------------------------------------------------
save_plot <- function(plot_fn, base_path, width, height, png_res = 300) {
  pdf_path <- paste0(base_path, ".pdf")
  png_path <- paste0(base_path, ".png")

  message("Writing ", pdf_path)
  pdf(pdf_path, width = width, height = height)
  plot_fn()
  dev.off()

  message("Writing ", png_path)
  png(png_path, width = width * png_res, height = height * png_res, res = png_res)
  plot_fn()
  dev.off()

  message(sprintf("  PDF %s (%s bytes)", pdf_path, format(file.size(pdf_path), big.mark = ",")))
  message(sprintf("  PNG %s (%s bytes)", png_path, format(file.size(png_path), big.mark = ",")))
}

figs_dir <- here::here("manuscript/figs")

# All chromosomes: tall portrait (25 rows)
save_plot(plot_main,      file.path(figs_dir, "ideogram_main"),      width = 7, height = 12)

# chr8 full: single-chromosome strip, landscape-ish
save_plot(plot_chr8_full, file.path(figs_dir, "ideogram_chr8"),      width = 7, height = 3)

# chr8 zoomed: same dimensions as full for easy alignment in layout software
save_plot(plot_chr8_zoom, file.path(figs_dir, "ideogram_chr8_zoom"), width = 7, height = 3)
