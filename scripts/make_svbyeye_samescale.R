#!/usr/bin/env Rscript
# SVbyEye same-scale "sandwich" figure: HG002 paternal/maternal assembly vs
# reference, with benchmark-region and large-exclusion annotation overlays.
#
# Reconstruction note: this reimplements the "same-scale redesign" developed
# 2026-06-17 in a scratch worktree (bright-oak-sj5e, scratch/ideogram-explore/)
# that was lost during a repo cleanup on 2026-08-10/11 (see
# scripts/make_ideogram_heatmap.R header for the full recovery story). The
# design decisions below come from session notes preserved outside the repo,
# not from the original code, and some implementation details are new
# (original plotAVA()/addAnnotation() call patterns were not recoverable):
#   - Exclusion regrouping: the red "large excluded regions" track is the
#     union of segdups + satellites + tandem-repeats + flanks + gaps only.
#     Seven other exclusion categories (self-discrep, consecutive-svs,
#     dipcall-pav_discrep-{smvar,stvar}, dipcall-bugs-T2TACE,
#     HG002Q100-errors, HG002-mosaic, pav-inversions) describe
#     benchmarking-tool limitations or HG002-assembly-specific error calls,
#     not reference/assembly structure, and were deliberately dropped from
#     this track (2026-06-17 decision).
#   - Same-scale layout: every panel shares the same bp-per-inch. Full-width
#     rows are windowed to [0, chr1 length] regardless of the chromosome's
#     actual length (shorter chromosomes simply end before the row's right
#     edge); chromosomes shorter than half of chr1 are paired two-per-row,
#     each windowed to [0, chr1 length / 2].
#   - Each row is a 3-track "sandwich": HG002 paternal ribbon / reference
#     ideogram row (with benchmark + exclusion annotation) / HG002 maternal
#     ribbon, colored by alignment direction (forward/inverted).
#   - chrX is an exception: HG002 is male (XY), so chrX has no paternal-origin
#     homolog. Its panel is a 2-track ref+maternal sandwich; no paternal
#     ribbon is fabricated.
#
# AI Disclosure: Developed with assistance from Claude (Anthropic).
#
# Usage:
#   Rscript scripts/make_svbyeye_samescale.R [REF] [chroms]
#   REF: GRCh38 (default) | GRCh37 | CHM13v2.0
#   chroms: comma-separated (default: chr6,chr8,chr15,chr20,chrX -- the
#           2026-06-17 main-text subset). Requires PAFs already generated via
#           scripts/prep_svbyeye_pafs.sh and BEDs via
#           scripts/prep_svbyeye_beds.sh.
# Output: figures/svbyeye_main_<ref>.{pdf,png}

suppressPackageStartupMessages({
  library(SVbyEye)
  library(ggplot2)
  library(GenomicRanges)
  library(patchwork)
  library(here)
})

args <- commandArgs(trailingOnly = TRUE)
REF <- if (length(args) >= 1) args[[1]] else "GRCh38"
CHROMS <- if (length(args) >= 2) strsplit(args[[2]], ",")[[1]] else c("chr6", "chr8", "chr15", "chr20", "chrX")

data_root <- Sys.getenv("Q100_DATA_ROOT", unset = here::here())
svb_dir <- here::here("results", "svbyeye", REF)
figs_dir <- here::here("figures")

fai_path <- file.path(data_root, "resources", "references", paste0(REF, ".fa.gz.fai"))
fai <- read.table(fai_path, header = FALSE, sep = "\t", colClasses = c("character", "numeric", rep("NULL", 3)))
names(fai) <- c("chrom", "length")
if (!grepl("^chr", fai$chrom[1])) fai$chrom <- paste0("chr", fai$chrom)
chrom_lengths <- setNames(fai$length, fai$chrom)

chr1_len <- chrom_lengths[["chr1"]]
HALF <- chr1_len / 2

# Colorblind-safe direction palette (orange forward / blue inverted), per the
# 2026-06-17 publication-quality checklist -- default green/blue is not safe
# for deuteranopia.
DIRECTION_COLORS <- c("+" = "#E69F00", "-" = "#0072B2")
BENCH_COLOR <- c("Benchmark regions" = "#54278F")
EXCL_COLOR <- c("Large excluded regions" = "#B2182B")
EXCL_MIN_SIZE <- 1e4 # >=10kb filter applied here, at plot time (not in prep)

read_bench_excl <- function(chrom) {
  bench <- read.table(file.path(svb_dir, "v5_benchmark_all.bed"), col.names = c("chrom", "start", "end"))
  bench <- bench[bench$chrom == chrom, ]
  bench_gr <- GRanges(bench$chrom, IRanges(bench$start + 1L, bench$end), type = "Benchmark regions")

  excl <- read.table(file.path(svb_dir, "excl_large_all.bed"), col.names = c("chrom", "start", "end"))
  excl <- excl[excl$chrom == chrom & (excl$end - excl$start) >= EXCL_MIN_SIZE, ]
  excl_gr <- GRanges(excl$chrom, IRanges(excl$start + 1L, excl$end), type = "Large excluded regions")

  list(bench = bench_gr, excl = excl_gr)
}

# Build one chromosome's sandwich panel (3 rows: PAT/ref/MAT, or 2 rows for
# chrX which has no paternal homolog), windowed to a shared xlim so bp-per-
# inch matches every other panel at the same width class.
build_panel <- function(chrom, xlim_max, show_legend = FALSE) {
  mat_path <- file.path(svb_dir, chrom, "ref_mat.named.paf")
  pat_path <- file.path(svb_dir, chrom, "ref_pat.named.paf")
  has_pat <- file.exists(pat_path)

  mat <- readPaf(mat_path, include.paf.tags = FALSE)
  if (has_pat) {
    pat <- readPaf(pat_path, include.paf.tags = FALSE)
    paf <- rbind(mat, pat)
    order <- c("HG002_PAT", chrom, "HG002_MAT")
  } else {
    paf <- mat
    order <- c(chrom, "HG002_MAT")
  }

  p <- plotAVA(paf, seqnames.order = order, color.by = "direction") +
    scale_fill_manual(values = DIRECTION_COLORS, name = "Alignment\ndirection") +
    scale_color_manual(values = DIRECTION_COLORS, guide = "none")

  be <- read_bench_excl(chrom)
  p <- addAnnotation(p, annot.gr = be$bench, coordinate.space = "target", shape = "rectangle",
    fill.by = "type", color.palette = BENCH_COLOR, annotation.level = 0.05
  )
  p <- addAnnotation(p, annot.gr = be$excl, coordinate.space = "target", shape = "rectangle",
    fill.by = "type", color.palette = EXCL_COLOR, annotation.level = 0.12
  )

  p <- p + coord_cartesian(xlim = c(0, xlim_max), expand = FALSE) +
    scale_x_continuous(labels = scales::label_number(scale = 1e-6, suffix = " Mb")) +
    ggtitle(chrom) +
    theme(plot.title = element_text(size = 9, face = "bold"), axis.title.x = element_blank())

  if (!show_legend) p <- p + theme(legend.position = "none")
  p
}

# --- Panel layout: full-width rows for chroms >= half of chr1; paired
# half-width rows (two per row) for the rest, in the order given. ----------
is_full <- sapply(CHROMS, function(ch) chrom_lengths[[ch]] >= HALF)
full_chroms <- CHROMS[is_full]
half_chroms <- CHROMS[!is_full]

message("Full-width panels (xlim = 0..", format(chr1_len, big.mark = ","), "): ", paste(full_chroms, collapse = ", "))
message("Half-width panels (xlim = 0..", format(HALF, big.mark = ","), "): ", paste(half_chroms, collapse = ", "))

rows <- list()
legend_used <- FALSE
for (ch in full_chroms) {
  rows[[length(rows) + 1]] <- build_panel(ch, chr1_len, show_legend = !legend_used)
  legend_used <- TRUE
}
if (length(half_chroms) > 0) {
  half_panels <- lapply(half_chroms, function(ch) {
    p <- build_panel(ch, HALF, show_legend = !legend_used)
    legend_used <<- TRUE
    p
  })
  # Odd number of half-width chroms: pad the last row with a blank panel.
  if (length(half_panels) %% 2 == 1) half_panels[[length(half_panels) + 1]] <- patchwork::plot_spacer()
  for (i in seq(1, length(half_panels), by = 2)) {
    rows[[length(rows) + 1]] <- wrap_plots(half_panels[i:min(i + 1, length(half_panels))], nrow = 1)
  }
}

combined <- wrap_plots(rows, ncol = 1) +
  plot_annotation(title = sprintf("HG002 vs %s: assembly alignment, same scale", REF))

n_rows <- length(rows)
out_base <- file.path(figs_dir, paste0("svbyeye_main_", tolower(REF)))
message("Writing ", out_base, ".pdf / .png")
ggsave(paste0(out_base, ".pdf"), combined, width = 7, height = 1.6 * n_rows + 0.6, limitsize = FALSE)
ggsave(paste0(out_base, ".png"), combined, width = 7, height = 1.6 * n_rows + 0.6, dpi = 300, limitsize = FALSE)
message("Done.")
