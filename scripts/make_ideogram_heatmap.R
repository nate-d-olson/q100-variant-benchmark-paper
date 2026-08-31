#!/usr/bin/env Rscript
# Region-coverage heatmap ideogram figures for the Q100 variant benchmark manuscript.
#
# Reconstruction note: this reimplements the "region-coverage heatmap" redesign
# developed 2026-06-05 to 2026-06-17 in a scratch worktree (bright-oak-sj5e,
# scratch/ideogram-explore/) that was lost during a repo cleanup on 2026-08-10/11
# (worktree removal partially failed, but the gitignored scratch/ directory was
# never committed and is not recoverable). The design decisions below are taken
# from session notes preserved outside the repo, not from the original code:
#   - Drop per-variant density tracks and the HP+TR+SD+MAP "Difficult" union band
#     (uninformative solid band at genome scale; see prior make_ideogram.R).
#   - Drop log2 fold-change tracks (flat for smvar, noisy for stvar).
#   - Replace with per-benchmark-version 100kb region-coverage heatmaps: each
#     bin colored by the fraction of the bin covered by that benchmark's BED.
#   - Render all target chromosomes in a single plotKaryotype() call so they
#     share one bp scale (chr1 visually wider than chr22, not independently
#     rescaled per panel).
#   - GRCh38: 3 tracks (v5.0q smvar, v4.2.1 smvar, v5.0q stvar). v0.6 stvar is
#     intentionally omitted here (GRCh38 has no native previous SV benchmark;
#     the 2026-06-17 session dropped it from this figure by request).
#   - GRCh37: 4 tracks, adding v0.6 stvar (GRCh37-native, no liftOver needed).
#   - Main-text panel: subset to chr1, chr8, chr9, chr13, chr19 (2026-06-17
#     recommendation) with larger fonts and a simplified legend.
#
# AI Disclosure: Developed with assistance from Claude (Anthropic).
#
# Usage: Rscript scripts/make_ideogram_heatmap.R
# Set Q100_DATA_ROOT to point at a checkout with resources/ populated when
# running from a worktree that doesn't have pipeline outputs generated.
#
# Outputs:
#   figures/ideogram_main.{pdf,png}               - main-text subset (GRCh38)
#   figures/ideogram_genomewide_grch38.{pdf,png}  - supplemental, all autosomes
#   figures/ideogram_genomewide_grch37.{pdf,png}  - supplemental, all autosomes

suppressPackageStartupMessages({
  library(karyoploteR)
  library(GenomicRanges)
  library(here)
})

BIN_SIZE <- 1e5 # 100 kb
RAMP <- colorRampPalette(c("#c6dbef", "#08306b"))(101)
ZERO_COL <- "grey90"

data_root <- Sys.getenv("Q100_DATA_ROOT", unset = here::here())
res_dir <- file.path(data_root, "resources")
bmk_dir <- file.path(res_dir, "benchmarksets")
figs_dir <- here::here("figures")

# --- Chromosome sizes (from .fai, avoids a BSgenome/network dependency) ------
read_fai_lengths <- function(fai_path, chr_prefix) {
  fai <- read.table(fai_path, header = FALSE, sep = "\t", colClasses = c("character", "numeric", rep("NULL", 3)))
  names(fai) <- c("chrom", "length")
  if (chr_prefix && !grepl("^chr", fai$chrom[1])) fai$chrom <- paste0("chr", fai$chrom)
  setNames(fai$length, fai$chrom)
}

grch38_lengths <- read_fai_lengths(file.path(res_dir, "references", "GRCh38.fa.gz.fai"), chr_prefix = FALSE)
grch37_lengths <- read_fai_lengths(file.path(res_dir, "references", "GRCh37.fa.gz.fai"), chr_prefix = TRUE)

autosomes <- paste0("chr", 1:22)
main_text_chroms <- c("chr1", "chr8", "chr9", "chr13", "chr19")

# --- Track definitions --------------------------------------------------------
# Each track = one benchmark BED rendered as a 100kb coverage-fraction heatmap.
grch38_tracks <- list(
  list(label = "v5.0q smvar", bed = file.path(bmk_dir, "v5.0q_GRCh38_smvar_benchmark.bed"), chr_prefix = FALSE),
  list(label = "v4.2.1 smvar", bed = file.path(bmk_dir, "v4.2.1_GRCh38_smvar_benchmark.bed"), chr_prefix = FALSE),
  list(label = "v5.0q stvar", bed = file.path(bmk_dir, "v5.0q_GRCh38_stvar_benchmark.bed"), chr_prefix = FALSE)
)

grch37_tracks <- list(
  list(label = "v5.0q smvar", bed = file.path(bmk_dir, "v5.0q_GRCh37_smvar_benchmark.bed"), chr_prefix = TRUE),
  list(label = "v4.2.1 smvar", bed = file.path(bmk_dir, "v4.2.1_GRCh37_smvar_benchmark.bed"), chr_prefix = TRUE),
  list(label = "v5.0q stvar", bed = file.path(bmk_dir, "v5.0q_GRCh37_stvar_benchmark.bed"), chr_prefix = TRUE),
  list(label = "v0.6 stvar", bed = file.path(bmk_dir, "v0.6_GRCh37_stvar_benchmark.bed"), chr_prefix = TRUE)
)

load_bed_gr <- function(path, chr_prefix) {
  df <- read.table(
    path,
    header = FALSE, sep = "\t",
    col.names = c("chrom", "start", "end"),
    colClasses = c("character", "integer", "integer")
  )
  if (chr_prefix && !grepl("^chr", df$chrom[1])) df$chrom <- paste0("chr", df$chrom)
  GRanges(seqnames = df$chrom, ranges = IRanges(start = df$start + 1L, end = df$end))
}

# Fraction of each bin covered by bench_gr (0-1), vectorized over bins_gr.
bin_coverage_fraction <- function(bench_gr, bins_gr) {
  out <- rep(0, length(bins_gr))
  hits <- findOverlaps(bins_gr, bench_gr)
  if (length(hits) == 0) {
    return(out)
  }
  inter <- pintersect(bins_gr[queryHits(hits)], bench_gr[subjectHits(hits)])
  covered <- tapply(width(inter), queryHits(hits), sum)
  idx <- as.integer(names(covered))
  out[idx] <- pmin(1, covered / width(bins_gr)[idx])
  out
}

frac_to_color <- function(frac) {
  ifelse(frac <= 0, ZERO_COL, RAMP[pmax(1, pmin(101, round(frac * 100) + 1))])
}

make_bins <- function(chrom_lengths, chroms) {
  tileGenome(chrom_lengths[chroms], tilewidth = BIN_SIZE, cut.last.tile.in.chrom = TRUE)
}

# --- Plot params ---------------------------------------------------------------
# data1height/topmargin/bottommargin follow the 6/17 tuning notes: reduced
# margins so more chromosomes fit on a US-Letter page at the same bp scale.
make_plot_params <- function(n_tracks, leftmargin = 0.12) {
  pp <- getDefaultPlotParams(plot.type = 1)
  pp$data1height <- 12 * n_tracks + 4 * (n_tracks - 1)
  pp$ideogramheight <- 5
  pp$leftmargin <- leftmargin
  pp$data1inmargin <- 2
  pp$data1outmargin <- 10
  pp$topmargin <- 15
  pp$bottommargin <- 5
  pp
}

# Render one region-coverage-heatmap ideogram.
plot_heatmap_ideogram <- function(chroms, chrom_lengths, tracks, genome, cex.label = 0.5,
                                   leftmargin = 0.16, main = NULL) {
  n <- length(tracks)
  kp <- plotKaryotype(
    genome = genome, chromosomes = chroms, plot.type = 1,
    plot.params = make_plot_params(n, leftmargin), cex = 0.6, cex.main = 0.9,
    main = main
  )
  gap <- 0.01
  track_h <- (1 - gap * (n - 1)) / n
  for (i in seq_along(tracks)) {
    trk <- tracks[[i]]
    r0 <- (i - 1) * (track_h + gap)
    r1 <- r0 + track_h
    bench_gr <- load_bed_gr(trk$bed, trk$chr_prefix)
    bins <- make_bins(chrom_lengths, chroms)
    frac <- bin_coverage_fraction(bench_gr, bins)
    kpDataBackground(kp, r0 = r0, r1 = r1, color = "white")
    kpPlotRegions(kp, bins, col = frac_to_color(frac), border = NA, r0 = r0, r1 = r1)
    kpAddLabels(kp, labels = trk$label, r0 = r0, r1 = r1, cex = cex.label, label.margin = 0.02)
  }
  invisible(kp)
}

# Legend as its own standalone plot (NOT overlaid on the karyoplot device):
# plotKaryotype ignores par(fig=...)/par(mfrow) and always claims the full
# device (a documented gotcha from the original 2026-06 session notes), so a
# shared bottom legend on the same device gets painted over by the last
# chromosome row. The original scratch work worked around this the same way:
# render the plot and the legend as separate images and composite afterwards
# (there, with PIL; here, with ImageMagick, since no vector PDF stacking tool
# is installed -- the composited PDF is a raster wrap, same tradeoff as before).
draw_legend <- function() {
  par(mar = c(1.2, 3, 0.2, 3))
  plot.new()
  n <- length(RAMP)
  rect(seq(0, 1, length.out = n + 1)[-(n + 1)], 0.3, seq(0, 1, length.out = n + 1)[-1], 0.9,
    col = RAMP, border = NA
  )
  rect(0, 0.3, 0.02, 0.9, col = ZERO_COL, border = NA)
  text(0, 0.05, "0%", cex = 0.8, adj = c(0, 0))
  text(1, 0.05, "100%", cex = 0.8, adj = c(1, 0))
  text(0.5, 0.05, "Fraction of 100kb bin covered by benchmark region", cex = 0.8, adj = c(0.5, 0))
}

composite_with_legend <- function(main_png, legend_png, out_png, out_pdf) {
  status <- system2("magick", c(main_png, legend_png, "-append", out_png))
  if (status != 0) stop("ImageMagick composite failed for ", out_png)
  status <- system2("magick", c(out_png, out_pdf))
  if (status != 0) stop("ImageMagick PDF wrap failed for ", out_pdf)
}

save_plot <- function(plot_fn, base_path, width, height, png_res = 300, legend_height = 0.6) {
  pdf_path <- paste0(base_path, ".pdf")
  png_path <- paste0(base_path, ".png")
  main_tmp <- paste0(base_path, "_main_tmp.png")
  legend_tmp <- paste0(base_path, "_legend_tmp.png")

  png(main_tmp, width = width * png_res, height = (height - legend_height) * png_res, res = png_res)
  plot_fn()
  dev.off()

  png(legend_tmp, width = width * png_res, height = legend_height * png_res, res = png_res)
  draw_legend()
  dev.off()

  message("Writing ", png_path, " and ", pdf_path)
  composite_with_legend(main_tmp, legend_tmp, png_path, pdf_path)
  file.remove(main_tmp, legend_tmp)

  message(sprintf("  PDF %s (%s bytes)", pdf_path, format(file.size(pdf_path), big.mark = ",")))
  message(sprintf("  PNG %s (%s bytes)", png_path, format(file.size(png_path), big.mark = ",")))
}

# --- Genome-wide supplemental figures -------------------------------------
save_plot(
  function() {
    plot_heatmap_ideogram(autosomes, grch38_lengths, grch38_tracks,
      genome = "hg38", cex.label = 0.35, leftmargin = 0.12,
      main = "HG002 Q100 Variant Benchmark - GRCh38 Region Coverage"
    )
  },
  file.path(figs_dir, "ideogram_genomewide_grch38"),
  width = 7, height = 10
)

save_plot(
  function() {
    plot_heatmap_ideogram(autosomes, grch37_lengths, grch37_tracks,
      genome = "hg19", cex.label = 0.35, leftmargin = 0.12,
      main = "HG002 Q100 Variant Benchmark - GRCh37 Region Coverage"
    )
  },
  file.path(figs_dir, "ideogram_genomewide_grch37"),
  width = 7, height = 10
)

# --- Main-text subset figure -------------------------------------------------
# chr1 (broad baseline), chr8 (HG002 inversion / stvar contrast), chr9 (shared
# 9q12 heterochromatic gap), chr13 (acrocentric p-arm gap), chr19 (gene-dense,
# compact) -- per the 2026-06-17 recommendation.
save_plot(
  function() {
    plot_heatmap_ideogram(main_text_chroms, grch38_lengths, grch38_tracks,
      genome = "hg38", cex.label = 0.7, leftmargin = 0.2,
      main = "HG002 Q100 Variant Benchmark"
    )
  },
  file.path(figs_dir, "ideogram_main"),
  width = 7, height = 4.5
)
