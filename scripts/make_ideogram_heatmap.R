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
#   - Capture each chromosome as a vector grob with gridGraphics, then compose
#     rows with widths proportional to chromosome length. Rows are capped at
#     chr1's span so shorter chromosomes can share a row without changing scale.
#   - GRCh38: 3 tracks (v5.0q smvar, v4.2.1 smvar, v5.0q stvar). v0.6 stvar is
#     intentionally omitted here (GRCh38 has no native previous SV benchmark;
#     the 2026-06-17 session dropped it from this figure by request).
#   - GRCh37: 4 tracks, adding v0.6 stvar (GRCh37-native, no liftOver needed).
#   - Main-text panel: chr1 and chr8 heatmaps plus a 0-15 Mb SVbyEye view of
#     chr8 highlighting the large HG002 inversion excluded from the benchmark.
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
  library(grid)
  library(gridExtra)
  library(gridGraphics)
  library(SVbyEye)
  library(ggplot2)
  library(here)
})

BIN_SIZE <- 1e5 # 100 kb
COVERAGE_LOG_FLOOR <- 1e-4
RAMP <- colorRampPalette(c("#c6dbef", "#08306b"))(101)
ZERO_COL <- "grey90"

data_root <- Sys.getenv("Q100_DATA_ROOT", unset = here::here())
res_dir <- file.path(data_root, "resources")
bmk_dir <- file.path(res_dir, "benchmarksets")
svb_grch38_dir <- here::here("results", "svbyeye", "GRCh38")
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
main_text_chroms <- c("chr1", "chr8")


# --- Track definitions --------------------------------------------------------
# Each track = one benchmark BED rendered as a 100kb coverage-fraction heatmap.
grch38_tracks <- list(
  list(label = "v4.2.1 smvar", bed = file.path(bmk_dir, "v4.2.1_GRCh38_smvar_benchmark.bed"), chr_prefix = FALSE),
  list(label = "v5.0q smvar", bed = file.path(bmk_dir, "v5.0q_GRCh38_smvar_benchmark.bed"), chr_prefix = FALSE),
  list(label = "v5.0q stvar", bed = file.path(bmk_dir, "v5.0q_GRCh38_stvar_benchmark.bed"), chr_prefix = FALSE)
)

grch37_tracks <- list(
  list(label = "v4.2.1 smvar", bed = file.path(bmk_dir, "v4.2.1_GRCh37_smvar_benchmark.bed"), chr_prefix = TRUE),
  list(label = "v5.0q smvar", bed = file.path(bmk_dir, "v5.0q_GRCh37_smvar_benchmark.bed"), chr_prefix = TRUE),
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

coverage_to_color_scale <- function(frac) {
  pmax(
    0,
    pmin(
      1,
      -log10(pmax(1 - frac, COVERAGE_LOG_FLOOR)) /
        -log10(COVERAGE_LOG_FLOOR)
    )
  )
}

frac_to_color <- function(frac) {
  scaled <- coverage_to_color_scale(frac)
  ifelse(frac <= 0, ZERO_COL, RAMP[pmax(1, pmin(101, round(scaled * 100) + 1))])
}

make_bins <- function(chrom_lengths, chroms) {
  tileGenome(chrom_lengths[chroms], tilewidth = BIN_SIZE, cut.last.tile.in.chrom = TRUE)
}

# --- Grid-based vector layout -------------------------------------------------
LABEL_SPACE_BP <- 18e6
TRACK_LABEL_SPACE_BP <- 40e6

chromosome_panel_span <- function(chrom, chrom_lengths, last_chrom) {
  unname(chrom_lengths[chrom]) + LABEL_SPACE_BP +
    if (identical(chrom, last_chrom)) TRACK_LABEL_SPACE_BP else 0
}

make_plot_params <- function(n_tracks, leftmargin, rightmargin = 0.01) {
  pp <- getDefaultPlotParams(plot.type = 1)
  pp$data1height <- 12 * n_tracks + 4 * (n_tracks - 1)
  pp$ideogramheight <- 5
  pp$leftmargin <- leftmargin
  pp$rightmargin <- rightmargin
  pp$data1inmargin <- 2
  pp$data1outmargin <- 6
  pp$topmargin <- 5
  pp$bottommargin <- 3
  pp
}

prepare_tracks <- function(tracks) {
  lapply(tracks, function(trk) {
    trk$gr <- load_bed_gr(trk$bed, trk$chr_prefix)
    trk
  })
}

plot_chromosome_heatmap <- function(chrom, chrom_lengths, tracks, genome,
                                    show_track_labels = FALSE,
                                    track_label_side = "right",
                                    cex.label = 0.5, cex.chrom = 0.8) {
  n <- length(tracks)
  track_label_space <- if (show_track_labels) TRACK_LABEL_SPACE_BP else 0
  panel_span <- unname(chrom_lengths[chrom]) + LABEL_SPACE_BP + track_label_space
  left_space <- LABEL_SPACE_BP +
    if (show_track_labels && track_label_side == "left") TRACK_LABEL_SPACE_BP else 0
  right_space <-
    if (show_track_labels && track_label_side == "right") TRACK_LABEL_SPACE_BP else 0
  chromosome_labels <- function(karyoplot, ...) {
    kpAddChromosomeNames(
      karyoplot,
      cex = cex.chrom,
      font = 2,
      xoffset = 0.005
    )
  }
  kp <- plotKaryotype(
    genome = genome,
    chromosomes = chrom,
    plot.type = 1,
    plot.params = make_plot_params(
      n,
      leftmargin = left_space / panel_span,
      rightmargin = if (right_space > 0) right_space / panel_span else 0.01
    ),
    labels.plotter = chromosome_labels
  )
  gap <- 0.01
  track_h <- (1 - gap * (n - 1)) / n
  bins <- make_bins(chrom_lengths, chrom)

  for (i in seq_along(tracks)) {
    trk <- tracks[[i]]
    r0 <- (n - i) * (track_h + gap)
    r1 <- r0 + track_h
    frac <- bin_coverage_fraction(trk$gr, bins)
    kpDataBackground(kp, r0 = r0, r1 = r1, color = "white")
    kpPlotRegions(kp, bins, col = frac_to_color(frac), border = NA, r0 = r0, r1 = r1)
    if (show_track_labels) {
      kpAddLabels(
        kp,
        labels = trk$label,
        side = track_label_side,
        r0 = r0,
        r1 = r1,
        cex = cex.label,
        label.margin = 0.01
      )
    }
  }
  invisible(kp)
}

make_chromosome_grob <- function(chrom, chrom_lengths, tracks, genome,
                                 show_track_labels = FALSE,
                                 track_label_side = "right",
                                 cex.label = 0.5, cex.chrom = 0.8) {
  echoGrob(
    function() {
      par(mar = rep(0, 4))
      plot_chromosome_heatmap(
        chrom,
        chrom_lengths,
        tracks,
        genome,
        show_track_labels = show_track_labels,
        track_label_side = track_label_side,
        cex.label = cex.label,
        cex.chrom = cex.chrom
      )
    },
    prefix = paste0("ideogram-", chrom)
  )
}

pack_chromosome_rows <- function(chroms, chrom_lengths) {
  last_chrom <- tail(chroms, 1)
  max_span <- max(chrom_lengths[chroms]) + LABEL_SPACE_BP
  rows <- list()
  current <- character()
  current_span <- 0

  for (chrom in chroms) {
    chrom_span <- chromosome_panel_span(chrom, chrom_lengths, last_chrom)
    if (length(current) > 0 && current_span + chrom_span > max_span) {
      rows[[length(rows) + 1]] <- current
      current <- character()
      current_span <- 0
    }
    current <- c(current, chrom)
    current_span <- current_span + chrom_span
  }
  rows[[length(rows) + 1]] <- current
  rows
}

make_legend_grob <- function(fontsize = 8) {
  n <- length(RAMP)
  x_breaks <- seq(0.08, 0.92, length.out = n + 1)
  tick_coverage <- c(0, 0.9, 0.99, 0.999, 1)
  tick_x <- 0.08 + 0.84 * coverage_to_color_scale(tick_coverage)
  tick_labels <- c("0%", "90%", "99%", "99.9%", "100%")

  grobTree(
    textGrob(
      "Fraction covered per 100 kb bin\n(log scale toward 100%)",
      x = 0.5, y = 0.78,
      gp = gpar(fontsize = fontsize)
    ),
    rectGrob(
      x = x_breaks[-(n + 1)],
      y = 0.38,
      width = diff(x_breaks),
      height = 0.10,
      just = c("left", "center"),
      gp = gpar(fill = RAMP, col = NA)
    ),
    rectGrob(
      x = x_breaks[1],
      y = 0.38,
      width = diff(x_breaks)[1],
      height = 0.10,
      just = c("left", "center"),
      gp = gpar(fill = ZERO_COL, col = NA)
    ),
    segmentsGrob(
      x0 = tick_x, x1 = tick_x,
      y0 = 0.32, y1 = 0.34,
      gp = gpar(col = "black", lwd = 0.5)
    ),
    textGrob(
      tick_labels,
      x = tick_x, y = 0.24,
      just = c("centre", "top"),
      gp = gpar(fontsize = fontsize)
    )
  )
}

build_ideogram_figure <- function(chroms, chrom_lengths, tracks, genome,
                                  title = NULL, cex.label = 0.5,
                                  cex.chrom = 0.8) {
  rows <- pack_chromosome_rows(chroms, chrom_lengths)
  last_chrom <- tail(chroms, 1)
  max_span <- max(chrom_lengths[chroms]) + LABEL_SPACE_BP
  prepared_tracks <- prepare_tracks(tracks)

  row_loads <- vapply(rows, function(row) {
    sum(vapply(
      row,
      chromosome_panel_span,
      numeric(1),
      chrom_lengths = chrom_lengths,
      last_chrom = last_chrom
    ))
  }, numeric(1))
  legend_row <- which.max(max_span - row_loads)

  row_grobs <- lapply(seq_along(rows), function(row_index) {
    row <- rows[[row_index]]
    chromosome_grobs <- lapply(row, function(chrom) {
      make_chromosome_grob(
        chrom,
        chrom_lengths,
        prepared_tracks,
        genome,
        show_track_labels = identical(chrom, last_chrom),
        cex.label = cex.label,
        cex.chrom = cex.chrom
      )
    })
    chromosome_widths <- vapply(
      row,
      chromosome_panel_span,
      numeric(1),
      chrom_lengths = chrom_lengths,
      last_chrom = last_chrom
    )
    remaining_width <- max_span - sum(chromosome_widths)

    if (row_index == legend_row && remaining_width > 0) {
      chromosome_grobs <- c(
        chromosome_grobs,
        list(make_legend_grob(fontsize = if (length(chroms) > 10) 6.5 else 8))
      )
      chromosome_widths <- c(chromosome_widths, remaining_width)
    } else if (remaining_width > 0) {
      chromosome_grobs <- c(chromosome_grobs, list(nullGrob()))
      chromosome_widths <- c(chromosome_widths, remaining_width)
    }

    arrangeGrob(
      grobs = chromosome_grobs,
      nrow = 1,
      widths = unit(chromosome_widths, "null"),
      padding = unit(0, "pt")
    )
  })

  figure <- arrangeGrob(
    grobs = row_grobs,
    ncol = 1,
    heights = unit(rep(1, length(row_grobs)), "null"),
    padding = unit(0, "pt")
  )
  if (!is.null(title)) {
    figure <- arrangeGrob(
      textGrob(title, gp = gpar(fontsize = 11)),
      figure,
      ncol = 1,
      heights = unit.c(unit(0.28, "in"), unit(1, "null")),
      padding = unit(0, "pt")
    )
  }
  figure
}

make_chr8_svbyeye_grob <- function(xmax = 15e6) {
  direction_colors <- c("+" = "#E69F00", "-" = "#0072B2")
  annotation_colors <- c(
    "Benchmark regions" = "#54278F",
    "Large excluded regions" = "#B2182B"
  )
  paf <- rbind(
    readPaf(
      file.path(svb_grch38_dir, "chr8", "ref_mat.named.paf"),
      include.paf.tags = FALSE
    ),
    readPaf(
      file.path(svb_grch38_dir, "chr8", "ref_pat.named.paf"),
      include.paf.tags = FALSE
    )
  )

  bench <- read.table(
    file.path(svb_grch38_dir, "v5_benchmark_all.bed"),
    col.names = c("chrom", "start", "end")
  )
  bench <- bench[bench$chrom == "chr8" & bench$start < xmax, ]
  excl <- read.table(
    file.path(svb_grch38_dir, "excl_large_all.bed"),
    col.names = c("chrom", "start", "end")
  )
  excl <- excl[
    excl$chrom == "chr8" & excl$start < xmax & (excl$end - excl$start) >= 1e4,
  ]
  annotations <- c(
    GRanges(
      "chr8",
      IRanges(bench$start + 1L, pmin(bench$end, xmax)),
      type = "Benchmark regions",
      reference_track = "chr8"
    ),
    GRanges(
      "chr8",
      IRanges(excl$start + 1L, pmin(excl$end, xmax)),
      type = "Large excluded regions",
      reference_track = "chr8"
    )
  )

  plot <- plotAVA(
    paf,
    seqnames.order = c("HG002_PAT", "chr8", "HG002_MAT"),
    color.by = "direction"
  ) +
    scale_fill_manual(
      values = direction_colors,
      breaks = c("-", "+"),
      labels = c("Inverted", "Forward"),
      name = "Alignment direction"
    ) +
    scale_color_manual(values = direction_colors, guide = "none")
  plot <- addAnnotation(
    plot,
    annot.gr = annotations,
    coordinate.space = "target",
    shape = "rectangle",
    fill.by = "type",
    color.palette = annotation_colors,
    annotation.level = 0,
    y.label.id = "reference_track"
  )
  for (i in seq_along(plot$scales$scales)) {
    if (any(grepl("^(fill|colour|color)", plot$scales$scales[[i]]$aesthetics))) {
      plot$scales$scales[[i]]$guide <- "none"
    }
  }
  plot <- plot +
    coord_cartesian(xlim = c(0, xmax), expand = FALSE) +
    scale_x_continuous(
      breaks = seq(0, xmax, by = 5e6),
      labels = scales::label_number(scale = 1e-6, suffix = " Mb"),
      minor_breaks = NULL
    ) +
    labs(title = "chr8: 0-15 Mb", x = NULL) +
    theme(
      plot.title = element_text(size = 9, face = "bold"),
      axis.text = element_text(size = 7),
      legend.position = "none",
      plot.margin = margin(4, 14, 2, 5)
    )

  ggplotGrob(plot)
}

make_svbyeye_legend_grob <- function(fontsize = 7.5) {
  item_x <- c(0.31, 0.66)
  row_y <- c(0.67, 0.25)

  grobTree(
    textGrob(
      c("Region annotation", "Direction"),
      x = 0.02, y = row_y,
      just = c("left", "center"),
      gp = gpar(fontsize = fontsize, fontface = "bold")
    ),
    rectGrob(
      x = rep(item_x, 2),
      y = rep(row_y, each = 2),
      width = 0.035, height = 0.20,
      just = c("left", "center"),
      gp = gpar(
        fill = c("#54278F", "#B2182B", "#0072B2", "#E69F00"),
        col = NA
      )
    ),
    textGrob(
      c("Benchmark regions", "Large excluded regions", "Inverted", "Forward"),
      x = rep(item_x + 0.045, 2),
      y = rep(row_y, each = 2),
      just = c("left", "center"),
      gp = gpar(fontsize = fontsize)
    )
  )
}

add_panel_label <- function(grob, label) {
  grobTree(
    grob,
    textGrob(
      label,
      x = 0.008, y = 0.98,
      just = c("left", "top"),
      gp = gpar(fontsize = 13, fontface = "bold")
    )
  )
}

build_main_figure <- function() {
  prepared_tracks <- prepare_tracks(grch38_tracks)
  chr1_grob <- make_chromosome_grob(
    "chr1",
    grch38_lengths,
    prepared_tracks,
    genome = "hg38",
    show_track_labels = TRUE,
    track_label_side = "left",
    cex.label = 0.7,
    cex.chrom = 1
  )
  chr8_grob <- make_chromosome_grob(
    "chr8",
    grch38_lengths,
    prepared_tracks,
    genome = "hg38",
    cex.label = 0.7,
    cex.chrom = 1
  )
  # These widths align the chr8 chromosome body with chr1 while preserving the
  # same bp-per-inch scale for the two ideograms.
  second_row <- arrangeGrob(
    nullGrob(),
    chr8_grob,
    add_panel_label(make_chr8_svbyeye_grob(), "B"),
    nrow = 1,
    widths = unit(c(0.131, 0.529, 0.34), "null"),
    padding = unit(0, "pt")
  )
  legend_row <- arrangeGrob(
    make_legend_grob(fontsize = 7.5),
    make_svbyeye_legend_grob(fontsize = 7.5),
    nrow = 1,
    widths = unit(c(0.44, 0.56), "null"),
    padding = unit(0, "pt")
  )

  arrangeGrob(
    add_panel_label(chr1_grob, "A"),
    second_row,
    legend_row,
    ncol = 1,
    heights = unit(c(1, 1, 0.45), "null"),
    padding = unit(0, "pt")
  )
}

save_plot <- function(figure, base_path, width, height, png_res = 300) {
  pdf_path <- paste0(base_path, ".pdf")
  png_path <- paste0(base_path, ".png")
  message("Writing ", png_path, " and ", pdf_path)

  pdf(pdf_path, width = width, height = height, useDingbats = FALSE)
  grid.newpage()
  grid.draw(figure)
  dev.off()

  png(png_path, width = width * png_res, height = height * png_res, res = png_res)
  grid.newpage()
  grid.draw(figure)
  dev.off()

  message(sprintf("  PDF %s (%s bytes)", pdf_path, format(file.size(pdf_path), big.mark = ",")))
  message(sprintf("  PNG %s (%s bytes)", png_path, format(file.size(png_path), big.mark = ",")))
}

# --- Genome-wide supplemental figures -----------------------------------------
save_plot(
  build_ideogram_figure(
    autosomes,
    grch38_lengths,
    grch38_tracks,
    genome = "hg38",
    title = "HG002 Q100 Variant Benchmark - GRCh38 Region Coverage",
    cex.label = 0.35,
    cex.chrom = 0.75
  ),
  file.path(figs_dir, "ideogram_genomewide_grch38"),
  width = 7,
  height = 8.5
)

save_plot(
  build_ideogram_figure(
    autosomes,
    grch37_lengths,
    grch37_tracks,
    genome = "hg19",
    title = "HG002 Q100 Variant Benchmark - GRCh37 Region Coverage",
    cex.label = 0.35,
    cex.chrom = 0.75
  ),
  file.path(figs_dir, "ideogram_genomewide_grch37"),
  width = 7,
  height = 8.5
)

# --- Main-text figure ----------------------------------------------------------
# Panel A contains chr1 and chr8 ideograms. Panel B is a 0-15 Mb SVbyEye view
# of the large chr8 inversion. Both legends are collected in the bottom row.
save_plot(
  build_main_figure(),
  file.path(figs_dir, "ideogram_main"),
  width = 7,
  height = 4.6
)
