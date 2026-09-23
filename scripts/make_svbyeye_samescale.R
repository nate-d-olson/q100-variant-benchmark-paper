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
#   - Same-scale layout: every panel shares the same bp-per-inch, but each
#     coordinate axis ends at that chromosome's length. Panel widths are
#     proportional to chromosome length and each row is padded to chr1 length;
#     chromosomes shorter than half of chr1 are paired two-per-row.
#   - Each row is a 3-track "sandwich": HG002 paternal ribbon / reference
#     ideogram row / HG002 maternal ribbon, colored by alignment direction
#     (forward/inverted). Benchmark and exclusion annotations are centered on
#     the reference row rather than offset toward either HG002 assembly track.
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
# "main" (default) = the 2026-06-17 main-text subset -> svbyeye_main_<ref>;
# "supplement" (pass explicitly for the full 24-chromosome set) ->
# svbyeye_supplement_<ref>, so the two don't overwrite each other.
OUT_LABEL <- if (length(args) >= 3) args[[3]] else "main"

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
ANNOTATION_COLORS <- c(BENCH_COLOR, EXCL_COLOR)
LEGEND_COLORS <- c(
  "Forward alignment (+)" = unname(DIRECTION_COLORS[["+"]]),
  "Inverted alignment (-)" = unname(DIRECTION_COLORS[["-"]]),
  ANNOTATION_COLORS
)
EXCL_MIN_SIZE <- 1e4 # >=10kb filter applied here, at plot time (not in prep)

read_bench_excl <- function(chrom) {
  bench <- read.table(file.path(svb_dir, "v5_benchmark_all.bed"), col.names = c("chrom", "start", "end"))
  bench <- bench[bench$chrom == chrom, ]
  bench_gr <- GRanges(
    bench$chrom,
    IRanges(bench$start + 1L, bench$end),
    type = "Benchmark regions",
    reference_track = chrom
  )

  excl <- read.table(file.path(svb_dir, "excl_large_all.bed"), col.names = c("chrom", "start", "end"))
  excl <- excl[excl$chrom == chrom & (excl$end - excl$start) >= EXCL_MIN_SIZE, ]
  excl_gr <- GRanges(
    excl$chrom,
    IRanges(excl$start + 1L, excl$end),
    type = "Large excluded regions",
    reference_track = chrom
  )

  list(bench = bench_gr, excl = excl_gr)
}

# Build one chromosome's sandwich panel (3 rows: PAT/ref/MAT, or 2 rows for
# chrX/chrY which each have only one homolog in a male sample). The axis ends at
# the chromosome length; row-level width allocation supplies the common scale.
build_panel <- function(chrom) {
  mat_path <- file.path(svb_dir, chrom, "ref_mat.named.paf")
  pat_path <- file.path(svb_dir, chrom, "ref_pat.named.paf")
  has_mat <- file.exists(mat_path)
  has_pat <- file.exists(pat_path)

  if (has_mat && has_pat) {
    paf <- rbind(readPaf(mat_path, include.paf.tags = FALSE), readPaf(pat_path, include.paf.tags = FALSE))
    order <- c("HG002_PAT", chrom, "HG002_MAT")
  } else if (has_mat) {
    paf <- readPaf(mat_path, include.paf.tags = FALSE)
    order <- c(chrom, "HG002_MAT")
  } else if (has_pat) {
    paf <- readPaf(pat_path, include.paf.tags = FALSE)
    order <- c("HG002_PAT", chrom)
  } else {
    stop("No PAF found for ", chrom, " (neither mat nor pat)")
  }

  p <- plotAVA(paf, seqnames.order = order, color.by = "direction") +
    scale_fill_manual(values = DIRECTION_COLORS, name = "Alignment direction") +
    scale_color_manual(values = DIRECTION_COLORS, guide = "none")

  be <- read_bench_excl(chrom)
  annotations <- c(be$bench, be$excl)
  p <- addAnnotation(
    p,
    annot.gr = annotations,
    coordinate.space = "target",
    shape = "rectangle",
    fill.by = "type",
    color.palette = ANNOTATION_COLORS,
    annotation.level = 0,
    y.label.id = "reference_track"
  )

  # Suppress guides at the scale level because ggnewscale guides can survive
  # theme(legend.position = "none") inside nested patchwork layouts.
  for (i in seq_along(p$scales$scales)) {
    if (any(grepl("^(fill|colour|color)", p$scales$scales[[i]]$aesthetics))) {
      p$scales$scales[[i]]$guide <- "none"
    }
  }

  p + coord_cartesian(xlim = c(0, chrom_lengths[[chrom]]), expand = FALSE) +
    scale_x_continuous(
      breaks = seq(0, chrom_lengths[[chrom]], by = 50e6),
      labels = scales::label_number(scale = 1e-6, suffix = " Mb"),
      minor_breaks = NULL
    ) +
    ggtitle(chrom) +
    theme(
      plot.title = element_text(size = 9, face = "bold"),
      axis.title.x = element_blank(),
      legend.position = "none"
    )
}

# --- Panel layout: full-width rows for chroms >= half of chr1; paired
# half-width rows (two per row) for the rest, in the order given. ----------
is_full <- sapply(CHROMS, function(ch) chrom_lengths[[ch]] >= HALF)
full_chroms <- CHROMS[is_full]
half_chroms <- CHROMS[!is_full]

message("Full-width panels: ", paste(full_chroms, collapse = ", "))
message("Paired half-width panels: ", paste(half_chroms, collapse = ", "))

# Pad every row to chr1 length so a fixed genomic distance has the same physical
# width in every panel while each axis displays only its chromosome coordinates.
build_row <- function(chroms) {
  widths <- unname(chrom_lengths[chroms]) / chr1_len
  panels <- lapply(chroms, build_panel)
  remaining_width <- 1 - sum(widths)

  if (remaining_width > 0) {
    panels[[length(panels) + 1]] <- plot_spacer()
    widths <- c(widths, remaining_width)
  }

  wrap_plots(panels, nrow = 1, widths = widths)
}

rows <- lapply(full_chroms, build_row)
if (length(half_chroms) > 0) {
  half_row_starts <- seq(1, length(half_chroms), by = 2)
  rows <- c(rows, lapply(half_row_starts, function(i) {
    build_row(half_chroms[i:min(i + 1, length(half_chroms))])
  }))
}

# Extract one compact legend rather than collecting the ggnewscale-generated
# guides from every nested chromosome panel.
legend_source <- ggplot(
  data.frame(
    x = seq_along(LEGEND_COLORS),
    y = 1,
    category = factor(names(LEGEND_COLORS), levels = names(LEGEND_COLORS))
  ),
  aes(x, y, fill = category)
) +
  geom_tile() +
  scale_fill_manual(
    values = LEGEND_COLORS,
    breaks = names(LEGEND_COLORS),
    name = NULL,
    guide = guide_legend(nrow = 1, byrow = TRUE)
  ) +
  theme_void() +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.justification = "center",
    legend.text = element_text(size = 8),
    legend.key.width = grid::unit(0.8, "lines")
  )

legend_table <- ggplotGrob(legend_source)
legend_grob <- legend_table$grobs[[which(legend_table$layout$name == "guide-box-top")]]
legend_row <- wrap_elements(full = legend_grob)

combined <- wrap_plots(
  c(list(legend_row), rows),
  ncol = 1,
  heights = c(0.18, rep(1, length(rows)))
) +
  plot_annotation(title = sprintf("HG002 vs %s: assembly alignment, same scale", REF))

n_rows <- length(rows)
out_base <- file.path(figs_dir, paste0("svbyeye_", OUT_LABEL, "_", tolower(REF)))
message("Writing ", out_base, ".pdf / .png")
ggsave(paste0(out_base, ".pdf"), combined, width = 7, height = 1.6 * n_rows + 0.6, limitsize = FALSE)
ggsave(paste0(out_base, ".png"), combined, width = 7, height = 1.6 * n_rows + 0.6, dpi = 300, limitsize = FALSE)
message("Done.")
