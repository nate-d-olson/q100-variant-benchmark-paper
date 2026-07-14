#!/usr/bin/env python3
# AI Disclosure: This script was developed with assistance from Claude (Anthropic)
"""
make_chr8_figure.py — multi-panel Chr8 synteny figure using plotsr
=============================================================

Generates a two-panel figure comparing GRCh38 chr8 to HG002 maternal (MAT)
and paternal (PAT) haplotypes using plotsr for synteny visualization.

Panel A  Full Chr8 synteny: REF <-> MAT <-> PAT with structural annotations
         and markers for the PAV-excluded inversion region.
Panel B  Zoomed view of the largest inversion on the PAT haplotype, showing
         the region excluded from the v5.0q benchmark (chr8:8.2-12.2 Mb).

Track colours: REF (blue), MAT (green), PAT (red). Ribbon colours follow plotsr
defaults: grey = syntenic, orange = inversion, light green = translocation,
cyan = duplication. The plotsr legend is cropped from both panels; describe
these colours in the figure caption.

Dependencies (conda env: plotsr.yaml):
  plotsr, syri, minimap2, matplotlib, Pillow

Snakemake usage
---------------
  Called by workflow/rules/chr8_synteny.smk via the chr8_make_figure rule.
  Inversion coordinates are pre-computed by the chr8_find_inversion rule
  (find_chr8_inversion.py) and passed via --coords JSON file.

  Required SyRI inputs (consecutive genome pairs for plotsr):
    --rm  ref_matsyri.out  (REF <-> MAT)
    --mp  mat_patsyri.out  (MAT <-> PAT)

  The PAV exclusion region (chr8:8,237,843-12,234,345) is hardcoded based
  on the PAV callset inversion coordinates.
"""

import argparse
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec


# ---------------------------------------------------------------------------
# CONFIG — adjust these to tune the figure appearance
#
# To rerun after editing:
#   mamba activate q100-smk
#   snakemake --sdm conda --cores 4 chr8_synteny
#
# To check readability at true manuscript size, open the PDF (not PNG):
#   open results/chr8_synteny/chr8_figure.pdf
# Then View → Actual Size in Preview. PDFs carry physical dimensions so
# Preview scales them correctly, whereas PNG is shown at screen pixels.
# ---------------------------------------------------------------------------

# --- Font sizes (points) ---
PLOTSR_FONTSIZE = 14  # plotsr axis labels and tick labels (both panels)
PANEL_LABEL_FONTSIZE = 18  # bold A / B panel labels
ZOOM_TITLE_FONTSIZE = 14  # "Zoomed: REF … | Excluded: …" subtitle above Panel B

# --- plotsr canvas dimensions (inches) passed to plotsr before compositing ---
# Taller panels give plotsr more room and reduce the legend crop fraction needed.
# Width is kept fixed; only height matters here — the final figure width is set
# by FIG_WIDTH_INCHES below.
PANEL_A_HEIGHT = 2  # full chr8 overview panel
PANEL_B_HEIGHT = 2  # zoomed inversion panel

# --- Legend crop ---
# plotsr draws its legend at the top of every rendered PNG; we crop it off.
# Each value is the fraction of that panel's pixel height to remove from the top.
# Increase if a fragment of the legend is still visible; decrease if the
# chromosome y-axis label ("chr8") gets clipped.
PANEL_A_CROP = 0.0
PANEL_B_CROP = 0.0

# --- Gridlines ---
# Fraction of the axis height that each shortened vertical dashed gridline spans
# (centered on the chromosome tracks). 0.0 = invisible, 1.0 = full axis height.
GRIDLINE_SPAN = 1.0

# --- Intermediate panel rasterization DPI ---
# plotsr panels are rasterized from memory at this DPI before being embedded in
# the composite figure. Higher = sharper chromosome tracks in the final outputs.
# 600 is good for print; 300 is fine for screen-only use and renders faster.
# This is independent of the final PNG DPI (controlled by --dpi, default 300).
INTERMEDIATE_DPI = 300

# --- Assembled figure dimensions (inches) ---
# FIG_WIDTH_INCHES should match the target journal column width (e.g. 7.0 for a
# two-column figure, 3.5 for one-column). At 300 DPI, 7 in = 2100 px wide.
FIG_WIDTH_INCHES = 7
FIG_HEIGHT_INCHES = 4


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--ref", required=True, help="REF chromosome-length (.cl) or FASTA file"
    )
    p.add_argument(
        "--mat", required=True, help="MAT chromosome-length (.cl) or FASTA file"
    )
    p.add_argument(
        "--pat", required=True, help="PAT chromosome-length (.cl) or FASTA file"
    )
    p.add_argument(
        "--rm",
        required=True,
        dest="syri_rm",
        help="SyRI output: REF vs MAT (ref_matsyri.out)",
    )
    p.add_argument(
        "--mp",
        required=True,
        dest="syri_mp",
        help="SyRI output: MAT vs PAT (mat_patsyri.out)",
    )
    p.add_argument(
        "--coords",
        required=True,
        metavar="FILE",
        help="JSON file with pre-computed inversion coordinates (from find_chr8_inversion.py)",
    )
    p.add_argument("--chrom", default="chr8", help="Chromosome name (default: chr8)")
    p.add_argument(
        "--out",
        default="chr8_figure",
        help="Output basename (.pdf and .png appended, default: chr8_figure)",
    )
    p.add_argument(
        "--dpi", type=int, default=300, help="DPI for PNG output (default: 300)"
    )
    p.add_argument(
        "--cfg", dest="plotsr_cfg", metavar="FILE", help="Path to plotsr config file"
    )
    return p.parse_args()


# ---------------------------------------------------------------------------
# Write plotsr input files
# ---------------------------------------------------------------------------


def _ft_tag(filepath: str) -> str:
    """Return plotsr file-type tag based on extension."""
    return "cl" if filepath.endswith(".cl") else "fa"


def write_genomes(path: str, ref: str, mat: str, pat: str) -> None:
    with open(path, "w") as fh:
        fh.write("#file\tname\ttags\n")
        fh.write(f"{ref}\tREF\tft:{_ft_tag(ref)};lw:1.5;lc:#1f77b4\n")
        fh.write(f"{mat}\tMAT\tft:{_ft_tag(mat)};lw:1.5;lc:#2ca02c\n")
        fh.write(f"{pat}\tPAT\tft:{_ft_tag(pat)};lw:1.5;lc:#d62728\n")


def write_markers(
    path: str,
    chrom: str,
    inv_start: int,
    inv_end: int,
    genome_id: str = "PAT",
    excl_start: int | None = None,
    excl_end: int | None = None,
) -> None:
    """Mark inversion breakpoints and optional excluded region boundaries."""
    with open(path, "w") as fh:
        fh.write("#chr\tstart\tend\tgenome_id\ttags\n")
        # Inversion breakpoints on PAT (triangles, hidden text)
        _hidden = "tt: ;tp:0.01;ts:1;tf:Arial;tc:white"
        fh.write(
            f"{chrom}\t{inv_start}\t{inv_start + 1}\t{genome_id}\t"
            f"mt:v;mc:#d62728;ms:4;{_hidden}\n"
        )
        fh.write(
            f"{chrom}\t{inv_end}\t{inv_end + 1}\t{genome_id}\t"
            f"mt:v;mc:#d62728;ms:4;{_hidden}\n"
        )
        # Excluded region boundaries on REF (diamonds)
        if excl_start is not None and excl_end is not None:
            fh.write(
                f"{chrom}\t{excl_start}\t{excl_start + 1}\tREF\t"
                f"mt:D;mc:#555555;ms:3;{_hidden}\n"
            )
            fh.write(
                f"{chrom}\t{excl_end}\t{excl_end + 1}\tREF\t"
                f"mt:D;mc:#555555;ms:3;{_hidden}\n"
            )


# ---------------------------------------------------------------------------
# Run plotsr
# ---------------------------------------------------------------------------


def _post_process_axes(fig, gridline_span: float = 0.8) -> None:
    """Tidy plotsr-rendered axes:
    - Drop the hardcoded "Reference Chromosome ID" y-axis label
    - Replace full-height vertical dashed gridlines with shorter ones
      centered on the chromosome tracks (gridline_span is the fraction of
      the data area height that the dashed line spans).
    """
    ymin = 0.5 - gridline_span / 2
    ymax = 0.5 + gridline_span / 2
    for ax in fig.axes:
        if ax.get_ylabel() == "Reference Chromosome ID":
            ax.set_ylabel("")
        # Capture tick positions, disable the full-axis grid, then redraw
        # each gridline as a short vertical segment.
        xticks = [
            t for t in ax.get_xticks() if ax.get_xlim()[0] <= t <= ax.get_xlim()[1]
        ]
        if xticks:
            ax.xaxis.grid(False)
            for x in xticks:
                ax.axvline(
                    x,
                    ymin=ymin,
                    ymax=ymax,
                    color="lightgray",
                    linestyle="--",
                    linewidth=0.5,
                    zorder=0,
                )


def run_plotsr(
    genomes: str,
    sr_files: list,
    intermediate_pdf: str,
    markers: str = None,
    region: str = None,
    width: float = 7,
    height: float = 2,
    fontsize: int = 8,
    cfg: str | None = None,
) -> plt.Figure:
    """Render a plotsr panel via the Python API.

    Returns the matplotlib Figure (axes post-processed, figure NOT closed).
    The caller owns the figure; call _fig_to_array() then plt.close() on it.

    Also saves a PDF of the panel to intermediate_pdf for inspection.
    Compositing does NOT read from that file — it uses the in-memory figure.
    """
    import argparse as _argparse
    from plotsr.scripts.plotsr import plotsr as _plotsr_run

    captured: dict = {}
    _orig_savefig = plt.Figure.savefig

    def _intercept(self, *args, **kwargs):  # noqa: ANN001
        captured["fig"] = self

    plt.Figure.savefig = _intercept
    sr_handles: list = []
    genomes_h = None
    markers_h = None
    cfg_h = None
    log_h = None
    log_path = Path(intermediate_pdf).with_suffix(".plotsr.log")
    try:
        sr_handles = [open(p) for p in sr_files]
        genomes_h = open(genomes)
        markers_h = open(markers) if markers else None
        cfg_h = open(cfg) if cfg else None
        log_h = open(log_path, "w")
        ns = _argparse.Namespace(
            sr=sr_handles,
            bp=None,
            genomes=genomes_h,
            markers=markers_h,
            tracks=None,
            chrord=None,
            chrname=None,
            o=intermediate_pdf,
            itx=False,
            chr=None,
            reg=region,
            rtr=False,
            nosyn=False,
            noinv=False,
            notr=False,
            nodup=False,
            s=10000,
            cfg=cfg_h,
            R=False,
            f=fontsize,
            H=height,
            W=width,
            S=0.5,
            d=300,
            b="agg",
            v=False,
            log="WARN",
            logfin=log_h,
        )
        print(f"  plotsr(reg={region}, W={width}, H={height}, f={fontsize})")
        _plotsr_run(ns)
    finally:
        plt.Figure.savefig = _orig_savefig
        for h in sr_handles:
            h.close()
        if genomes_h:
            genomes_h.close()
        if markers_h:
            markers_h.close()
        if cfg_h:
            cfg_h.close()
        if log_h:
            log_h.close()

    fig = captured.get("fig")
    if fig is None:
        raise RuntimeError("plotsr did not produce a figure")

    _post_process_axes(fig, gridline_span=GRIDLINE_SPAN)

    # Save intermediate PDF for inspection — compositing uses the in-memory fig.
    fig.savefig(intermediate_pdf, bbox_inches="tight", pad_inches=0.01)
    print(f"  -> wrote {intermediate_pdf}")

    return fig  # caller closes after _fig_to_array()


def _fig_to_array(fig: plt.Figure, dpi: int, crop_frac: float) -> np.ndarray:
    """Rasterize a matplotlib figure to a numpy array and close it.

    Renders at `dpi` using the agg backend (no file I/O), then crops the top
    `crop_frac` fraction to remove the plotsr legend.
    """
    fig.set_dpi(dpi)
    fig.canvas.draw()
    arr = np.asarray(fig.canvas.buffer_rgba())[:, :, :3]  # drop alpha
    plt.close(fig)
    crop_px = int(arr.shape[0] * crop_frac)
    return arr[crop_px:, :, :]


# ---------------------------------------------------------------------------
# Assemble multi-panel figure
# ---------------------------------------------------------------------------


def _panel_label(ax, letter: str) -> None:
    ax.text(
        -0.04,
        1.02,
        letter,
        transform=ax.transAxes,
        fontsize=PANEL_LABEL_FONTSIZE,
        fontweight="bold",
        va="bottom",
        ha="left",
    )


def assemble_figure(
    full_arr: np.ndarray,
    zoom_arr: np.ndarray,
    out_base: str,
    inv_start: int,
    inv_end: int,
    chrom: str,
    ref_start: int,
    ref_end: int,
    excl_start: int | None = None,
    excl_end: int | None = None,
    dpi: int = 300,
) -> None:
    fig = plt.figure(figsize=(FIG_WIDTH_INCHES, FIG_HEIGHT_INCHES))
    gs = GridSpec(
        2,
        1,
        figure=fig,
        hspace=0.0,
        left=0.0,
        right=1.0,
        top=1.0,
        bottom=0.00,
        height_ratios=[1.0, 1.0],
    )

    # -- Panel A: full synteny ------------------------------------------------
    ax_a = fig.add_subplot(gs[0])
    ax_a.imshow(full_arr, aspect="equal", interpolation="none")
    ax_a.axis("off")
    _panel_label(ax_a, "A")

    # -- Panel B: zoomed inversion --------------------------------------------
    ax_b = fig.add_subplot(gs[1])
    ax_b.imshow(zoom_arr, aspect="equal", interpolation="none")
    ax_b.axis("on")
    _panel_label(ax_b, "B")

    #    ax_b.set_title(
    #        f"Zoomed: REF {ref_start / 1e6:.1f} – {ref_end / 1e6:.1f} Mb"
    #        f"{excl_label}",
    #        fontsize=ZOOM_TITLE_FONTSIZE,
    #        pad=4,
    #    )

    # -- Save -----------------------------------------------------------------
    Path(out_base).parent.mkdir(parents=True, exist_ok=True)
    for ext, kwargs in [("pdf", {}), ("png", {"dpi": dpi})]:
        path = f"{out_base}.{ext}"
        fig.savefig(path, bbox_inches="tight", **kwargs)
        print(f"  Saved: {path}")

    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    args = parse_args()

    # Validate inputs
    missing = [
        f
        for f in [args.ref, args.mat, args.pat, args.syri_rm, args.syri_mp, args.coords]
        if not Path(f).exists()
    ]
    if missing:
        print("[ERROR] Missing input files:")
        for m in missing:
            print(f"  {m}")
        sys.exit(1)

    with open(args.coords) as fh:
        coords = json.load(fh)

    inv_start = coords["qry_start"]
    inv_end = coords["qry_end"]
    ref_inv_start = coords["ref_start"]
    ref_inv_end = coords["ref_end"]
    z_start = coords["zoom_ref_start"]
    z_end = coords["zoom_ref_end"]

    # Excluded region from PAV callset (hardcoded)
    excl_start = 8_237_843
    excl_end = 12_234_345

    print(
        "[INFO] Inversion coordinates from coords file: "
        f"REF {args.chrom}:{ref_inv_start}-{ref_inv_end}, "
        f"PAT {args.chrom}:{inv_start}-{inv_end}, "
        f"size {coords['size']:,} bp"
    )
    print(f"[INFO] Excluded region (PAV): REF {args.chrom}:{excl_start}-{excl_end}")

    sr_files = [args.syri_rm, args.syri_mp]

    # Determine working directory for intermediate plotsr files
    out_dir = Path(args.out).parent
    out_dir.mkdir(parents=True, exist_ok=True)

    genomes_path = str(out_dir / "genomes.txt")
    markers_path = str(out_dir / "markers.bed")
    full_pdf = str(out_dir / "plotsr_chr8_full.pdf")
    zoom_pdf = str(out_dir / "plotsr_chr8_zoom.pdf")

    # 1. Write plotsr support files
    print("\n[1/4] Writing plotsr input files ...")
    write_genomes(genomes_path, args.ref, args.mat, args.pat)
    write_markers(
        markers_path,
        args.chrom,
        inv_start,
        inv_end,
        excl_start=excl_start,
        excl_end=excl_end,
    )

    # 2. Panel A: full chr8 synteny (compact height), capped at 100 Mb
    print(f"\n[2/4] Generating full {args.chrom} synteny panel ...")
    fig_full = run_plotsr(
        genomes=genomes_path,
        sr_files=sr_files,
        intermediate_pdf=full_pdf,
        markers=markers_path,
        region=f"REF:{args.chrom}:1-100000000",
        width=FIG_WIDTH_INCHES,
        height=PANEL_A_HEIGHT,
        fontsize=PLOTSR_FONTSIZE,
        cfg=args.plotsr_cfg,
    )

    # 3. Panel B: zoomed inversion
    print(f"\n[3/4] Generating zoomed inversion panel ({z_start} – {z_end} bp) ...")
    fig_zoom = run_plotsr(
        genomes=genomes_path,
        sr_files=sr_files,
        intermediate_pdf=zoom_pdf,
        markers=markers_path,
        region=f"REF:{args.chrom}:{z_start}-{z_end}",
        width=FIG_WIDTH_INCHES,
        height=PANEL_B_HEIGHT,
        fontsize=PLOTSR_FONTSIZE,
        cfg=args.plotsr_cfg,
    )

    # Rasterize panels from memory at INTERMEDIATE_DPI — no file roundtrip.
    print(f"\n  Rasterizing panels at {INTERMEDIATE_DPI} DPI ...")
    full_arr = _fig_to_array(fig_full, dpi=INTERMEDIATE_DPI, crop_frac=PANEL_A_CROP)
    zoom_arr = _fig_to_array(fig_zoom, dpi=INTERMEDIATE_DPI, crop_frac=PANEL_B_CROP)

    # 4. Assemble figure
    print(f"\n[4/4] Assembling multi-panel figure -> {args.out}.pdf / .png ...")
    assemble_figure(
        full_arr=full_arr,
        zoom_arr=zoom_arr,
        out_base=args.out,
        inv_start=inv_start,
        inv_end=inv_end,
        chrom=args.chrom,
        ref_start=ref_inv_start,
        ref_end=ref_inv_end,
        excl_start=excl_start,
        excl_end=excl_end,
        dpi=args.dpi,
    )
    print("\nDone.")


if __name__ == "__main__":
    main()
