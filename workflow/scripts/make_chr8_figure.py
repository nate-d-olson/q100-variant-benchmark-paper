#!/usr/bin/env python3
"""
make_chr8_figure.py — multi-panel Chr8 synteny figure using plotsr
=============================================================

Panel A  Full Chr8: REF <-> MAT <-> PAT synteny (plotsr)
Panel B  Zoomed view of the large inversion on the PAT haplotype
Panel C  Structural-annotation colour legend

Dependencies (install via conda/bioconda):
  conda install -c bioconda plotsr syri minimap2
  conda install matplotlib Pillow

Standalone usage
----------------
  python make_chr8_figure.py \\
      --ref  ref_chr8.cl  --mat  mat_chr8.cl  --pat  pat_chr8.cl \\
      --rm   syri/ref_matsyri.out  --mp  syri/mat_patsyri.out \\
      --chrom chr8 \\
      --inv-start 50000000  --inv-end 58000000 \\
      --out  results/chr8_figure

Snakemake usage
---------------
  Called by workflow/rules/chr8_synteny.smk via the make_figure rule.
  Configure paths and inversion coordinates in config/config.yaml under
  the `chr8_synteny` key.
"""

import argparse
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from PIL import Image


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--ref", required=True, help="REF chromosome-length (.cl) or FASTA file")
    p.add_argument("--mat", required=True, help="MAT chromosome-length (.cl) or FASTA file")
    p.add_argument("--pat", required=True, help="PAT chromosome-length (.cl) or FASTA file")
    p.add_argument(
        "--rm",
        required=True,
        dest="syri_rm",
        help="SyRI output: REF vs MAT  (ref_matsyri.out)",
    )
    p.add_argument(
        "--mp",
        required=True,
        dest="syri_mp",
        help="SyRI output: MAT vs PAT  (mat_patsyri.out)",
    )
    p.add_argument("--chrom", default="chr8", help="Chromosome name (default: chr8)")
    p.add_argument(
        "--inv-start",
        type=int,
        default=50_000_000,
        help="Inversion start on PAT haplotype, bp (default: 50000000)",
    )
    p.add_argument(
        "--inv-end",
        type=int,
        default=58_000_000,
        help="Inversion end on PAT haplotype, bp (default: 58000000)",
    )
    p.add_argument(
        "--zoom-padding",
        type=int,
        default=2_000_000,
        help="Flanking bp around inversion in zoom panel (default: 2000000)",
    )
    p.add_argument(
        "--out",
        default="chr8_figure",
        help="Output basename (.pdf and .png appended, default: chr8_figure)",
    )
    p.add_argument("--dpi", type=int, default=300, help="DPI for PNG output (default: 300)")
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
    path: str, chrom: str, inv_start: int, inv_end: int, genome_id: str = "PAT"
) -> None:
    """Mark inversion breakpoints as downward triangles on the PAT haplotype."""
    with open(path, "w") as fh:
        fh.write("#chr\tstart\tend\tgenome_id\ttags\n")
        fh.write(
            f"{chrom}\t{inv_start}\t{inv_start + 1}\t{genome_id}\t"
            "mt:v;mc:#d62728;ms:4;tt:INV_start;tp:0.04;ts:7;tf:Arial;tc:#d62728\n"
        )
        fh.write(
            f"{chrom}\t{inv_end}\t{inv_end + 1}\t{genome_id}\t"
            "mt:v;mc:#d62728;ms:4;tt:INV_end;tp:0.04;ts:7;tf:Arial;tc:#d62728\n"
        )


# ---------------------------------------------------------------------------
# Run plotsr
# ---------------------------------------------------------------------------


def run_plotsr(
    genomes: str,
    sr_files: list,
    output: str,
    markers: str = None,
    region: str = None,
    width: float = 8,
    height: float = 6,
    fontsize: int = 8,
) -> None:
    cmd = ["plotsr"]
    for sr in sr_files:
        cmd += ["--sr", sr]
    cmd += [
        "--genomes",
        genomes,
        "-o",
        output,
        "-W",
        str(width),
        "-H",
        str(height),
        "-f",
        str(fontsize),
        "-S",
        "0.5",
    ]
    if markers:
        cmd += ["--markers", markers]
    if region:
        cmd += ["--reg", region]

    print("  $", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print("[STDOUT]", result.stdout)
        print("[STDERR]", result.stderr)
        raise RuntimeError(f"plotsr exited with code {result.returncode}")
    print(f"  -> wrote {output}")


# ---------------------------------------------------------------------------
# Assemble multi-panel figure
# ---------------------------------------------------------------------------


def _panel_label(ax, letter: str) -> None:
    ax.text(
        -0.04,
        1.02,
        letter,
        transform=ax.transAxes,
        fontsize=15,
        fontweight="bold",
        va="bottom",
        ha="left",
    )


def assemble_figure(
    full_png: str,
    zoom_png: str,
    out_base: str,
    inv_start: int,
    inv_end: int,
    chrom: str,
    dpi: int = 300,
) -> None:

    fig = plt.figure(figsize=(15, 11))
    gs = GridSpec(
        2,
        2,
        figure=fig,
        hspace=0.30,
        wspace=0.06,
        left=0.03,
        right=0.97,
        top=0.92,
        bottom=0.03,
        height_ratios=[1.1, 0.9],
        width_ratios=[1.6, 1.0],
    )

    # -- Panel A: full synteny (left, spans both rows) -----------------------
    ax_a = fig.add_subplot(gs[:, 0])
    ax_a.imshow(Image.open(full_png), aspect="equal", interpolation="bilinear")
    ax_a.axis("off")
    _panel_label(ax_a, "A")
    ax_a.set_title(
        f"Chromosome 8 synteny  —  REF  ·  MAT  ·  PAT",
        fontsize=11,
        pad=7,
    )

    # -- Panel B: zoomed inversion (top-right) --------------------------------
    ax_b = fig.add_subplot(gs[0, 1])
    ax_b.imshow(Image.open(zoom_png), aspect="equal", interpolation="bilinear")
    ax_b.axis("off")
    _panel_label(ax_b, "B")
    ax_b.set_title(
        f"Zoomed: PAT large inversion ({chrom})\n"
        f"{inv_start / 1e6:.1f} – {inv_end / 1e6:.1f} Mb  "
        f"({(inv_end - inv_start) / 1e6:.1f} Mb)",
        fontsize=10,
        pad=7,
    )

    # -- Panel C: legend (bottom-right) ---------------------------------------
    ax_c = fig.add_subplot(gs[1, 1])
    ax_c.set_xlim(0, 1)
    ax_c.set_ylim(0, 1)
    ax_c.axis("off")
    _panel_label(ax_c, "C")
    ax_c.set_title("Structural annotation legend", fontsize=10, pad=7)

    # SyRI / plotsr default annotation colours
    legend_patches = [
        mpatches.Patch(fc="#91d6e6", ec="#555", label="Syntenic (SYN)"),
        mpatches.Patch(fc="#fdbe85", ec="#555", label="Inversion (INV)"),
        mpatches.Patch(fc="#c2b2d9", ec="#555", label="Translocation (TRA)"),
        mpatches.Patch(fc="#c4e5b8", ec="#555", label="Inv. translocation (INVTR)"),
        mpatches.Patch(fc="#fccde5", ec="#555", label="Duplication (DUP)"),
        mpatches.Patch(fc="#d9d9d9", ec="#555", label="Inv. duplication (INVDP)"),
    ]
    ax_c.legend(
        handles=legend_patches,
        loc="upper left",
        fontsize=8.5,
        frameon=True,
        framealpha=0.9,
        edgecolor="#aaaaaa",
        title="SyRI annotation types",
        title_fontsize=9,
    )

    # Inversion callout box
    size_mb = (inv_end - inv_start) / 1e6
    ax_c.text(
        0.05,
        0.38,
        f"[v]  PAT large inversion\n"
        f"   {chrom}:{inv_start / 1e6:.2f} – {inv_end / 1e6:.2f} Mb\n"
        f"   Size approx {size_mb:.1f} Mb",
        transform=ax_c.transAxes,
        fontsize=9,
        va="top",
        color="#d62728",
        bbox=dict(boxstyle="round,pad=0.5", fc="#fff5f5", ec="#d62728", lw=1.2, alpha=0.85),
    )

    # Genome-track colour swatches
    for i, (label, color) in enumerate([
        ("REF", "#1f77b4"),
        ("MAT", "#2ca02c"),
        ("PAT", "#d62728"),
    ]):
        y = 0.12 - i * 0.05
        ax_c.plot([0.05, 0.15], [y, y], color=color, lw=3, transform=ax_c.transAxes)
        ax_c.text(
            0.18, y, label, transform=ax_c.transAxes, fontsize=9, va="center",
            color=color, fontweight="bold",
        )

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
        for f in [args.ref, args.mat, args.pat, args.syri_rm, args.syri_mp]
        if not Path(f).exists()
    ]
    if missing:
        print("[ERROR] Missing input files:")
        for m in missing:
            print(f"  {m}")
        sys.exit(1)

    sr_files = [args.syri_rm, args.syri_mp]

    # Determine working directory for intermediate plotsr files
    out_dir = Path(args.out).parent
    out_dir.mkdir(parents=True, exist_ok=True)

    genomes_path = str(out_dir / "genomes.txt")
    markers_path = str(out_dir / "markers.bed")
    full_png = str(out_dir / "plotsr_chr8_full.png")
    zoom_png = str(out_dir / "plotsr_chr8_zoom.png")

    # 1. Write plotsr support files
    print("\n[1/4] Writing plotsr input files ...")
    write_genomes(genomes_path, args.ref, args.mat, args.pat)
    write_markers(markers_path, args.chrom, args.inv_start, args.inv_end)

    # 2. Panel A: full chr8 synteny
    print(f"\n[2/4] Generating full {args.chrom} synteny panel ...")
    run_plotsr(
        genomes=genomes_path,
        sr_files=sr_files,
        output=full_png,
        markers=markers_path,
        region=args.chrom,
        width=8,
        height=7,
        fontsize=8,
    )

    # 3. Panel B: zoomed inversion
    z_start = max(0, args.inv_start - args.zoom_padding)
    z_end = args.inv_end + args.zoom_padding
    print(f"\n[3/4] Generating zoomed inversion panel ({z_start} – {z_end} bp) ...")
    run_plotsr(
        genomes=genomes_path,
        sr_files=sr_files,
        output=zoom_png,
        markers=markers_path,
        region=f"REF:{args.chrom}:{z_start}:{z_end}",
        width=7,
        height=5,
        fontsize=8,
    )

    # 4. Assemble figure
    print(f"\n[4/4] Assembling multi-panel figure -> {args.out}.pdf / .png ...")
    assemble_figure(
        full_png=full_png,
        zoom_png=zoom_png,
        out_base=args.out,
        inv_start=args.inv_start,
        inv_end=args.inv_end,
        chrom=args.chrom,
        dpi=args.dpi,
    )
    print("\nDone.")


if __name__ == "__main__":
    main()
