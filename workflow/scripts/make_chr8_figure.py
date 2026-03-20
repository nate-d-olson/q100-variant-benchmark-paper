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

The colour legend for structural annotations (syntenic, inverted, translocated,
duplicated) is embedded in the plotsr output in Panel A.

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
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt
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


def _crop_top(img: Image.Image, frac: float = 0.22) -> Image.Image:
    """Crop the top portion of a plotsr image."""
    crop_top = int(img.height * frac)
    return img.crop((0, crop_top, img.width, img.height))


def assemble_figure(
    full_png: str,
    zoom_png: str,
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

    # Crop whitespace: remove gap between legend and tracks (full panel)
    # and remove legend entirely from zoom panel
    full_img = _crop_top(Image.open(full_png), frac=0.07)
    zoom_img = _crop_top(Image.open(zoom_png), frac=0.28)

    fig = plt.figure(figsize=(10, 8))
    gs = GridSpec(
        2,
        1,
        figure=fig,
        hspace=0.08,
        left=0.02,
        right=0.98,
        top=0.95,
        bottom=0.02,
        height_ratios=[1.0, 0.8],
    )

    # -- Panel A: full synteny ------------------------------------------------
    ax_a = fig.add_subplot(gs[0])
    ax_a.imshow(full_img, aspect="equal", interpolation="bilinear")
    ax_a.axis("off")
    _panel_label(ax_a, "A")

    # -- Panel B: zoomed inversion --------------------------------------------
    ax_b = fig.add_subplot(gs[1])
    ax_b.imshow(zoom_img, aspect="equal", interpolation="bilinear")
    ax_b.axis("off")
    _panel_label(ax_b, "B")

    excl_label = ""
    if excl_start is not None and excl_end is not None:
        excl_label = (
            f"  |  Excluded: {excl_start / 1e6:.1f} – {excl_end / 1e6:.1f} Mb"
        )
    ax_b.set_title(
        f"Zoomed: REF {ref_start / 1e6:.1f} – {ref_end / 1e6:.1f} Mb"
        f"{excl_label}",
        fontsize=9,
        pad=4,
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
    print(
        f"[INFO] Excluded region (PAV): "
        f"REF {args.chrom}:{excl_start}-{excl_end}"
    )

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
    write_markers(
        markers_path,
        args.chrom,
        inv_start,
        inv_end,
        excl_start=excl_start,
        excl_end=excl_end,
    )

    # 2. Panel A: full chr8 synteny (compact height)
    print(f"\n[2/4] Generating full {args.chrom} synteny panel ...")
    run_plotsr(
        genomes=genomes_path,
        sr_files=sr_files,
        output=full_png,
        markers=markers_path,
        region=None,
        width=10,
        height=3,
        fontsize=8,
    )

    # 3. Panel B: zoomed inversion
    print(f"\n[3/4] Generating zoomed inversion panel ({z_start} – {z_end} bp) ...")
    run_plotsr(
        genomes=genomes_path,
        sr_files=sr_files,
        output=zoom_png,
        markers=markers_path,
        region=f"REF:{args.chrom}:{z_start}-{z_end}",
        width=10,
        height=4,
        fontsize=8,
    )

    # 4. Assemble figure
    print(f"\n[4/4] Assembling multi-panel figure -> {args.out}.pdf / .png ...")
    assemble_figure(
        full_png=full_png,
        zoom_png=zoom_png,
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
