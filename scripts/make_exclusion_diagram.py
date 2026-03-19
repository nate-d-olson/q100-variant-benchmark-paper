#!/usr/bin/env python3
"""Generate a schematic diagram showing how exclusions are applied to
produce benchmark regions from dipcall assembly regions (dip.bed).

Illustrates:
  - Three exclusion categories (assembly-agnostic, assembly-intersect, ref-agnostic)
  - BED operations: slop, slopmerge, assembly-break filtering
  - Final subtraction to produce benchmark.bed

IGV Snapshot Alternative
------------------------
If a real-data IGV screenshot is preferred over this schematic, the following
GRCh38 regions have multiple overlapping exclusion types and are compact enough
for clear visualization. Load BED files from
resources/exclusions/v5.0q_GRCh38_smvar/ as separate colored tracks in IGV.

  chr22:12,423,690-12,503,690  (80 kb, 6 types: flanks, gaps, satellites,
                                segdups, svs-and-simple-repeats, tandem-repeats)
  chr15:19,710,254-19,816,277  (106 kb, 6 types — same set, near centromere)
  chr5:47,001,864-47,168,439   (166 kb, 6 types)
  chr16:37,914,782-37,971,525  (56 kb, 5 types — very compact)
  chr4:49,068,543-49,127,006   (58 kb, 5 types — compact with segdups)

Suggested IGV session setup:
  1. Load GRCh38 reference
  2. Load dip.bed (light blue) — shows assembly coverage
  3. Load each exclusion BED as a separate track with distinct colors:
     - segdups (purple), tandem-repeats (purple, lighter), satellites (purple, lightest)
     - gaps (orange), VDJ (orange, lighter), HG002Q100-errors (orange, lighter)
     - flanks (green), svs-and-simple-repeats (green, lighter)
  4. Load final benchmark BED (dark blue) — shows remaining regions after subtraction
  5. Navigate to one of the regions above
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import numpy as np
from pathlib import Path


# ---------- colour palette (colorblind-friendly) ----------
C = {
    "dip": "#B8C9E1",           # light blue — dip.bed
    "bench": "#4A90D9",         # blue — benchmark regions
    "asm_agnostic": "#E07B54",  # orange — asm-agnostic exclusions
    "asm_intersect": "#8E6BB0", # purple — asm-intersect exclusions
    "ref_agnostic": "#5AAA46",  # green — ref-agnostic exclusions
    "slop": "#FADCB0",          # light orange — slop buffer
    "break": "#D94E4E",         # red — assembly breaks
    "excluded": "#F0F0F0",      # light grey — excluded from benchmark
    "bg": "#FFFFFF",
}


def draw_bar(ax, start, end, y, height, color, edgecolor="grey",
             linewidth=0.5, alpha=1.0, label=None, zorder=2):
    """Draw a horizontal bar representing a genomic interval."""
    rect = mpatches.FancyBboxPatch(
        (start, y - height / 2), end - start, height,
        boxstyle="round,pad=0.002", facecolor=color,
        edgecolor=edgecolor, linewidth=linewidth, alpha=alpha, zorder=zorder,
    )
    ax.add_patch(rect)
    return rect


def add_bracket(ax, x1, x2, y, text, color="grey", fontsize=7):
    """Add a bracket annotation with label."""
    ax.annotate(
        "", xy=(x1, y), xytext=(x2, y),
        arrowprops=dict(arrowstyle="<->", color=color, lw=1),
    )
    ax.text((x1 + x2) / 2, y + 0.12, text, ha="center", va="bottom",
            fontsize=fontsize, color=color)


def panel_a(ax):
    """Panel A: BED operations — slop, slopmerge, assembly-break filtering."""

    ax.set_xlim(-0.5, 24)
    ax.set_ylim(-1, 8.5)
    ax.set_title("A. BED Operations for Exclusion Processing",
                 fontsize=10, fontweight="bold", loc="left", pad=8)

    # --- Row 1: Slop ---
    y = 7.0
    ax.text(-0.3, y, "slop", fontsize=9, fontweight="bold", va="center")
    # Original interval
    draw_bar(ax, 3, 7, y + 0.5, 0.4, C["asm_agnostic"], label="Original")
    ax.text(5, y + 0.5, "Original", ha="center", va="center", fontsize=6.5, color="white", fontweight="bold")
    # After slop
    draw_bar(ax, 0.5, 9.5, y - 0.5, 0.4, C["slop"])
    draw_bar(ax, 3, 7, y - 0.5, 0.4, C["asm_agnostic"])
    ax.text(5, y - 0.5, "+15 kb each side", ha="center", va="center", fontsize=6.5)
    add_bracket(ax, 0.5, 3, y - 0.1, "15 kb", color="#888")
    add_bracket(ax, 7, 9.5, y - 0.1, "15 kb", color="#888")

    # --- Row 2: Slopmerge ---
    y = 4.5
    ax.text(-0.3, y, "slopmerge", fontsize=9, fontweight="bold", va="center")
    # Two nearby intervals
    draw_bar(ax, 2, 5, y + 0.8, 0.4, C["asm_intersect"])
    draw_bar(ax, 7, 10, y + 0.8, 0.4, C["asm_intersect"])
    ax.text(3.5, y + 0.8, "A", ha="center", va="center", fontsize=7, color="white", fontweight="bold")
    ax.text(8.5, y + 0.8, "B", ha="center", va="center", fontsize=7, color="white", fontweight="bold")
    add_bracket(ax, 5, 7, y + 1.2, "gap < 10 kb", color="#888")

    # After slop (before merge)
    draw_bar(ax, 0, 12, y, 0.35, C["slop"], alpha=0.5)
    draw_bar(ax, 0, 7.5, y + 0.15, 0.15, C["asm_intersect"], alpha=0.3)
    draw_bar(ax, 4.5, 12, y - 0.15, 0.15, C["asm_intersect"], alpha=0.3)
    ax.text(6, y, "+15 kb slop → overlap", ha="center", va="center", fontsize=6.5)

    # After merge
    draw_bar(ax, 0, 12, y - 1.0, 0.4, C["asm_intersect"])
    ax.text(6, y - 1.0, "Merged single region", ha="center", va="center", fontsize=6.5, color="white", fontweight="bold")

    # --- Row 3: Assembly-break filtering ---
    y = 1.5
    ax.text(-0.3, y, "asm-intersect", fontsize=9, fontweight="bold", va="center")

    # dip.bed regions (baseline) with breaks
    draw_bar(ax, 1, 8, y + 0.8, 0.3, C["dip"])
    draw_bar(ax, 9, 14, y + 0.8, 0.3, C["dip"])
    # Break markers
    ax.plot([8, 8], [y + 0.5, y + 1.1], color=C["break"], lw=2, zorder=3)
    ax.plot([9, 9], [y + 0.5, y + 1.1], color=C["break"], lw=2, zorder=3)
    ax.text(8.5, y + 1.25, "break", ha="center", va="bottom", fontsize=6,
            color=C["break"], fontweight="bold")
    ax.text(4.5, y + 0.8, "dip.bed", ha="center", va="center", fontsize=6.5)

    # Exclusion regions (e.g. segdups)
    draw_bar(ax, 3, 6, y, 0.3, C["asm_intersect"], alpha=0.4)
    draw_bar(ax, 7.5, 11, y, 0.3, C["asm_intersect"], alpha=0.4)
    draw_bar(ax, 13, 15, y, 0.3, C["asm_intersect"], alpha=0.4)
    ax.text(4.5, y, "A", ha="center", va="center", fontsize=6.5)
    ax.text(9.25, y, "B ✓", ha="center", va="center", fontsize=6.5, fontweight="bold")
    ax.text(14, y, "C", ha="center", va="center", fontsize=6.5)

    # Result — only region B overlaps break
    draw_bar(ax, 7.5, 11, y - 1.0, 0.3, C["asm_intersect"])
    ax.text(9.25, y - 1.0, "Only B excluded", ha="center", va="center",
            fontsize=6.5, color="white", fontweight="bold")
    ax.text(4.5, y - 1.0, "A: no break overlap → kept", ha="center", va="center",
            fontsize=6, color="#888", style="italic")

    ax.set_axis_off()


def panel_b(ax):
    """Panel B: Exclusion categories and subtraction flow."""

    ax.set_xlim(-1, 26)
    ax.set_ylim(-1.5, 9)
    ax.set_title("B. Exclusion Categories and Benchmark Region Construction",
                 fontsize=10, fontweight="bold", loc="left", pad=8)

    # --- dip.bed row ---
    y = 8
    ax.text(-0.5, y, "dip.bed", fontsize=8, fontweight="bold", va="center")
    draw_bar(ax, 2, 25, y, 0.5, C["dip"])
    ax.text(13.5, y, "Dipcall assembly regions", ha="center", va="center",
            fontsize=7, color="#333")

    # --- Assembly-agnostic exclusions ---
    y = 6.2
    ax.text(-0.5, y, "Assembly-\nagnostic", fontsize=7.5, fontweight="bold",
            va="center", color=C["asm_agnostic"], linespacing=1.3)
    # gaps (with slop)
    draw_bar(ax, 3, 3.5, y + 0.3, 0.3, C["slop"])
    draw_bar(ax, 3.1, 3.4, y + 0.3, 0.3, C["asm_agnostic"])
    ax.text(3.25, y + 0.65, "gaps\n(+slop)", ha="center", va="bottom", fontsize=5.5)
    # VDJ
    draw_bar(ax, 8, 10, y + 0.3, 0.3, C["asm_agnostic"])
    ax.text(9, y + 0.65, "VDJ", ha="center", va="bottom", fontsize=5.5)
    # HG002Q100-errors
    draw_bar(ax, 15, 16.5, y + 0.3, 0.3, C["asm_agnostic"])
    ax.text(15.75, y + 0.65, "errors", ha="center", va="bottom", fontsize=5.5)
    # PAV-inversions
    draw_bar(ax, 20, 23, y + 0.3, 0.3, C["asm_agnostic"])
    ax.text(21.5, y + 0.65, "PAV inv.", ha="center", va="bottom", fontsize=5.5)
    # Description
    ax.text(13.5, y - 0.3, "Entire regions excluded regardless of assembly quality",
            ha="center", va="top", fontsize=6, color="#666", style="italic")

    # --- Assembly-intersect exclusions ---
    y = 4.2
    ax.text(-0.5, y, "Assembly-\nintersect", fontsize=7.5, fontweight="bold",
            va="center", color=C["asm_intersect"], linespacing=1.3)
    # segdups (with slopmerge, filtered to breaks)
    draw_bar(ax, 5, 7, y + 0.3, 0.3, C["asm_intersect"])
    ax.text(6, y + 0.65, "segdups\n(+slopmerge, breaks only)", ha="center",
            va="bottom", fontsize=5.5)
    # tandem-repeats (with slop, filtered to breaks)
    draw_bar(ax, 12, 13.5, y + 0.3, 0.3, C["asm_intersect"])
    ax.text(12.75, y + 0.65, "TR\n(+slop, breaks only)", ha="center",
            va="bottom", fontsize=5.5)
    # satellites (with slopmerge, filtered)
    draw_bar(ax, 18, 20, y + 0.3, 0.3, C["asm_intersect"])
    ax.text(19, y + 0.65, "satellites\n(+slopmerge, breaks only)", ha="center",
            va="bottom", fontsize=5.5)
    ax.text(13.5, y - 0.3, "Only excluded where assembly has alignment breaks",
            ha="center", va="top", fontsize=6, color="#666", style="italic")

    # --- Ref-agnostic exclusions ---
    y = 2.2
    ax.text(-0.5, y, "Ref-\nagnostic", fontsize=7.5, fontweight="bold",
            va="center", color=C["ref_agnostic"], linespacing=1.3)
    # flanks
    draw_bar(ax, 2, 2.5, y + 0.3, 0.3, C["ref_agnostic"])
    draw_bar(ax, 25, 25.5, y + 0.3, 0.3, C["ref_agnostic"], edgecolor="grey")
    ax.text(2.25, y + 0.65, "flanks", ha="center", va="bottom", fontsize=5.5)
    # svs-and-simple-repeats
    draw_bar(ax, 7, 7.8, y + 0.3, 0.3, C["ref_agnostic"])
    draw_bar(ax, 10.5, 11, y + 0.3, 0.3, C["ref_agnostic"])
    ax.text(9, y + 0.65, "SVs ∩ repeats", ha="center", va="bottom", fontsize=5.5)
    # consecutive-svs
    draw_bar(ax, 14, 14.5, y + 0.3, 0.3, C["ref_agnostic"])
    ax.text(14.25, y + 0.65, "consec.\nSVs", ha="center", va="bottom", fontsize=5.5)
    # self-discrep
    draw_bar(ax, 17, 17.3, y + 0.3, 0.3, C["ref_agnostic"])
    ax.text(17.15, y + 0.65, "self-\ndiscrep", ha="center", va="bottom", fontsize=5.5)
    ax.text(13.5, y - 0.3, "Derived from assembly alignments and variant calls (reference-independent)",
            ha="center", va="top", fontsize=6, color="#666", style="italic")

    # --- Subtraction arrow ---
    y_arrow = 0.8
    ax.annotate("", xy=(13.5, 0.55), xytext=(13.5, 1.2),
                arrowprops=dict(arrowstyle="->", color="#333", lw=1.5))
    ax.text(14.5, 0.9, "subtract all\nexclusions", ha="left", va="center",
            fontsize=6.5, color="#555", style="italic")

    # --- Benchmark regions ---
    y = 0
    ax.text(-0.5, y, "benchmark\n.bed", fontsize=8, fontweight="bold", va="center",
            color=C["bench"])
    # Draw full dip.bed background (excluded = grey)
    draw_bar(ax, 2, 25, y, 0.5, C["excluded"], edgecolor="#ccc", linewidth=0.3)
    # Overlay benchmark (non-excluded) segments — gaps where exclusions were
    bench_segs = [(2.5, 3), (3.5, 5), (7.8, 8), (11, 12), (13.5, 14),
                  (14.5, 15), (16.5, 17), (17.3, 18), (23, 25)]
    for s, e in bench_segs:
        draw_bar(ax, s, e, y, 0.5, C["bench"])

    ax.text(13.5, y - 0.6, "dip.bed − (assembly-agnostic ∪ assembly-intersect ∪ ref-agnostic) = benchmark.bed",
            ha="center", va="top", fontsize=7, fontweight="bold", color="#333")

    ax.set_axis_off()


def main():
    fig, axes = plt.subplots(2, 1, figsize=(12, 10),
                              gridspec_kw={"height_ratios": [1, 1], "hspace": 0.3})

    panel_a(axes[0])
    panel_b(axes[1])

    # Legend
    legend_elements = [
        mpatches.Patch(facecolor=C["dip"], edgecolor="grey", label="Dipcall assembly regions (dip.bed)"),
        mpatches.Patch(facecolor=C["bench"], edgecolor="grey", label="Final benchmark regions"),
        mpatches.Patch(facecolor=C["asm_agnostic"], edgecolor="grey", label="Assembly-agnostic exclusions (gaps, VDJ, errors, PAV inv.)"),
        mpatches.Patch(facecolor=C["asm_intersect"], edgecolor="grey", label="Assembly-intersect exclusions (segdups, TR, satellites)"),
        mpatches.Patch(facecolor=C["ref_agnostic"], edgecolor="grey", label="Ref-agnostic exclusions (flanks, SVs∩repeats, self-discrep)"),
        mpatches.Patch(facecolor=C["slop"], edgecolor="grey", label="Slop buffer (±15 kb)"),
        mpatches.Patch(facecolor=C["break"], edgecolor="grey", label="Assembly alignment break"),
    ]
    fig.legend(handles=legend_elements, loc="lower center", ncol=3,
               fontsize=7.5, frameon=True, edgecolor="#ccc",
               bbox_to_anchor=(0.5, -0.02))

    outpath = Path(__file__).resolve().parent.parent / "manuscript" / "figs" / "exclusion_diagram.png"
    fig.savefig(outpath, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"Saved to {outpath}")
    plt.close()


if __name__ == "__main__":
    main()
