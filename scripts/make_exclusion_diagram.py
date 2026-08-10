#!/usr/bin/env python3
"""Generate schematic diagrams illustrating how exclusions are applied to
produce benchmark regions from dipcall assembly regions (dip.bed).

Outputs two figures:
  - exclusion_bed_operations.png : BED operations (slop, slopmerge,
    assembly-break filtering)
  - exclusion_diagram.png        : Exclusion categories and benchmark region
    construction, with shaded vertical connectors linking each exclusion
    interval to the corresponding gap in the final benchmark track

IGV Snapshot Alternative
------------------------
If a real-data IGV screenshot is preferred over this schematic, the following
GRCh38 regions have multiple overlapping exclusion types and are compact enough
for clear visualization. Load BED files from
resources/exclusions/v5.0q_GRCh38_smvar/ as separate colored tracks in IGV.

  chr22:12,423,690-12,503,690  (80 kb, 6 types)
  chr15:19,710,254-19,816,277  (106 kb, 6 types)
  chr5:47,001,864-47,168,439   (166 kb, 6 types)
  chr16:37,914,782-37,971,525  (56 kb, 5 types)
  chr4:49,068,543-49,127,006   (58 kb, 5 types)
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path


# ---------- colour palette (colorblind-friendly) ----------
C = {
    "dip": "#B8C9E1",  # light blue — dip.bed
    "bench": "#4A90D9",  # blue — benchmark regions
    "asm_agnostic": "#E07B54",  # orange — asm-agnostic exclusions
    "asm_intersect": "#8E6BB0",  # purple — asm-intersect exclusions
    "ref_agnostic": "#5AAA46",  # green — ref-agnostic exclusions
    "slop": "#FADCB0",  # light orange — slop buffer
    "break": "#D94E4E",  # red — assembly breaks
    "excluded": "#F0F0F0",  # light grey — excluded from benchmark
    "bg": "#FFFFFF",
}


def draw_bar(
    ax,
    start,
    end,
    y,
    height,
    color,
    edgecolor="grey",
    linewidth=0.5,
    alpha=1.0,
    zorder=2,
):
    """Draw a horizontal bar representing a genomic interval."""
    rect = mpatches.FancyBboxPatch(
        (start, y - height / 2),
        end - start,
        height,
        boxstyle="round,pad=0.002",
        facecolor=color,
        edgecolor=edgecolor,
        linewidth=linewidth,
        alpha=alpha,
        zorder=zorder,
    )
    ax.add_patch(rect)
    return rect


def add_bracket(ax, x1, x2, y, text, color="grey", fontsize=10):
    """Add a bracket annotation with label."""
    ax.annotate(
        "",
        xy=(x1, y),
        xytext=(x2, y),
        arrowprops=dict(arrowstyle="<->", color=color, lw=1),
    )
    ax.text(
        (x1 + x2) / 2,
        y + 0.18,
        text,
        ha="center",
        va="bottom",
        fontsize=fontsize,
        color=color,
    )


def draw_bed_operations(ax):
    """BED operations — slop, slopmerge, assembly-break filtering."""

    ax.set_xlim(-2.5, 16.5)
    ax.set_ylim(-1, 8.5)

    label_fs = 13
    inline_fs = 10
    note_fs = 9

    # --- Row 1: Slop ---
    y = 7.0
    ax.text(-2.4, y, "slop", fontsize=label_fs, fontweight="bold", va="center")
    draw_bar(ax, 3, 7, y + 0.5, 0.45, C["asm_agnostic"])
    ax.text(
        5,
        y + 0.5,
        "Original",
        ha="center",
        va="center",
        fontsize=inline_fs,
        color="white",
        fontweight="bold",
    )
    draw_bar(ax, 0.5, 9.5, y - 0.5, 0.45, C["slop"])
    draw_bar(ax, 3, 7, y - 0.5, 0.45, C["asm_agnostic"])
    ax.text(
        5, y - 0.5, "+15 kb each side", ha="center", va="center", fontsize=inline_fs
    )
    add_bracket(ax, 0.5, 3, y - 0.05, "15 kb", color="#888", fontsize=note_fs)
    add_bracket(ax, 7, 9.5, y - 0.05, "15 kb", color="#888", fontsize=note_fs)

    # --- Row 2: Slopmerge ---
    y = 4.5
    ax.text(-2.4, y, "slopmerge", fontsize=label_fs, fontweight="bold", va="center")
    draw_bar(ax, 2, 5, y + 0.85, 0.45, C["asm_intersect"])
    draw_bar(ax, 7, 10, y + 0.85, 0.45, C["asm_intersect"])
    ax.text(
        3.5,
        y + 0.85,
        "A",
        ha="center",
        va="center",
        fontsize=inline_fs,
        color="white",
        fontweight="bold",
    )
    ax.text(
        8.5,
        y + 0.85,
        "B",
        ha="center",
        va="center",
        fontsize=inline_fs,
        color="white",
        fontweight="bold",
    )
    add_bracket(ax, 5, 7, y + 1.3, "gap < 10 kb", color="#888", fontsize=note_fs)

    draw_bar(ax, 0, 12, y, 0.4, C["slop"], alpha=0.5)
    draw_bar(ax, 0, 7.5, y + 0.18, 0.18, C["asm_intersect"], alpha=0.3)
    draw_bar(ax, 4.5, 12, y - 0.18, 0.18, C["asm_intersect"], alpha=0.3)
    ax.text(6, y, "+15 kb slop → overlap", ha="center", va="center", fontsize=inline_fs)

    draw_bar(ax, 0, 12, y - 1.05, 0.45, C["asm_intersect"])
    ax.text(
        6,
        y - 1.05,
        "Merged single region",
        ha="center",
        va="center",
        fontsize=inline_fs,
        color="white",
        fontweight="bold",
    )

    # --- Row 3: Assembly-break filtering ---
    y = 1.5
    ax.text(-2.4, y, "asm-intersect", fontsize=label_fs, fontweight="bold", va="center")

    draw_bar(ax, 1, 8, y + 0.85, 0.35, C["dip"])
    draw_bar(ax, 9, 14, y + 0.85, 0.35, C["dip"])
    ax.plot([8, 8], [y + 0.55, y + 1.15], color=C["break"], lw=2.5, zorder=3)
    ax.plot([9, 9], [y + 0.55, y + 1.15], color=C["break"], lw=2.5, zorder=3)
    ax.text(
        8.5,
        y + 1.32,
        "break",
        ha="center",
        va="bottom",
        fontsize=note_fs,
        color=C["break"],
        fontweight="bold",
    )
    ax.text(4.5, y + 0.85, "dip.bed", ha="center", va="center", fontsize=inline_fs)

    draw_bar(ax, 3, 6, y, 0.35, C["asm_intersect"], alpha=0.4)
    draw_bar(ax, 7.5, 11, y, 0.35, C["asm_intersect"], alpha=0.4)
    draw_bar(ax, 13, 15, y, 0.35, C["asm_intersect"], alpha=0.4)
    ax.text(4.5, y, "A", ha="center", va="center", fontsize=inline_fs)
    ax.text(
        9.25, y, "B ✓", ha="center", va="center", fontsize=inline_fs, fontweight="bold"
    )
    ax.text(14, y, "C", ha="center", va="center", fontsize=inline_fs)

    draw_bar(ax, 7.5, 11, y - 1.05, 0.35, C["asm_intersect"])
    ax.text(
        9.25,
        y - 1.05,
        "Only B excluded",
        ha="center",
        va="center",
        fontsize=inline_fs,
        color="white",
        fontweight="bold",
    )
    ax.text(
        4.5,
        y - 1.05,
        "A: no break overlap → kept",
        ha="center",
        va="center",
        fontsize=note_fs,
        color="#888",
        style="italic",
    )

    ax.set_axis_off()


def draw_categories(ax):
    """Exclusion categories with shaded connectors to benchmark track."""

    # x-axis spans 2-25 for genomic intervals; left margin holds track labels
    ax.set_xlim(-9, 26)
    ax.set_ylim(-1.5, 9.6)

    label_fs = 11
    interval_fs = 9
    note_fs = 9
    eq_fs = 11

    # Track y-positions
    y_dip = 8.6
    y_asm_agn = 6.5
    y_asm_int = 4.4
    y_ref_agn = 2.3
    y_bench = 0.0
    bench_top = y_bench + 0.3  # benchmark bar top edge

    # --- Exclusion intervals: (start, end) per track ---
    asm_agn_intervals = [
        (3.1, 3.4, "gaps (+slop)", C["slop"]),  # slop drawn separately
        (8, 10, "VDJ", None),
        (15, 16.5, "errors", None),
        (20, 23, "PAV inv.", None),
    ]
    asm_int_intervals = [
        (5, 7, "segdups\n(+slopmerge, breaks only)"),
        (12, 13.5, "TR\n(+slop, breaks only)"),
        (18, 20, "satellites\n(+slopmerge, breaks only)"),
    ]
    ref_agn_intervals = [
        (2, 2.5, "flanks"),
        (25, 25.5, None),  # right-edge flank, label shared
        (7, 7.8, "SVs ∩ repeats"),
        (10.5, 11, None),  # second SVs∩repeats segment
        (14, 14.5, "consec.\nSVs"),
        (17, 17.3, "self-\ndiscrep"),
    ]

    # --- Shaded vertical connectors (drawn first, behind everything else) ---
    # Each exclusion interval drops a translucent band down to the benchmark
    # track, visually linking it to the gap it carves out.
    def add_connector(s, e, y_top, color):
        ax.add_patch(
            mpatches.Rectangle(
                (s, bench_top),
                e - s,
                y_top - bench_top,
                facecolor=color,
                edgecolor="none",
                alpha=0.13,
                zorder=1,
            )
        )

    for s, e, _, _ in asm_agn_intervals:
        add_connector(s, e, y_asm_agn - 0.15, C["asm_agnostic"])
    # Include the gaps slop region for the connector (wider than core)
    add_connector(3, 3.5, y_asm_agn - 0.15, C["asm_agnostic"])
    for s, e, _ in asm_int_intervals:
        add_connector(s, e, y_asm_int - 0.15, C["asm_intersect"])
    for s, e, _ in ref_agn_intervals:
        add_connector(s, e, y_ref_agn - 0.15, C["ref_agnostic"])

    # --- dip.bed track ---
    ax.text(
        -9,
        y_dip,
        "Dipcall assembly\nregions (dip.bed)",
        fontsize=label_fs,
        fontweight="bold",
        va="center",
        ha="left",
        color="#333",
        linespacing=1.25,
    )
    draw_bar(ax, 2, 25, y_dip, 0.55, C["dip"])

    # --- Assembly-agnostic ---
    ax.text(
        -9,
        y_asm_agn,
        "Assembly-agnostic exclusions\n(gaps, VDJ, errors, PAV inv.)",
        fontsize=label_fs,
        fontweight="bold",
        va="center",
        ha="left",
        color=C["asm_agnostic"],
        linespacing=1.25,
    )
    # gaps with slop buffer
    draw_bar(ax, 3, 3.5, y_asm_agn + 0.3, 0.35, C["slop"])
    draw_bar(ax, 3.1, 3.4, y_asm_agn + 0.3, 0.35, C["asm_agnostic"])
    ax.text(
        3.25,
        y_asm_agn + 0.7,
        "gaps\n(+slop)",
        ha="center",
        va="bottom",
        fontsize=interval_fs - 1,
        linespacing=1.1,
    )
    # other intervals
    for s, e, label, _ in asm_agn_intervals[1:]:
        draw_bar(ax, s, e, y_asm_agn + 0.3, 0.35, C["asm_agnostic"])
        ax.text(
            (s + e) / 2,
            y_asm_agn + 0.7,
            label,
            ha="center",
            va="bottom",
            fontsize=interval_fs,
        )
    ax.text(
        13.5,
        y_asm_agn - 0.4,
        "Entire regions excluded regardless of assembly quality",
        ha="center",
        va="top",
        fontsize=note_fs,
        color="#666",
        style="italic",
    )

    # --- Assembly-intersect ---
    ax.text(
        -9,
        y_asm_int,
        "Assembly-intersect exclusions\n(segdups, TR, satellites)",
        fontsize=label_fs,
        fontweight="bold",
        va="center",
        ha="left",
        color=C["asm_intersect"],
        linespacing=1.25,
    )
    for s, e, label in asm_int_intervals:
        draw_bar(ax, s, e, y_asm_int + 0.3, 0.35, C["asm_intersect"])
        ax.text(
            (s + e) / 2,
            y_asm_int + 0.7,
            label,
            ha="center",
            va="bottom",
            fontsize=interval_fs - 1,
            linespacing=1.1,
        )
    ax.text(
        13.5,
        y_asm_int - 0.4,
        "Only excluded where assembly has alignment breaks",
        ha="center",
        va="top",
        fontsize=note_fs,
        color="#666",
        style="italic",
    )

    # --- Ref-agnostic ---
    ax.text(
        -9,
        y_ref_agn,
        "Ref-agnostic exclusions\n(flanks, SVs∩repeats, self-discrep)",
        fontsize=label_fs,
        fontweight="bold",
        va="center",
        ha="left",
        color=C["ref_agnostic"],
        linespacing=1.25,
    )
    for s, e, label in ref_agn_intervals:
        draw_bar(ax, s, e, y_ref_agn + 0.3, 0.35, C["ref_agnostic"])
        if label is not None:
            ax.text(
                (s + e) / 2,
                y_ref_agn + 0.7,
                label,
                ha="center",
                va="bottom",
                fontsize=interval_fs,
                linespacing=1.1,
            )
    ax.text(
        13.5,
        y_ref_agn - 0.4,
        "Derived from assembly alignments and variant calls (reference-independent)",
        ha="center",
        va="top",
        fontsize=note_fs,
        color="#666",
        style="italic",
    )

    # --- Subtraction arrow ---
    ax.annotate(
        "",
        xy=(13.5, 0.7),
        xytext=(13.5, 1.4),
        arrowprops=dict(arrowstyle="->", color="#333", lw=1.8),
    )
    ax.text(
        14.3,
        1.05,
        "subtract all\nexclusions",
        ha="left",
        va="center",
        fontsize=note_fs,
        color="#555",
        style="italic",
        linespacing=1.25,
    )

    # --- Benchmark regions ---
    ax.text(
        -9,
        y_bench,
        "Final benchmark\nregions",
        fontsize=label_fs,
        fontweight="bold",
        va="center",
        ha="left",
        color=C["bench"],
        linespacing=1.25,
    )
    draw_bar(ax, 2, 25, y_bench, 0.6, C["excluded"], edgecolor="#ccc", linewidth=0.3)
    bench_segs = [
        (2.5, 3),
        (3.5, 5),
        (7.8, 8),
        (11, 12),
        (13.5, 14),
        (14.5, 15),
        (16.5, 17),
        (17.3, 18),
        (23, 25),
    ]
    for s, e in bench_segs:
        draw_bar(ax, s, e, y_bench, 0.6, C["bench"])

    ax.text(
        13.5,
        y_bench - 0.7,
        "dip.bed − (assembly-agnostic ∪ assembly-intersect ∪ "
        "ref-agnostic) = benchmark.bed",
        ha="center",
        va="top",
        fontsize=eq_fs,
        fontweight="bold",
        color="#333",
    )

    ax.set_axis_off()


def main():
    figs_dir = Path(__file__).resolve().parent.parent / "figures"
    figs_dir.mkdir(parents=True, exist_ok=True)

    # --- Figure 1: BED operations ---
    fig_a, ax_a = plt.subplots(figsize=(11, 5.5))
    draw_bed_operations(ax_a)
    out_a = figs_dir / "exclusion_bed_operations.png"
    fig_a.savefig(out_a, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig_a)
    print(f"Saved to {out_a}")

    # --- Figure 2: Exclusion categories ---
    fig_b, ax_b = plt.subplots(figsize=(14, 6.5))
    draw_categories(ax_b)
    out_b = figs_dir / "exclusion_diagram.png"
    fig_b.savefig(out_b, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig_b)
    print(f"Saved to {out_b}")


if __name__ == "__main__":
    main()
