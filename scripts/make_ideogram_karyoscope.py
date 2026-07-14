#!/usr/bin/env python3
"""Alternative genome-view ideogram for the Q100 variant benchmark, rendered with
KaryoScope's "painted chromosome" SVG renderer.

AI Disclosure: Developed with assistance from Claude (Anthropic).

This is a visual *alternative* to ``scripts/make_ideogram.R`` (karyoploteR). It reuses
ONLY KaryoScope's decoupled renderer (``karyoscope.core.karyotype.render_karyotype``):
KaryoScope's k-mer annotation pipeline, its 17 GB database, KMC and the C++ helper are
NOT used. By default we feed the renderer the actual benchmark regions on GRCh38 (no
coverage binning; merged only at sub-pixel display resolution), draw telomere markers,
and rotate the result to landscape (chromosomes horizontal). ``--coverage-mode
presence|graded`` instead paints per-1 Mb-bin coverage.

Each chromosome is drawn as up to four parallel painted bars (one "contig" each,
all under a single haplotype so the bars carry no individual text label -- they are
identified by colour via the built-in legend):

  Difficult   SD10kb + TR10kb (large repeats)        grey  where covered
  v0.6 stvar  v0.6 GRCh37 -> liftOver GRCh38         teal  where covered
  v4.2.1 smvar v4.2.1 GRCh38                         orange where covered
  v5.0q       v5.0q smvar/stvar GRCh38               purple by category
                                                     (smvar-only / stvar-only / both)

v0.6 and v4.2.1 are omitted on chrX / chrY (they don't meaningfully cover the sex
chromosomes in HG002), so those rows show only Difficult + v5.0q.

Uncovered regions are labelled "Not in benchmark" (white) so every chromosome bar
renders at its true telomere-to-telomere length with white gaps at centromeres / hard
regions.

Usage:
  mamba run -n karyoscope python scripts/make_ideogram_karyoscope.py
  # fast iteration on one chromosome:
  mamba run -n karyoscope python scripts/make_ideogram_karyoscope.py --chroms chr1

Environment: see workflow/envs/karyoscope.yaml (provides karyoscope, drawsvg, cairosvg,
bedtools, ucsc-liftover, cairo). PDF/PNG conversion needs libcairo; SVG is always
written even if cairo is missing.

Output: manuscript/figs/ideogram_karyoscope.{svg,pdf,png}
"""

from __future__ import annotations

import argparse
import gzip
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from karyoscope.core.io.scaffold_map import MapRow
from karyoscope.core.karyotype import RenderInput, convert_svg, render_karyotype

# ============================================================================
# CONFIG
# ============================================================================
REPO = Path(__file__).resolve().parent.parent
RES_DIR = REPO / "resources"
BMK_DIR = RES_DIR / "benchmarksets"
STRAT_DIR = RES_DIR / "stratifications"
CHAIN_FILE = RES_DIR / "hg19ToHg38.over.chain.gz"
FIGS_DIR = REPO / "manuscript" / "figs"
OUTPUT_BASE = FIGS_DIR / "ideogram_karyoscope"

# Source for GRCh38 chromosome sizes (read from the VCF ##contig header).
CHROM_SIZES_VCF = BMK_DIR / "v5.0q_GRCh38_smvar_benchmark.vcf.gz"

# Benchmark + stratification BEDs (GRCh38 unless noted).
PATH_V421 = BMK_DIR / "v4.2.1_GRCh38_smvar_benchmark.bed"
PATH_SMVAR = BMK_DIR / "v5.0q_GRCh38_smvar_benchmark.bed"
PATH_STVAR = BMK_DIR / "v5.0q_GRCh38_stvar_benchmark.bed"
PATH_V06 = BMK_DIR / "v0.6_GRCh37_stvar_benchmark.bed"  # GRCh37, lifted below
# "Difficult" = large repeats only (segmental duplications + tandem repeats >= 10 kb).
# The full HP/TR/SD/MAP union is densely scattered (~20% of the genome in tiny pieces),
# so at this display resolution it paints almost every chromosome solid. The >=10 kb
# segdup + tandem-repeat regions (~8%) are localized (centromeres, acrocentric arms,
# large repeat arrays) and read as distinct hard-region blocks.
STRAT_BEDS = {
    "SD10kb": STRAT_DIR / "GRCh38_SD10kb.bed.gz",
    "TR10kb": STRAT_DIR / "GRCh38_TR10kb.bed.gz",
}

DEFAULT_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

# Per-chromosome track omissions. v0.6 and v4.2.1 don't meaningfully cover the sex
# chromosomes in HG002, so drop those bars on chrX / chrY (leaving Difficult + v5.0q).
EXCLUDE_TRACKS: dict[str, frozenset[str]] = {
    "chrX": frozenset({"v0.6", "v4.2.1"}),
    "chrY": frozenset({"v0.6", "v4.2.1"}),
}

# Render mode:
#   "regions"  : paint the actual benchmark intervals (no coverage binning) -- DEFAULT
#   "presence" : solid track colour per 1 Mb bin with >= COVERAGE_THRESHOLD covered
#   "graded"   : per 1 Mb bin, light->dark colour ramp by coverage fraction
COVERAGE_MODE = "regions"
# "regions" mode only: merge intervals separated by < MERGE_GAP bp before drawing.
# This is NOT coverage binning -- it only bridges gaps below the figure's pixel
# resolution (genome mode draws ~4 px/Mb, so 1 px ~= 250 kb) so the SVG stays small
# (raw "difficult" has ~4.8M intervals). Real gaps (centromeres) are far larger and
# are preserved.
MERGE_GAP = 50_000
# "presence"/"graded" modes only:
BIN_SIZE = 1_000_000
COVERAGE_THRESHOLD = 0.05  # min fraction of a 1 Mb bin covered to paint it
SHOW_TELOMERES = False  # draw KaryoScope telomere circles at chromosome ends
LANDSCAPE = True  # rotate the SVG 90deg so chromosomes lie horizontally
KEEP_LABELS_HORIZONTAL = True  # counter-rotate text so labels stay upright in landscape
SHOW_TITLE = True
SAMPLE_LABEL = "HG002 Q100 variant benchmark"
SEX = "male"  # HG002 is male; X/Y also appear because data is present

# Track order (left -> right within each chromosome). Internal keys; display labels
# and colours below. Palette matches scripts/make_ideogram.R.
TRACK_ORDER = ["difficult", "v0.6", "v4.2.1", "v5.0q"]
TRACK_LABELS = {
    "difficult": "Difficult >=10kb",
    "v0.6": "v0.6 stvar",
    "v4.2.1": "v4.2.1 smvar",
}
BASE_COLORS = {
    "difficult": "#888888",
    "v0.6": "#1B9E77",  # teal
    "v4.2.1": "#D95F02",  # orange
}
# v5.0q per-bin dominant category: internal key -> (display label, colour).
# (At 1 Mb resolution "both" almost always dominates; the only-categories surface
# only at finer bin sizes.)
V5_CATEGORIES = {
    "smvar_only": ("v5.0q smvar only", "#9E9AC8"),  # light purple
    "stvar_only": ("v5.0q stvar only", "#756BB1"),  # medium purple
    "both": ("v5.0q both", "#54278F"),  # dark purple
}
UNCOVERED_LABEL = "Not in benchmark"
UNCOVERED_COLOR = "#ffffff"

# Resolve external tools (conda env puts them on PATH; fall back to common paths).
BEDTOOLS = shutil.which("bedtools") or "/opt/homebrew/bin/bedtools"
LIFTOVER = shutil.which("liftOver")


# ============================================================================
# helpers
# ============================================================================
def log(msg: str) -> None:
    print(msg, file=sys.stderr, flush=True)


def run(cmd: str) -> None:
    """Run a shell command (bash, for process substitution etc.), raising on failure."""
    res = subprocess.run(
        cmd, shell=True, capture_output=True, text=True, executable="/bin/bash"
    )
    if res.returncode != 0:
        raise RuntimeError(
            f"command failed (rc={res.returncode}):\n{cmd}\n{res.stderr}"
        )


def read_chrom_sizes(vcf: Path, chroms: list[str]) -> dict[str, int]:
    """Read GRCh38 contig lengths from a VCF ##contig header (stops after the header)."""
    sizes: dict[str, int] = {}
    with gzip.open(vcf, "rt") as fh:
        for line in fh:
            if not line.startswith("##"):
                break
            if not line.startswith("##contig"):
                continue
            cid = _extract(line, "ID=")
            length = _extract(line, "length=")
            if cid in chroms and length is not None:
                sizes[cid] = int(length)
    missing = [c for c in chroms if c not in sizes]
    if missing:
        raise RuntimeError(f"chrom sizes missing from {vcf.name}: {missing}")
    return sizes


def _extract(contig_line: str, key: str) -> str | None:
    """Pull a ``key=value`` token out of a ``##contig=<...>`` header line."""
    try:
        rest = contig_line.split(key, 1)[1]
    except IndexError:
        return None
    return rest.split(",", 1)[0].rstrip(">\n").strip()


def write_genome_and_windows(
    sizes: dict[str, int], chroms: list[str], tmp: Path
) -> Path:
    """Write a bedtools genome file and tile it into BIN_SIZE windows. Returns windows BED."""
    genome = tmp / "genome.txt"
    genome.write_text("".join(f"{c}\t{sizes[c]}\n" for c in chroms))
    windows = tmp / "windows.bed"
    run(f"{BEDTOOLS} makewindows -g {genome} -w {BIN_SIZE} > {windows}")
    return windows


def prepare_tracks(tmp: Path) -> dict[str, Path]:
    """Build the GRCh38 BED for every painted track. Returns {track_key: bed_path}.

    Mirrors the bedtools/liftOver data prep in scripts/make_ideogram.R.
    """
    beds: dict[str, Path] = {}

    # Difficult = large repeats: merge(SD10kb, TR10kb) (see STRAT_BEDS).
    difficult = tmp / "difficult.bed"
    cat = "; ".join(f"gunzip -c {p}" for p in STRAT_BEDS.values())
    run(f"{{ {cat}; }} | sort -k1,1 -k2,2n | {BEDTOOLS} merge -i - > {difficult}")
    beds["difficult"] = difficult

    # v4.2.1 smvar -- GRCh38 native.
    beds["v4.2.1"] = PATH_V421

    # v5.0q category sub-tracks: smvar-only, stvar-only, both.
    smvar_only = tmp / "smvar_only.bed"
    stvar_only = tmp / "stvar_only.bed"
    both = tmp / "both.bed"
    run(f"{BEDTOOLS} subtract -a {PATH_SMVAR} -b {PATH_STVAR} > {smvar_only}")
    run(f"{BEDTOOLS} subtract -a {PATH_STVAR} -b {PATH_SMVAR} > {stvar_only}")
    run(f"{BEDTOOLS} intersect -a {PATH_SMVAR} -b {PATH_STVAR} > {both}")
    beds["smvar_only"] = smvar_only
    beds["stvar_only"] = stvar_only
    beds["both"] = both

    # v0.6 stvar -- GRCh37, lift to GRCh38.
    if LIFTOVER is None:
        log("WARNING: liftOver not found; skipping the v0.6 track.")
    else:
        v06 = lift_v06(tmp)
        if v06 is not None:
            beds["v0.6"] = v06
    return beds


def lift_v06(tmp: Path) -> Path | None:
    """Add a chr prefix to v0.6 (GRCh37), liftOver to GRCh38, sort + merge."""
    chr_bed = tmp / "v06_chr.bed"
    with PATH_V06.open() as src, chr_bed.open("w") as dst:
        for line in src:
            if not line.strip():
                continue
            f = line.split("\t")
            chrom = f[0] if f[0].startswith("chr") else f"chr{f[0]}"
            dst.write(f"{chrom}\t{f[1]}\t{f[2].rstrip()}\n")
    lifted = tmp / "v06_lifted.bed"
    unmapped = tmp / "v06_unmapped.bed"
    run(f"{LIFTOVER} {chr_bed} {CHAIN_FILE} {lifted} {unmapped}")
    merged = tmp / "v06_hg38.bed"
    run(f"sort -k1,1 -k2,2n {lifted} | {BEDTOOLS} merge -i - > {merged}")
    n_in = sum(1 for _ in chr_bed.open())
    n_out = sum(1 for _ in merged.open())
    log(f"v0.6 stvar: {n_in} intervals -> {n_out} after liftOver + merge")
    return merged


def coverage_fractions(windows: Path, track_bed: Path) -> dict[tuple[str, int], float]:
    """Fraction of each window covered by track_bed, keyed by (chrom, start)."""
    res = subprocess.run(
        f"{BEDTOOLS} coverage -a {windows} -b {track_bed}",
        shell=True,
        capture_output=True,
        text=True,
    )
    if res.returncode != 0:
        raise RuntimeError(f"bedtools coverage failed:\n{res.stderr}")
    out: dict[tuple[str, int], float] = {}
    for line in res.stdout.splitlines():
        f = line.split("\t")
        out[(f[0], int(f[1]))] = float(f[-1])
    return out


def load_regions(
    track_bed: Path, chroms: set[str], sizes: dict[str, int], merge_gap: int
) -> dict[str, list[tuple[int, int]]]:
    """Load actual track intervals per chromosome (regions mode).

    Intervals are merged at ``merge_gap`` (display resolution -- not coverage binning)
    and clipped to chromosome length. Returns ``{chrom: [(start, end), ...]}``.
    """
    res = subprocess.run(
        f"sort -k1,1 -k2,2n {track_bed} | {BEDTOOLS} merge -d {merge_gap} -i -",
        shell=True,
        capture_output=True,
        text=True,
        executable="/bin/bash",
    )
    if res.returncode != 0:
        raise RuntimeError(f"merge failed for {track_bed}:\n{res.stderr}")
    out: dict[str, list[tuple[int, int]]] = {}
    for line in res.stdout.splitlines():
        f = line.split("\t")
        chrom = f[0]
        if chrom not in chroms:
            continue
        start, end = int(f[1]), min(int(f[2]), sizes[chrom])
        if end > start:
            out.setdefault(chrom, []).append((start, end))
    return out


# ============================================================================
# colour ramp (graded mode)
# ============================================================================
def lighten(hex_color: str, amount: float) -> str:
    """Blend hex_color toward white. amount=0 -> original, amount=1 -> white."""
    r = int(hex_color[1:3], 16)
    g = int(hex_color[3:5], 16)
    b = int(hex_color[5:7], 16)
    r = round(r + (255 - r) * amount)
    g = round(g + (255 - g) * amount)
    b = round(b + (255 - b) * amount)
    return f"#{r:02x}{g:02x}{b:02x}"


# q1 (lowest coverage) -> q4 (highest); higher coverage = darker.
_GRADED_AMOUNTS = {1: 0.70, 2: 0.45, 3: 0.20, 4: 0.0}
_GRADED_BREAKS = [0.25, 0.50, 0.75]  # frac thresholds between q1..q4


def graded_bucket(frac: float) -> int:
    """Map a coverage fraction to a quantile bucket 1..4."""
    for i, brk in enumerate(_GRADED_BREAKS, start=1):
        if frac <= brk:
            return i
    return 4


# ============================================================================
# build renderer inputs
# ============================================================================
def build_render_inputs(
    chroms: list[str],
    sizes: dict[str, int],
    windows: Path,
    track_beds: dict[str, Path],
    coverage_mode: str,
) -> tuple[RenderInput, dict[str, str], list[str]]:
    """Construct the KaryoScope RenderInput, colour dict, and legend order."""
    if coverage_mode == "regions":
        return _build_regions(chroms, sizes, track_beds)
    return _build_coverage(chroms, sizes, windows, track_beds, coverage_mode)


def _present_tracks(track_beds: dict[str, Path]) -> list[str]:
    return [t for t in TRACK_ORDER if (t == "v5.0q" or t in track_beds)]


def _make_map_row(contig: str, chrom: str, sizes: dict[str, int]) -> MapRow:
    return MapRow(
        new_name=contig,
        original_name=contig,
        input_file="benchmark",
        hap="hap1",
        chromosome=chrom,
        flipped=False,
        length=sizes[chrom],
        stats="TT" if SHOW_TELOMERES else "",
    )


def _feature_order(
    present_tracks: list[str], coverage_mode: str, used_labels: set[str]
) -> list[str]:
    """Legend order: tracks left-to-right, v5.0q categories, then uncovered."""
    order: list[str] = []
    for track in present_tracks:
        if track == "v5.0q":
            order.extend(label for (label, _color) in V5_CATEGORIES.values())
        elif coverage_mode == "graded":
            order.extend(f"{TRACK_LABELS[track]} q{n}" for n in (1, 2, 3, 4))
        else:
            order.append(TRACK_LABELS[track])
    order.append(UNCOVERED_LABEL)
    return [f for f in order if f in used_labels]


def _build_regions(
    chroms: list[str], sizes: dict[str, int], track_beds: dict[str, Path]
) -> tuple[RenderInput, dict[str, str], list[str]]:
    """Region mode: paint the actual benchmark intervals (merged at MERGE_GAP)."""
    chroms_set = set(chroms)
    regions = {
        key: load_regions(track_beds[key], chroms_set, sizes, MERGE_GAP)
        for key in ("difficult", "v0.6", "v4.2.1", "smvar_only", "stvar_only", "both")
        if key in track_beds
    }

    colors: dict[str, str] = {UNCOVERED_LABEL: UNCOVERED_COLOR}
    used_labels: set[str] = {UNCOVERED_LABEL}
    map_rows: list[MapRow] = []
    binned: dict[str, list[tuple[int, int, str]]] = {}
    present_tracks = _present_tracks(track_beds)

    for chrom in chroms:
        for track in present_tracks:
            if track in EXCLUDE_TRACKS.get(chrom, frozenset()):
                continue
            contig = f"{chrom}__{track}"
            map_rows.append(_make_map_row(contig, chrom, sizes))
            # White full-length background first so the bar spans the whole
            # chromosome (sets length + telomere ends); regions drawn on top.
            intervals: list[tuple[int, int, str]] = [(0, sizes[chrom], UNCOVERED_LABEL)]
            cats = (
                ("both", "stvar_only", "smvar_only") if track == "v5.0q" else (track,)
            )
            for cat in cats:
                if cat not in regions:
                    continue
                if track == "v5.0q":
                    label, color = V5_CATEGORIES[cat]
                else:
                    label, color = TRACK_LABELS[track], BASE_COLORS[track]
                for start, end in regions[cat].get(chrom, []):
                    intervals.append((start, end, label))
                    colors[label] = color
                    used_labels.add(label)
            binned[contig] = intervals

    fo = _feature_order(present_tracks, "regions", used_labels)
    return RenderInput(map_rows=map_rows, binned_bed=binned), colors, fo


def _build_coverage(
    chroms: list[str],
    sizes: dict[str, int],
    windows: Path,
    track_beds: dict[str, Path],
    coverage_mode: str,
) -> tuple[RenderInput, dict[str, str], list[str]]:
    """Coverage mode (presence/graded): paint per-bin coverage of 1 Mb windows."""
    cov: dict[str, dict[tuple[str, int], float]] = {}
    binary_tracks = [t for t in ("difficult", "v0.6", "v4.2.1") if t in track_beds]
    for t in binary_tracks:
        cov[t] = coverage_fractions(windows, track_beds[t])
    for cat in ("smvar_only", "stvar_only", "both"):
        cov[cat] = coverage_fractions(windows, track_beds[cat])

    # Window list (chrom, start, end) restricted to requested chroms.
    win_by_chrom: dict[str, list[tuple[int, int]]] = {c: [] for c in chroms}
    for line in windows.read_text().splitlines():
        f = line.split("\t")
        if f[0] in win_by_chrom:
            win_by_chrom[f[0]].append((int(f[1]), int(f[2])))

    colors: dict[str, str] = {UNCOVERED_LABEL: UNCOVERED_COLOR}
    used_labels: set[str] = set()
    map_rows: list[MapRow] = []
    binned: dict[str, list[tuple[int, int, str]]] = {}
    present_tracks = _present_tracks(track_beds)

    for chrom in chroms:
        for track in present_tracks:
            if track in EXCLUDE_TRACKS.get(chrom, frozenset()):
                continue
            contig = f"{chrom}__{track}"
            map_rows.append(_make_map_row(contig, chrom, sizes))
            intervals: list[tuple[int, int, str]] = []
            for start, end in win_by_chrom[chrom]:
                label = _bin_label(track, chrom, start, cov, coverage_mode, colors)
                used_labels.add(label)
                intervals.append((start, end, label))
            binned[contig] = intervals

    fo = _feature_order(present_tracks, coverage_mode, used_labels)
    return RenderInput(map_rows=map_rows, binned_bed=binned), colors, fo


def _bin_label(
    track: str,
    chrom: str,
    start: int,
    cov: dict[str, dict[tuple[str, int], float]],
    coverage_mode: str,
    colors: dict[str, str],
) -> str:
    """Decide the feature label (and register its colour) for one bin of one track."""
    key = (chrom, start)
    if track == "v5.0q":
        fracs = {cat: cov[cat].get(key, 0.0) for cat in V5_CATEGORIES}
        cat, frac = max(fracs.items(), key=lambda kv: kv[1])
        if frac < COVERAGE_THRESHOLD:
            return UNCOVERED_LABEL
        label, color = V5_CATEGORIES[cat]
        colors[label] = color
        return label

    frac = cov[track].get(key, 0.0)
    if frac < COVERAGE_THRESHOLD:
        return UNCOVERED_LABEL
    base = BASE_COLORS[track]
    track_label = TRACK_LABELS[track]
    if coverage_mode == "graded":
        n = graded_bucket(frac)
        label = f"{track_label} q{n}"
        colors[label] = lighten(base, _GRADED_AMOUNTS[n])
        return label
    colors[track_label] = base
    return track_label


# ============================================================================
# landscape rotation
# ============================================================================
def rotate_svg_landscape(svg_path: Path) -> None:
    """Rotate a drawsvg SVG 90deg clockwise in place so the (vertical) chromosomes
    lie horizontally (chr1 ends up at the top). Swaps width/height and wraps all
    content in a rotation group, so PDF/PNG conversions inherit the rotation.
    """
    text = svg_path.read_text()
    m = re.search(r"<svg\b[^>]*>", text)
    if m is None:
        raise RuntimeError(f"could not find <svg> tag in {svg_path}")
    open_tag = m.group(0)
    w = int(re.search(r'width="(\d+)"', open_tag).group(1))
    h = int(re.search(r'height="(\d+)"', open_tag).group(1))
    new_tag = re.sub(r'width="\d+"', f'width="{h}"', open_tag)
    new_tag = re.sub(r'height="\d+"', f'height="{w}"', new_tag)
    new_tag = re.sub(r'viewBox="[^"]*"', f'viewBox="0 0 {h} {w}"', new_tag)
    # 90deg CW: point (x, y) -> (h - y, x). Columns (chromosomes) map top-to-bottom.
    g_open = f'<g transform="translate({h},0) rotate(90)">'
    text = text.replace(open_tag, f"{new_tag}\n{g_open}", 1)
    text = text.replace("</svg>", "</g>\n</svg>", 1)
    if KEEP_LABELS_HORIZONTAL:
        text = _counter_rotate_texts(text)
    svg_path.write_text(text)


def _counter_rotate_texts(svg_text: str) -> str:
    """Keep the per-chromosome labels upright after the 90deg group rotation.

    Only the chromosome labels are counter-rotated (identified by
    ``text-anchor="middle"`` and not bold). The legend (``text-anchor="start"``)
    and the bold title are laid out as vertical stacks, so they read fine rotated
    as a block; counter-rotating them in place would pile their entries on top of
    one another. Scale-bar labels already carry a transform and are left alone.
    """

    def repl(match: re.Match) -> str:
        attrs = match.group(1)
        if "transform=" in attrs:
            return match.group(0)
        if 'text-anchor="middle"' not in attrs or "font-weight" in attrs:
            return match.group(0)
        xm = re.search(r'\bx="(-?[\d.]+)"', attrs)
        ym = re.search(r'\by="(-?[\d.]+)"', attrs)
        if not (xm and ym):
            return match.group(0)
        return f'<text{attrs} transform="rotate(-90 {xm.group(1)} {ym.group(1)})">'

    return re.sub(r"<text\b([^>]*)>", repl, svg_text)


# ============================================================================
# main
# ============================================================================
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--chroms",
        default=",".join(DEFAULT_CHROMS),
        help="Comma-separated chromosomes (default: chr1..chr22,chrX,chrY).",
    )
    ap.add_argument(
        "--coverage-mode",
        choices=("regions", "presence", "graded"),
        default=COVERAGE_MODE,
        help="regions (actual intervals, no binning), presence (solid 1 Mb bins), "
        "or graded (1 Mb coverage ramp).",
    )
    ap.add_argument("--output-base", default=str(OUTPUT_BASE), help="Output path stem.")
    args = ap.parse_args()

    chroms = [c.strip() for c in args.chroms.split(",") if c.strip()]
    out_base = Path(args.output_base)
    out_base.parent.mkdir(parents=True, exist_ok=True)

    log(f"Reading chromosome sizes from {CHROM_SIZES_VCF.name} ...")
    sizes = read_chrom_sizes(CHROM_SIZES_VCF, chroms)

    with tempfile.TemporaryDirectory(prefix="ks_ideogram_") as td:
        tmp = Path(td)
        log("Tiling windows ...")
        windows = write_genome_and_windows(sizes, chroms, tmp)
        log("Preparing track BEDs (bedtools + liftOver) ...")
        track_beds = prepare_tracks(tmp)
        log(f"Building render inputs ({args.coverage_mode} mode) ...")
        render_input, colors, feature_order = build_render_inputs(
            chroms, sizes, windows, track_beds, args.coverage_mode
        )

    n_contigs = len(render_input.map_rows)
    n_rects = sum(len(v) for v in render_input.binned_bed.values())
    log(
        f"Built {n_contigs} painted bars across {len(chroms)} chromosomes "
        f"({n_rects} rectangles)."
    )
    log(f"Legend features: {feature_order}")

    # seed_human_chromosomes=False so only the requested chroms appear (no empty
    # columns when rendering a subset for fast iteration).
    svg_path = out_base.with_suffix(".svg")
    log(f"Rendering {svg_path} ...")
    render_karyotype(
        [render_input],
        colors=colors,
        mode="genome",
        sex=SEX,
        background_color="white",
        seed_human_chromosomes=False,
        output_path=svg_path,
        sample_label=SAMPLE_LABEL,
        database_id=None,
        feature_set_label="benchmark coverage",
        smoothed=False,
        show_title=SHOW_TITLE,
        show_legend=True,
        feature_order=feature_order,
    )

    if LANDSCAPE:
        log("Rotating to landscape (chromosomes horizontal) ...")
        rotate_svg_landscape(svg_path)

    for ext in (".pdf", ".png"):
        target = out_base.with_suffix(ext)
        try:
            convert_svg(svg_path, target)
            log(f"Wrote {target} ({target.stat().st_size:,} bytes)")
        except Exception as exc:  # noqa: BLE001 - cairo may be unavailable
            log(f"WARNING: could not write {target.name}: {exc}")
    log(f"Wrote {svg_path} ({svg_path.stat().st_size:,} bytes)")
    log("Done.")


if __name__ == "__main__":
    main()
