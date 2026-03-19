#!/usr/bin/env python3
"""Extract SV benchmarking metrics from truvari bench+refine output.

Reads post-refine truvari output for both SV callsets and generates three
CSV files matching the format expected by use_case_evaluation.qmd.

Requires: bcftools, bedtools (available in truvari conda env)

AI Disclosure: Developed with assistance from Claude Opus 4.6 (Anthropic).
"""

import csv
import subprocess
import sys
from pathlib import Path

PROJ_DIR = Path(__file__).resolve().parent.parent
STVAR_DIR = PROJ_DIR / "results" / "use_case" / "stvar"
STRAT_DIR = PROJ_DIR / "resources" / "stratifications"

CALLSETS = ["baylor_ont", "dragen"]
CONTEXTS = ["HP", "TR", "SD", "MAP"]
SIZE_BINS = [
    (50, 100, "50-99bp"),
    (100, 300, "100-299bp"),
    (300, 1000, "300-999bp"),
    (1000, 10000, "1-10kb"),
    (10000, float("inf"), ">=10kb"),
]


def get_refine_vcfs(run_dir: Path) -> dict[str, Path]:
    """Return paths to post-refine TP/FP/FN VCFs.

    Truvari refine updates tp-base.vcf.gz, fp.vcf.gz, fn.vcf.gz in-place.
    Note: refine.base.vcf.gz contains ALL base variants (TP+FN), not just TPs.
    """
    return {
        "TP": run_dir / "tp-base.vcf.gz",
        "FP": run_dir / "fp.vcf.gz",
        "FN": run_dir / "fn.vcf.gz",
    }


def count_vcf_records(vcf_path: Path) -> int:
    """Count records in a VCF file."""
    result = subprocess.run(
        ["bcftools", "view", "-H", str(vcf_path)],
        capture_output=True,
        text=True,
        check=True,
    )
    return len(result.stdout.strip().split("\n")) if result.stdout.strip() else 0


def count_intersecting_variants(vcf_path: Path, bed_path: Path) -> int:
    """Count VCF records overlapping a BED file using bedtools."""
    result = subprocess.run(
        ["bedtools", "intersect", "-u", "-a", str(vcf_path), "-b", str(bed_path)],
        capture_output=True,
        text=True,
        check=True,
    )
    # Filter to data lines (skip headers)
    lines = [l for l in result.stdout.strip().split("\n") if l and not l.startswith("#")]
    return len(lines)


def get_svtype_counts(vcf_path: Path) -> dict[str, int]:
    """Count variants by SVTYPE from a VCF."""
    result = subprocess.run(
        ["bcftools", "query", "-f", "%INFO/SVTYPE\n", str(vcf_path)],
        capture_output=True,
        text=True,
        check=True,
    )
    counts: dict[str, int] = {}
    for line in result.stdout.strip().split("\n"):
        if line:
            svtype = line.strip()
            counts[svtype] = counts.get(svtype, 0) + 1
    return counts


def get_svtype_size_counts(vcf_path: Path) -> dict[tuple[str, str], int]:
    """Count variants by SVTYPE and size bin."""
    result = subprocess.run(
        ["bcftools", "query", "-f", "%INFO/SVTYPE\t%INFO/SVLEN\n", str(vcf_path)],
        capture_output=True,
        text=True,
        check=True,
    )
    counts: dict[tuple[str, str], int] = {}
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        parts = line.strip().split("\t")
        if len(parts) < 2:
            continue
        svtype = parts[0]
        try:
            svlen = abs(int(parts[1]))
        except ValueError:
            svlen = 50  # default for missing SVLEN
        if svlen == 0:
            svlen = 50

        size_bin = ">=10kb"  # default
        for lo, hi, label in SIZE_BINS:
            if lo <= svlen < hi:
                size_bin = label
                break

        key = (svtype, size_bin)
        counts[key] = counts.get(key, 0) + 1
    return counts


def compute_metrics(tp: int, fp: int, fn: int) -> tuple[float, float]:
    """Compute recall and precision."""
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    return recall, precision


def main() -> None:
    # Verify all run directories exist
    for callset in CALLSETS:
        run_dir = STVAR_DIR / f"{callset}_v5.0q"
        if not run_dir.exists():
            print(f"ERROR: {run_dir} does not exist. Run scripts/run_sv_use_case.sh first.", file=sys.stderr)
            sys.exit(1)

    # === 1. Stratified metrics ===
    strat_rows: list[dict[str, str | int | float]] = []

    for callset in CALLSETS:
        run_dir = STVAR_DIR / f"{callset}_v5.0q"
        vcfs = get_refine_vcfs(run_dir)

        # Overall counts from VCF records (consistent with per-context counts)
        tp_overall = count_vcf_records(vcfs["TP"])
        fp_overall = count_vcf_records(vcfs["FP"])
        fn_overall = count_vcf_records(vcfs["FN"])
        recall, precision = compute_metrics(tp_overall, fp_overall, fn_overall)
        strat_rows.append({
            "callset": callset, "context": "Overall",
            "tp": tp_overall, "fp": fp_overall, "fn": fn_overall,
            "recall": round(recall, 4), "precision": round(precision, 4),
        })

        # Per-context counts via bedtools intersect
        for ctx in CONTEXTS:
            bed = STRAT_DIR / f"GRCh38_{ctx}.bed.gz"
            tp = count_intersecting_variants(vcfs["TP"], bed)
            fp = count_intersecting_variants(vcfs["FP"], bed)
            fn = count_intersecting_variants(vcfs["FN"], bed)
            recall, precision = compute_metrics(tp, fp, fn)
            strat_rows.append({
                "callset": callset, "context": ctx,
                "tp": tp, "fp": fp, "fn": fn,
                "recall": round(recall, 4), "precision": round(precision, 4),
            })

    strat_path = STVAR_DIR / "stratified_metrics.csv"
    with open(strat_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["callset", "context", "tp", "fp", "fn", "recall", "precision"])
        writer.writeheader()
        writer.writerows(strat_rows)
    print(f"Wrote {strat_path}")

    # === 2. SV type metrics ===
    svtype_rows: list[dict[str, str | int | float]] = []

    for callset in CALLSETS:
        run_dir = STVAR_DIR / f"{callset}_v5.0q"
        vcfs = get_refine_vcfs(run_dir)
        tp_counts = get_svtype_counts(vcfs["TP"])
        fp_counts = get_svtype_counts(vcfs["FP"])
        fn_counts = get_svtype_counts(vcfs["FN"])

        all_types = sorted(set(tp_counts) | set(fp_counts) | set(fn_counts))
        for svtype in all_types:
            tp = tp_counts.get(svtype, 0)
            fp = fp_counts.get(svtype, 0)
            fn = fn_counts.get(svtype, 0)
            recall, precision = compute_metrics(tp, fp, fn)
            svtype_rows.append({
                "callset": callset, "svtype": svtype,
                "tp": tp, "fp": fp, "fn": fn,
                "recall": round(recall, 4), "precision": round(precision, 4),
            })

    svtype_path = STVAR_DIR / "svtype_metrics.csv"
    with open(svtype_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["callset", "svtype", "tp", "fp", "fn", "recall", "precision"])
        writer.writeheader()
        writer.writerows(svtype_rows)
    print(f"Wrote {svtype_path}")

    # === 3. SV type x size bin counts ===
    size_rows: list[dict[str, str | int]] = []

    for callset in CALLSETS:
        run_dir = STVAR_DIR / f"{callset}_v5.0q"
        vcfs = get_refine_vcfs(run_dir)
        for category, vcf_path in vcfs.items():
            counts = get_svtype_size_counts(vcf_path)
            for (svtype, size_bin), count in sorted(counts.items()):
                size_rows.append({
                    "callset": callset, "category": category,
                    "svtype": svtype, "size_bin": size_bin,
                    "count": count,
                })

    size_path = STVAR_DIR / "svtype_size_counts.csv"
    with open(size_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["callset", "category", "svtype", "size_bin", "count"])
        writer.writeheader()
        writer.writerows(size_rows)
    print(f"Wrote {size_path}")


if __name__ == "__main__":
    main()
