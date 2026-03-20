#!/usr/bin/env python3
"""Extract small variant benchmarking metrics from vcfeval ROC output files.

Parses ROC TSV files produced by rtg vcfeval with --roc-regions and --roc-subset
flags. Produces CSV files for use by analysis/use_case_evaluation.qmd.

ROC file types:
  - 8-column (weighted, snp, non_snp): includes FN, precision, sensitivity, f_measure
  - 4-column (per-region: hp, tr, sd, map): includes only TP, FP, TP_call

AI Disclosure: Developed with assistance from Claude Opus 4.6 (Anthropic).
"""

import csv
import gzip
import subprocess
import sys
from pathlib import Path

PROJ_DIR = Path(__file__).resolve().parent.parent
SMVAR_DIR = PROJ_DIR / "results" / "use_case" / "smvar"

BENCH_VERSIONS = ["v4.2.1", "v5.0q"]

# ROC files with 8 columns (full metrics)
FULL_ROC_FILES = {
    "Overall": "weighted_roc.tsv.gz",
    "SNP": "snp_roc.tsv.gz",
    "NON_SNP": "non_snp_roc.tsv.gz",
}

# ROC files with 4 columns (TP + FP only)
REGION_ROC_FILES = {
    "HP": "hp_roc.tsv.gz",
    "TR": "tr_roc.tsv.gz",
    "SD": "sd_roc.tsv.gz",
    "MAP": "map_roc.tsv.gz",
}

# Size bin definitions for post-hoc counting from TP/FP/FN VCFs
SIZE_BINS = [
    (0, 1, "SNP"),
    (1, 5, "1-4bp"),
    (5, 15, "5-14bp"),
    (15, 50, "15-49bp"),
]


def parse_roc_last_row(roc_path: Path) -> dict[str, float | None]:
    """Parse the last data row from a vcfeval ROC TSV file.

    Returns a dict with keys matching the column headers.
    For 4-column files, FN/precision/sensitivity/f_measure will be None.
    """
    last_line = ""
    total_baseline = None
    with gzip.open(roc_path, "rt") as f:
        for line in f:
            if line.startswith("#total baseline variants:"):
                total_baseline = int(line.split(":")[1].strip())
            if not line.startswith("#"):
                last_line = line.strip()

    if not last_line:
        print(f"WARNING: No data rows in {roc_path}", file=sys.stderr)
        return {}

    fields = last_line.split("\t")

    if len(fields) >= 8:
        # 8-column format: score, tp_baseline, fp, tp_call, fn, precision, sensitivity, f_measure
        return {
            "tp_baseline": float(fields[1]),
            "fp": float(fields[2]),
            "tp_call": float(fields[3]),
            "fn": float(fields[4]),
            "precision": float(fields[5]),
            "sensitivity": float(fields[6]),
            "f_measure": float(fields[7]),
            "total_baseline": total_baseline,
        }
    elif len(fields) >= 4:
        # 4-column format: score, tp_baseline, fp, tp_call
        return {
            "tp_baseline": float(fields[1]),
            "fp": float(fields[2]),
            "tp_call": float(fields[3]),
            "fn": None,
            "precision": None,
            "sensitivity": None,
            "f_measure": None,
            "total_baseline": total_baseline,
        }
    else:
        print(f"WARNING: Unexpected column count in {roc_path}: {len(fields)}", file=sys.stderr)
        return {}


def get_size_bin_counts(vcf_path: Path) -> dict[str, int]:
    """Count variants by size bin from a VCF file."""
    result = subprocess.run(
        ["bcftools", "query", "-f", "%REF\t%ALT\n", str(vcf_path)],
        capture_output=True,
        text=True,
        check=True,
    )
    counts: dict[str, int] = {}
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        parts = line.split("\t")
        ref_len = len(parts[0])
        alt_len = len(parts[1])
        asize = abs(alt_len - ref_len)

        size_bin = ">=50bp"
        for lo, hi, label in SIZE_BINS:
            if lo <= asize < hi:
                size_bin = label
                break

        counts[size_bin] = counts.get(size_bin, 0) + 1
    return counts


def main() -> None:
    # Verify run directories exist
    for bench in BENCH_VERSIONS:
        run_dir = SMVAR_DIR / f"roche_{bench}"
        if not run_dir.exists():
            print(
                f"ERROR: {run_dir} does not exist. Run scripts/run_smvar_use_case.sh first.",
                file=sys.stderr,
            )
            sys.exit(1)

    # === 1. ROC metrics (all subsets) ===
    roc_rows: list[dict[str, str | float | None]] = []

    for bench in BENCH_VERSIONS:
        run_dir = SMVAR_DIR / f"roche_{bench}"

        # 8-column files (full metrics)
        for subset, filename in FULL_ROC_FILES.items():
            roc_path = run_dir / filename
            if not roc_path.exists():
                print(f"WARNING: {roc_path} not found, skipping", file=sys.stderr)
                continue
            metrics = parse_roc_last_row(roc_path)
            if metrics:
                roc_rows.append({"bench_version": bench, "subset": subset, **metrics})

        # 4-column files (per-region)
        for subset, filename in REGION_ROC_FILES.items():
            roc_path = run_dir / filename
            if not roc_path.exists():
                print(f"WARNING: {roc_path} not found, skipping", file=sys.stderr)
                continue
            metrics = parse_roc_last_row(roc_path)
            if metrics:
                roc_rows.append({"bench_version": bench, "subset": subset, **metrics})

    roc_path = SMVAR_DIR / "vcfeval_roc_metrics.csv"
    fieldnames = [
        "bench_version",
        "subset",
        "tp_baseline",
        "fp",
        "tp_call",
        "fn",
        "precision",
        "sensitivity",
        "f_measure",
        "total_baseline",
    ]
    with open(roc_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(roc_rows)
    print(f"Wrote {roc_path}")

    # === 2. Size bin counts (post-hoc from TP/FP/FN VCFs) ===
    size_rows: list[dict[str, str | int]] = []

    for bench in BENCH_VERSIONS:
        run_dir = SMVAR_DIR / f"roche_{bench}"
        vcf_map = {
            "TP": run_dir / "tp-baseline.vcf.gz",
            "FP": run_dir / "fp.vcf.gz",
            "FN": run_dir / "fn.vcf.gz",
        }
        for category, vcf_path in vcf_map.items():
            counts = get_size_bin_counts(vcf_path)
            for size_bin, count in sorted(counts.items()):
                size_rows.append({
                    "bench_version": bench,
                    "category": category,
                    "size_bin": size_bin,
                    "count": count,
                })

    size_path = SMVAR_DIR / "size_bin_counts.csv"
    with open(size_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["bench_version", "category", "size_bin", "count"])
        writer.writeheader()
        writer.writerows(size_rows)
    print(f"Wrote {size_path}")


if __name__ == "__main__":
    main()
