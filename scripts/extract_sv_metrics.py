#!/usr/bin/env python3
"""Extract SV benchmarking metrics from truvari bench+refine output.

Reads post-refine truvari output for both SV callsets and generates three
CSV files matching the format expected by use_case_evaluation.qmd.

Requires: truvari, pandas (available in truvari conda env)

AI Disclosure: Developed with assistance from Claude Opus 4.6 (Anthropic).
"""

import csv
import json
import sys
from pathlib import Path

import pandas as pd
import truvari

PROJ_DIR = Path(__file__).resolve().parent.parent
STVAR_DIR = PROJ_DIR / "results" / "use_case" / "stvar"
STRAT_DIR = PROJ_DIR / "resources" / "stratifications"

CALLSETS = ["ont-sniffles", "ont-verkko"]
CONTEXTS = ["HP", "TR", "SD", "MAP"]


def get_refine_vcfs(run_dir: Path) -> dict[str, Path]:
    """Return paths to post-refine TP/FP/FN VCFs."""
    return {
        "TP": run_dir / "tp-base.vcf.gz",
        "FP": run_dir / "fp.vcf.gz",
        "FN": run_dir / "fn.vcf.gz",
    }


def main() -> None:
    # Verify all run directories exist
    for callset in CALLSETS:
        run_dir = STVAR_DIR / f"{callset}_v5.0q"
        if not run_dir.exists():
            print(f"ERROR: {run_dir} does not exist. Run scripts/run_sv_use_case.sh first.", file=sys.stderr)
            sys.exit(1)

    # === 1. Stratified metrics (Overall + per-context) ===
    strat_rows: list[dict[str, str | int | float]] = []

    for callset in CALLSETS:
        run_dir = STVAR_DIR / f"{callset}_v5.0q"

        # Overall metrics from summary.json
        with open(run_dir / "summary.json") as f:
            summary = json.load(f)

        precision, recall, f1 = truvari.performance_metrics(
            summary["TP-base"],
            summary["TP-call"],
            summary["FN"],
            summary["FP"]
        )

        strat_rows.append({
            "callset": callset,
            "context": "Overall",
            "tp": summary["TP-call"],
            "fp": summary["FP"],
            "fn": summary["FN"],
            "recall": round(recall, 4),
            "precision": round(precision, 4),
        })

        # Per-context metrics from stratify text files
        for context in CONTEXTS:
            strat_file = run_dir / f"stratify_{context}.txt"
            strat_df = pd.read_csv(
                strat_file,
                sep='\t',
                names=['chrom', 'start', 'end', 'tpbase', 'tp', 'fn', 'fp']
            )

            # Sum across all regions in this context
            totals = strat_df[["tpbase", "tp", "fn", "fp"]].sum()

            precision, recall, f1 = truvari.performance_metrics(
                int(totals["tpbase"]),
                int(totals["tp"]),
                int(totals["fn"]),
                int(totals["fp"])
            )

            strat_rows.append({
                "callset": callset,
                "context": context,
                "tp": int(totals["tp"]),
                "fp": int(totals["fp"]),
                "fn": int(totals["fn"]),
                "recall": round(recall, 4),
                "precision": round(precision, 4),
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

        # Load VCFs via Truvari API
        tp_df = truvari.vcf_to_df(str(vcfs["TP"]))
        fp_df = truvari.vcf_to_df(str(vcfs["FP"]))
        fn_df = truvari.vcf_to_df(str(vcfs["FN"]))

        # Get all SV types present in any category
        all_types = sorted(set(tp_df["svtype"]) | set(fp_df["svtype"]) | set(fn_df["svtype"]))

        for svtype in all_types:
            tp = len(tp_df[tp_df["svtype"] == svtype])
            fp = len(fp_df[fp_df["svtype"] == svtype])
            fn = len(fn_df[fn_df["svtype"] == svtype])

            precision, recall, f1 = truvari.performance_metrics(tp, tp, fn, fp)

            svtype_rows.append({
                "callset": callset,
                "svtype": svtype,
                "tp": tp,
                "fp": fp,
                "fn": fn,
                "recall": round(recall, 4),
                "precision": round(precision, 4),
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
            # Load VCF and add size_bin column
            df = truvari.vcf_to_df(str(vcf_path))
            df["size_bin"] = df["svlen"].apply(
                lambda x: truvari.get_sizebin(abs(x) if x else 0)
            )

            # Count by (svtype, size_bin)
            counts = df.groupby(["svtype", "size_bin"]).size().reset_index(name="count")
            counts["callset"] = callset
            counts["category"] = category

            size_rows.extend(counts.to_dict("records"))

    size_path = STVAR_DIR / "svtype_size_counts.csv"
    with open(size_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["callset", "category", "svtype", "size_bin", "count"])
        writer.writeheader()
        writer.writerows(size_rows)
    print(f"Wrote {size_path}")


if __name__ == "__main__":
    main()
