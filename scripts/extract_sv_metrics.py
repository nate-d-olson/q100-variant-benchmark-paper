#!/usr/bin/env python3
"""Extract SV benchmarking metrics from truvari bench+refine output.

Reads post-refine truvari output for 4 callset configurations (HiFi/ONT ×
Sniffles1/Sniffles2) from the data/ directory and generates three CSV files.

Runs `truvari stratify` for per-context metrics using stratification BEDs.

Requires: truvari, pandas (available in truvari conda env)

AI Disclosure: Developed with assistance from Claude Opus 4.6 (Anthropic).
"""

import csv
import json
import sys
from pathlib import Path

import pandas as pd
import truvari
from truvari.stratify import stratify_main

PROJ_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = PROJ_DIR / "data"
STRAT_DIR = PROJ_DIR / "resources" / "stratifications"
OUTPUT_DIR = PROJ_DIR / "results" / "use_case" / "stvar"

CALLSET_CONFIG: dict[tuple[str, str], str] = {
    ("HiFi", "sniffles1"): "GRCh38_HG002_T~T2TQ100v1.1_Q~HiFi-Snifpub-sniffles1_TR~none_GRCh38_HG002-T2TQ100v1.1-dipz2k_stvar-excluded",
    ("HiFi", "sniffles2"): "GRCh38_HG002_T~T2TQ100v1.1_Q~HiFi-Snifpub-sniffles2_TR~none_GRCh38_HG002-T2TQ100v1.1-dipz2k_stvar-excluded",
    ("ONT", "sniffles1"): "GRCh38_HG002_T~T2TQ100v1.1_Q~ONT-Snifpub-sniffles1_TR~none_GRCh38_HG002-T2TQ100v1.1-dipz2k_stvar-excluded",
    ("ONT", "sniffles2"): "GRCh38_HG002_T~T2TQ100v1.1_Q~ONT-Snifpub-sniffles2_TR~none_GRCh38_HG002-T2TQ100v1.1-dipz2k_stvar-excluded",
}

CONTEXTS = ["HP", "TR", "SD", "MAP"]

# Map truvari's fine-grained size bins to manuscript-level bins
SIZEBIN_MAP: dict[str, str] = {
    "[50,100)": "[50,100)",
    "[100,200)": "[100,300)",
    "[200,300)": "[100,300)",
    "[300,400)": "[300,1000)",
    "[400,600)": "[300,1000)",
    "[600,800)": "[300,1000)",
    "[800,1k)": "[300,1000)",
    "[1k,2.5k)": "[1000,10000)",
    "[2.5k,5k)": "[1000,10000)",
    ">=5k": ">=10000",
}


def find_truvari_dir(top_dir: Path) -> Path:
    """Find the truvari output subdirectory (skip phab_bench/)."""
    for p in top_dir.iterdir():
        if p.is_dir() and p.name != "phab_bench":
            if (p / "refine.variant_summary.json").exists():
                return p
    raise FileNotFoundError(f"No truvari output dir with refine.variant_summary.json in {top_dir}")


def run_truvari_stratify(bed_path: Path, truvari_dir: Path, output_path: Path) -> None:
    """Run truvari stratify on a truvari output directory via Python API."""
    stratify_main([
        str(bed_path),
        str(truvari_dir),
        "-o", str(output_path),
    ])


def get_refine_vcfs(run_dir: Path) -> dict[str, Path]:
    """Return paths to post-refine TP/FP/FN VCFs."""
    return {
        "TP": run_dir / "tp-base.vcf.gz",
        "FP": run_dir / "fp.vcf.gz",
        "FN": run_dir / "fn.vcf.gz",
    }


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Resolve all truvari dirs upfront
    truvari_dirs: dict[tuple[str, str], Path] = {}
    for (platform, callset), dirname in CALLSET_CONFIG.items():
        top_dir = DATA_DIR / dirname
        if not top_dir.exists():
            print(f"ERROR: {top_dir} does not exist.", file=sys.stderr)
            sys.exit(1)
        truvari_dirs[(platform, callset)] = find_truvari_dir(top_dir)

    # === 1. Stratified metrics (Overall + per-context) ===
    strat_rows: list[dict[str, str | int | float]] = []

    for (platform, callset), run_dir in truvari_dirs.items():
        # Overall metrics from refine.variant_summary.json
        with open(run_dir / "refine.variant_summary.json") as f:
            summary = json.load(f)

        strat_rows.append({
            "platform": platform,
            "callset": callset,
            "context": "Overall",
            "tp": summary["TP-comp"],
            "fp": summary["FP"],
            "fn": summary["FN"],
            "recall": round(summary["recall"], 4),
            "precision": round(summary["precision"], 4),
        })

        # Per-context: run truvari stratify
        for context in CONTEXTS:
            bed_path = STRAT_DIR / f"GRCh38_{context}.bed.gz"
            if not bed_path.exists():
                print(f"WARNING: {bed_path} not found, skipping {context}", file=sys.stderr)
                continue

            strat_out = OUTPUT_DIR / f"stratify_{platform}_{callset}_{context}.txt"
            run_truvari_stratify(bed_path, run_dir, strat_out)

            strat_df = pd.read_csv(
                strat_out, sep="\t",
                names=["chrom", "start", "end", "tpbase", "tp", "fn", "fp"],
            )
            totals = strat_df[["tpbase", "tp", "fn", "fp"]].sum()

            precision, recall, f1 = truvari.performance_metrics(
                int(totals["tpbase"]),
                int(totals["tp"]),
                int(totals["fn"]),
                int(totals["fp"]),
            )

            strat_rows.append({
                "platform": platform,
                "callset": callset,
                "context": context,
                "tp": int(totals["tp"]),
                "fp": int(totals["fp"]),
                "fn": int(totals["fn"]),
                "recall": round(recall, 4),
                "precision": round(precision, 4),
            })

    strat_path = OUTPUT_DIR / "stratified_metrics.csv"
    with open(strat_path, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=["platform", "callset", "context", "tp", "fp", "fn", "recall", "precision"]
        )
        writer.writeheader()
        writer.writerows(strat_rows)
    print(f"Wrote {strat_path}")

    # === 2. SV type metrics ===
    svtype_rows: list[dict[str, str | int | float]] = []

    for (platform, callset), run_dir in truvari_dirs.items():
        vcfs = get_refine_vcfs(run_dir)

        tp_df = truvari.vcf_to_df(str(vcfs["TP"]))
        fp_df = truvari.vcf_to_df(str(vcfs["FP"]))
        fn_df = truvari.vcf_to_df(str(vcfs["FN"]))

        all_types = sorted(set(tp_df["svtype"]) | set(fp_df["svtype"]) | set(fn_df["svtype"]))

        for svtype in all_types:
            tp = len(tp_df[tp_df["svtype"] == svtype])
            fp = len(fp_df[fp_df["svtype"] == svtype])
            fn = len(fn_df[fn_df["svtype"] == svtype])

            precision, recall, f1 = truvari.performance_metrics(tp, tp, fn, fp)

            svtype_rows.append({
                "platform": platform,
                "callset": callset,
                "svtype": svtype,
                "tp": tp,
                "fp": fp,
                "fn": fn,
                "recall": round(recall, 4) if recall is not None else 0.0,
                "precision": round(precision, 4) if precision is not None else 0.0,
            })

    svtype_path = OUTPUT_DIR / "svtype_metrics.csv"
    with open(svtype_path, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=["platform", "callset", "svtype", "tp", "fp", "fn", "recall", "precision"]
        )
        writer.writeheader()
        writer.writerows(svtype_rows)
    print(f"Wrote {svtype_path}")

    # === 3. SV type x size bin counts ===
    size_rows: list[dict[str, str | int]] = []

    for (platform, callset), run_dir in truvari_dirs.items():
        vcfs = get_refine_vcfs(run_dir)

        for category, vcf_path in vcfs.items():
            df = truvari.vcf_to_df(str(vcf_path))
            df["size_bin"] = df["svlen"].apply(
                lambda x: SIZEBIN_MAP.get(
                    truvari.get_sizebin(abs(x) if x else 0), "other"
                )
            )
            df = df[df["size_bin"] != "other"]

            counts = df.groupby(["svtype", "size_bin"]).size().reset_index(name="count")
            counts["platform"] = platform
            counts["callset"] = callset
            counts["category"] = category

            size_rows.extend(counts.to_dict("records"))

    size_path = OUTPUT_DIR / "svtype_size_counts.csv"
    with open(size_path, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=["platform", "callset", "category", "svtype", "size_bin", "count"]
        )
        writer.writeheader()
        writer.writerows(size_rows)
    print(f"Wrote {size_path}")


if __name__ == "__main__":
    main()
