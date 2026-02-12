#!/usr/bin/env python3
# AI Disclosure: This script was developed with assistance from Claude (Anthropic).
"""
Count variants per exclusion from Parquet variant table.

Reads the Parquet variant table (already size-filtered and type-classified),
parses region_ids for EXCL_* tags, and counts variants per exclusion.
Merges with per-exclusion BED metrics to produce a combined impact table.
"""

import csv
import logging
from collections import Counter
from typing import Dict

import pandas as pd


def count_exclusion_variants(
    parquet_path: str,
    exclusion_tsvs: list,
    output_csv: str,
    bench_type: str,
    excl_name_mapping: Dict[str, str],
    log_path: str,
) -> None:
    """Count variants per exclusion from Parquet table.

    Size filtering and variant type classification are already done
    in the Parquet table.

    Args:
        parquet_path: Path to variants.parquet
        exclusion_tsvs: Paths to per-exclusion BED metric TSV files
        output_csv: Path to output CSV
        bench_type: "smvar" or "stvar"
        excl_name_mapping: Maps EXCL_* IDs to canonical names
        log_path: Path to log file
    """
    is_smvar = bench_type == "smvar"

    logging.basicConfig(
        filename=log_path,
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )

    logging.info("Counting exclusion variants for %s", bench_type)
    logging.info("Input: %s", parquet_path)
    logging.info("Exclusion mapping: %s", excl_name_mapping)

    # Read only needed columns from Parquet
    df = pd.read_parquet(parquet_path, columns=["region_ids", "var_type"])
    logging.info("Loaded %d variants", len(df))

    # Filter to variants with region_ids
    # Empty strings are already converted to NA in generate_variant_parquet.py
    df = df[df["region_ids"].notna()].copy()

    # Explode region_ids and filter to EXCL_* only
    df["region_ids"] = df["region_ids"].str.split(",")
    df_exploded = df.explode("region_ids")
    df_exploded["region_ids"] = df_exploded["region_ids"].str.strip()
    df_excl = df_exploded[df_exploded["region_ids"].str.startswith("EXCL_")].copy()

    logging.info("Found %d variant-exclusion pairs", len(df_excl))

    # Count by (excl_id, var_type)
    counts_df = (
        df_excl.groupby(["region_ids", "var_type"]).size().reset_index(name="count")
    )

    # Build counts dict: {excl_id: {var_type: count}}
    counts: Dict[str, Counter] = {eid: Counter() for eid in excl_name_mapping}
    for _, row in counts_df.iterrows():
        eid = row["region_ids"]
        if eid in counts:
            counts[eid][row["var_type"]] = row["count"]

    # Read per-exclusion BED metrics
    bed_metrics: Dict[str, Dict[str, str]] = {}
    for tsv_path in exclusion_tsvs:
        with open(tsv_path, "r") as tf:
            for line in tf:
                parts = line.strip().split("\t")
                if len(parts) >= 5:
                    bed_metrics[parts[0]] = {
                        "exclusion_bp": parts[1],
                        "dip_intersect_bp": parts[2],
                        "pct_of_exclusion": parts[3],
                        "pct_of_dip": parts[4],
                    }

    # Write output
    with open(output_csv, "w", newline="") as out_f:
        writer = csv.writer(out_f)
        if is_smvar:
            writer.writerow(
                [
                    "exclusion",
                    "excl_id",
                    "exclusion_bp",
                    "dip_intersect_bp",
                    "pct_of_exclusion",
                    "pct_of_dip",
                    "total_variants",
                    "snv_count",
                    "indel_count",
                ]
            )
        else:
            writer.writerow(
                [
                    "exclusion",
                    "excl_id",
                    "exclusion_bp",
                    "dip_intersect_bp",
                    "pct_of_exclusion",
                    "pct_of_dip",
                    "total_variants",
                    "ins_count",
                    "del_count",
                ]
            )

        for excl_id, canon_name in sorted(
            excl_name_mapping.items(), key=lambda x: x[1]
        ):
            metrics = bed_metrics.get(canon_name, {})
            c = counts.get(excl_id, Counter())
            total = sum(c.values())

            if is_smvar:
                writer.writerow(
                    [
                        canon_name,
                        excl_id,
                        metrics.get("exclusion_bp", 0),
                        metrics.get("dip_intersect_bp", 0),
                        metrics.get("pct_of_exclusion", 0),
                        metrics.get("pct_of_dip", 0),
                        total,
                        c.get("SNV", 0),
                        c.get("INDEL", 0),
                    ]
                )
            else:
                writer.writerow(
                    [
                        canon_name,
                        excl_id,
                        metrics.get("exclusion_bp", 0),
                        metrics.get("dip_intersect_bp", 0),
                        metrics.get("pct_of_exclusion", 0),
                        metrics.get("pct_of_dip", 0),
                        total,
                        c.get("INS", 0),
                        c.get("DEL", 0),
                    ]
                )

    logging.info("Output written to %s", output_csv)


if __name__ == "__main__":
    benchmark = snakemake.wildcards.benchmark
    bench_type = "stvar" if "stvar" in benchmark else "smvar"
    excl_name_mapping = dict(snakemake.params.excl_name_mapping)

    count_exclusion_variants(
        parquet_path=snakemake.input.variant_table,
        exclusion_tsvs=[str(p) for p in snakemake.input.exclusion_tsvs],
        output_csv=snakemake.output.csv,
        bench_type=bench_type,
        excl_name_mapping=excl_name_mapping,
        log_path=snakemake.log[0],
    )
