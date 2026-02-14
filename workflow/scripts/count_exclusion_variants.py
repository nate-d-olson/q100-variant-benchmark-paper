#!/usr/bin/env python3
# AI Disclosure: This script was developed with assistance from Claude (Anthropic).
"""
Count variants per exclusion from Parquet variant table.

Reads the Parquet variant table (with boolean excl_* columns), counts variants
per exclusion, and merges with per-exclusion BED metrics to produce a combined
impact table.

The Parquet table has boolean columns named excl_<name> for each exclusion
(e.g., excl_flanks, excl_satellites). This script reads those columns directly
instead of parsing comma-separated region_ids strings.
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
    excl_column_mapping: Dict[str, str],
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
        excl_column_mapping: Maps excl_* column names to canonical names
            (e.g., {"excl_flanks": "flanks", "excl_satellites": "satellites"})
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
    logging.info("Exclusion column mapping: %s", excl_column_mapping)

    # Read exclusion boolean columns + var_type from Parquet
    excl_cols = list(excl_column_mapping.keys())
    columns_to_read = ["var_type"] + excl_cols
    logging.info("Reading columns: %s", columns_to_read)

    df = pd.read_parquet(parquet_path, columns=columns_to_read)
    logging.info("Loaded %d variants", len(df))

    # Count variants per exclusion column, broken down by var_type
    counts: Dict[str, Counter] = {col: Counter() for col in excl_cols}
    for excl_col in excl_cols:
        if excl_col not in df.columns:
            logging.warning("Exclusion column %s not found in Parquet", excl_col)
            continue
        # Filter to variants in this exclusion
        in_excl = df[df[excl_col]]
        for var_type, count in in_excl["var_type"].value_counts().items():
            counts[excl_col][var_type] = count
        logging.info(
            "  %s: %d variants (%s)",
            excl_col,
            len(in_excl),
            dict(counts[excl_col]),
        )

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
                    "excl_column",
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
                    "excl_column",
                    "exclusion_bp",
                    "dip_intersect_bp",
                    "pct_of_exclusion",
                    "pct_of_dip",
                    "total_variants",
                    "ins_count",
                    "del_count",
                ]
            )

        for excl_col, canon_name in sorted(
            excl_column_mapping.items(), key=lambda x: x[1]
        ):
            metrics = bed_metrics.get(canon_name, {})
            c = counts.get(excl_col, Counter())
            total = sum(c.values())

            if is_smvar:
                writer.writerow(
                    [
                        canon_name,
                        excl_col,
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
                        excl_col,
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
    excl_column_mapping = dict(snakemake.params.excl_column_mapping)

    count_exclusion_variants(
        parquet_path=snakemake.input.variant_table,
        exclusion_tsvs=[str(p) for p in snakemake.input.exclusion_tsvs],
        output_csv=snakemake.output.csv,
        bench_type=bench_type,
        excl_column_mapping=excl_column_mapping,
        log_path=snakemake.log[0],
    )
