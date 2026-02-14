#!/usr/bin/env python3
# AI Disclosure: This script was developed with assistance from Claude (Anthropic)
# for code generation and debugging.
"""
Compute exclusion interaction (upset-style) decomposition.

For each unique combination of overlapping exclusions, computes:
- bases_bp: bases in dip.bed excluded by exactly this combination
- variant_count: size-filtered variants excluded by exactly this combination

Reads the Parquet variant table (with boolean excl_* columns).
"""

import csv
import os
import subprocess
import tempfile
from collections import Counter
from pathlib import Path
from typing import Dict

import pandas as pd


def compute_base_interactions(
    dip_bed: str,
    benchmark_bed: str,
    exclusion_beds: list[str],
    exclusion_names: list[str],
    log,
) -> Dict[str, int]:
    """Compute base-level exclusion interactions using bedtools multiinter.

    Returns dict mapping pipe-delimited exclusion combination to base count.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        # Compute excluded regions (dip.bed minus benchmark.bed)
        excluded_bed = os.path.join(tmpdir, "excluded.bed")
        cmd = f"bedtools subtract -a {dip_bed} -b {benchmark_bed} > {excluded_bed}"
        result = subprocess.run(
            cmd, shell=True, capture_output=True, text=True, executable="/bin/bash"
        )
        if result.returncode != 0:
            error_msg = f"bedtools subtract failed: {result.stderr}"
            log.write(f"ERROR computing excluded regions: {error_msg}\n")
            raise RuntimeError(error_msg)

        # Run bedtools multiinter on all exclusion BEDs
        beds_str = " ".join(exclusion_beds)
        names_str = " ".join(exclusion_names)
        multiinter_out = os.path.join(tmpdir, "multiinter.bed")

        cmd = (
            f"bedtools multiinter -i {beds_str} -names {names_str} "
            f"| bedtools intersect -a - -b {excluded_bed} "
            f"> {multiinter_out}"
        )
        result = subprocess.run(
            cmd, shell=True, capture_output=True, text=True, executable="/bin/bash"
        )
        if result.returncode != 0:
            error_msg = (
                f"bedtools multiinter/intersect pipeline failed: {result.stderr}"
            )
            log.write(f"ERROR in multiinter: {error_msg}\n")
            raise RuntimeError(error_msg)

        # Parse multiinter output
        # Format: chrom start end num_overlaps names_list [per-bed columns]
        # The 5th column onward contains the exclusion names
        combo_bases: Dict[str, int] = Counter()
        with open(multiinter_out, "r") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 5:
                    continue
                _chrom, start, end = parts[0], int(parts[1]), int(parts[2])
                bp = end - start
                # Column 4 is num_overlaps, column 5 is comma-separated names
                names = sorted(parts[4].split(","))
                combo_key = "|".join(names)
                combo_bases[combo_key] += bp

        log.write(f"Found {len(combo_bases)} base-level combinations\n")
        return dict(combo_bases)


def compute_variant_interactions(
    parquet_path: str,
    excl_column_mapping: Dict[str, str],
    log,
) -> Dict[str, int]:
    """Compute variant-level exclusion interactions from Parquet variant table.

    Uses boolean excl_* columns to determine which exclusions overlap each
    variant, then counts unique combinations.

    Returns dict mapping pipe-delimited exclusion combination to variant count.
    """
    combo_counts: Dict[str, int] = Counter()

    excl_cols = list(excl_column_mapping.keys())
    df = pd.read_parquet(parquet_path, columns=excl_cols)
    log.write(f"Loaded {len(df)} variants from Parquet\n")
    log.write(f"Exclusion columns: {excl_cols}\n")

    # For each variant, find which exclusion columns are True
    for idx, row in df.iterrows():
        active_excls = sorted([
            excl_column_mapping[col]
            for col in excl_cols
            if col in df.columns and row.get(col, False)
        ])

        if not active_excls:
            continue

        combo_key = "|".join(active_excls)
        combo_counts[combo_key] += 1

    log.write(f"Found {len(combo_counts)} variant-level combinations\n")
    return dict(combo_counts)


def main():
    benchmark = snakemake.wildcards.benchmark
    excl_column_mapping = dict(snakemake.params.excl_column_mapping)

    dip_bed = str(snakemake.input.dip_bed)
    benchmark_bed = str(snakemake.input.benchmark_bed)
    exclusion_beds = list(snakemake.input.exclusion_beds)
    parquet_path = str(snakemake.input.variant_table)
    output_csv = str(snakemake.output.csv)
    log_path = str(snakemake.log[0])

    # Derive canonical names from bed paths (filename without .bed extension)
    exclusion_names = [Path(p).stem for p in exclusion_beds]

    with open(log_path, "w") as log:
        log.write(f"Computing exclusion interactions for {benchmark}\n")
        log.write(f"Exclusion BEDs: {exclusion_beds}\n")
        log.write(f"Exclusion column mapping: {excl_column_mapping}\n\n")

        # Base-level interactions
        log.write("Computing base-level interactions...\n")
        base_combos = compute_base_interactions(
            dip_bed, benchmark_bed, exclusion_beds, exclusion_names, log
        )

        # Variant-level interactions
        log.write("Computing variant-level interactions...\n")
        var_combos = compute_variant_interactions(parquet_path, excl_column_mapping, log)

        # Merge and write output
        all_combos = set(base_combos.keys()) | set(var_combos.keys())

        with open(output_csv, "w", newline="") as out_f:
            writer = csv.writer(out_f)
            writer.writerow(
                [
                    "exclusion_combination",
                    "n_exclusions",
                    "bases_bp",
                    "variant_count",
                ]
            )

            for combo in sorted(all_combos):
                n_excl = len(combo.split("|"))
                writer.writerow(
                    [
                        combo,
                        n_excl,
                        base_combos.get(combo, 0),
                        var_combos.get(combo, 0),
                    ]
                )

        log.write(f"\nOutput: {len(all_combos)} combinations written\n")


if __name__ == "__main__":
    main()
