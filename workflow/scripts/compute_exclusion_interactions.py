#!/usr/bin/env python3
# AI Disclosure: This script was developed with assistance from Claude (Anthropic)
# for code generation and debugging.
"""
Compute exclusion interaction (upset-style) decomposition.

For each unique combination of overlapping exclusions, computes:
- bases_bp: bases in dip.bed excluded by exactly this combination
- variant_count: size-filtered variants excluded by exactly this combination
"""

import csv
import os
import subprocess
import tempfile
from collections import Counter
from pathlib import Path
from typing import Dict


def _find_col(fieldnames: list[str], candidates: set[str]) -> str | None:
    """Find column name handling bcftools [n] prefix."""
    for f in fieldnames:
        norm = f.split("]", 1)[1] if "]" in f else f
        if norm in candidates:
            return f
    return None


def _compute_variant_size(
    row: dict, ref_col: str, alt_col: str, svlen_col: str | None
) -> int:
    """Compute variant size in bp."""
    if svlen_col:
        svlen_val = row.get(svlen_col, ".")
        if svlen_val and svlen_val != ".":
            try:
                return abs(int(svlen_val))
            except ValueError:
                pass
    ref_len = int(row.get(ref_col, "1"))
    alt_len = int(row.get(alt_col, "1"))
    return max(ref_len, alt_len) - 1


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
            log.write(f"ERROR computing excluded regions: {result.stderr}\n")
            return {}

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
            log.write(f"ERROR in multiinter: {result.stderr}\n")
            return {}

        # Parse multiinter output
        # Format: chrom start end num_overlaps names_list [per-bed columns]
        # The 5th column onward contains the exclusion names
        combo_bases: Dict[str, int] = Counter()
        with open(multiinter_out, "r") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 5:
                    continue
                chrom, start, end = parts[0], int(parts[1]), int(parts[2])
                bp = end - start
                # Column 4 is num_overlaps, column 5 is comma-separated names
                names = sorted(parts[4].split(","))
                combo_key = "|".join(names)
                combo_bases[combo_key] += bp

        log.write(f"Found {len(combo_bases)} base-level combinations\n")
        return dict(combo_bases)


def compute_variant_interactions(
    variant_table: str,
    bench_type: str,
    excl_name_mapping: Dict[str, str],
    log,
) -> Dict[str, int]:
    """Compute variant-level exclusion interactions from variant table.

    Returns dict mapping pipe-delimited exclusion combination to variant count.
    """
    is_smvar = bench_type == "smvar"
    combo_counts: Dict[str, int] = Counter()

    with open(variant_table, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        fieldnames = reader.fieldnames or []

        region_col = _find_col(fieldnames, {"INFO/REGION_IDS", "REGION_IDS"})
        ref_col = _find_col(fieldnames, {"STRLEN(REF)"})
        alt_col = _find_col(fieldnames, {"STRLEN(ALT)"})
        svlen_col = _find_col(fieldnames, {"INFO/SVLEN", "SVLEN"})

        missing_cols = []
        if not region_col:
            missing_cols.append("INFO/REGION_IDS or REGION_IDS")
        if not ref_col:
            missing_cols.append("STRLEN(REF)")
        if not alt_col:
            missing_cols.append("STRLEN(ALT)")

        if missing_cols:
            msg = (
                f"Missing required columns in variant table {variant_table}: "
                + ", ".join(missing_cols)
            )
            log.write(f"ERROR: {msg}\n")
            raise ValueError(msg)
        for row in reader:
            region_ids_str = row.get(region_col, "")
            if not region_ids_str or region_ids_str == ".":
                continue

            var_size = _compute_variant_size(row, ref_col, alt_col, svlen_col)
            if is_smvar and var_size >= 50:
                continue
            if not is_smvar and var_size < 50:
                continue

            region_ids = [r.strip() for r in region_ids_str.split(",")]
            excl_ids = sorted([r for r in region_ids if r.startswith("EXCL_")])

            if not excl_ids:
                continue

            # Map EXCL_* IDs to canonical names
            canon_names = sorted([excl_name_mapping.get(eid, eid) for eid in excl_ids])
            combo_key = "|".join(canon_names)
            combo_counts[combo_key] += 1

    log.write(f"Found {len(combo_counts)} variant-level combinations\n")
    return dict(combo_counts)


def main():
    benchmark = snakemake.wildcards.benchmark
    bench_type = "stvar" if "stvar" in benchmark else "smvar"
    excl_name_mapping = dict(snakemake.params.excl_name_mapping)

    dip_bed = str(snakemake.input.dip_bed)
    benchmark_bed = str(snakemake.input.benchmark_bed)
    exclusion_beds = list(snakemake.input.exclusion_beds)
    variant_table = str(snakemake.input.variant_table)
    output_csv = str(snakemake.output.csv)
    log_path = str(snakemake.log[0])

    # Derive canonical names from bed paths (filename without .bed extension)
    exclusion_names = [Path(p).stem for p in exclusion_beds]

    with open(log_path, "w") as log:
        log.write(f"Computing exclusion interactions for {benchmark}\n")
        log.write(f"Bench type: {bench_type}\n")
        log.write(f"Exclusion BEDs: {exclusion_beds}\n\n")

        # Base-level interactions
        log.write("Computing base-level interactions...\n")
        base_combos = compute_base_interactions(
            dip_bed, benchmark_bed, exclusion_beds, exclusion_names, log
        )

        # Variant-level interactions
        log.write("Computing variant-level interactions...\n")
        var_combos = compute_variant_interactions(
            variant_table, bench_type, excl_name_mapping, log
        )

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
