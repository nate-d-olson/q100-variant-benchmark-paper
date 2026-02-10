#!/usr/bin/env python3
# AI Disclosure: This script was developed with assistance from Claude (Anthropic)
# for code generation and debugging.
"""
Count variants per exclusion from variant table REGION_IDS annotations.

Reads the variant table TSV, parses REGION_IDS for EXCL_* tags, and counts
variants per exclusion with size filtering (smvar <50bp, stvar >=50bp).
Merges with per-exclusion BED metrics to produce a combined impact table.
"""

import csv
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
    """Compute variant size in bp.

    For SVs with SVLEN, use abs(SVLEN).
    Otherwise use max(len(REF), len(ALT)) - 1.
    """
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


def _classify_sv_type(
    row: dict, svtype_col: str | None, ref_col: str, alt_col: str
) -> str:
    """Classify SV as insertion or deletion."""
    if svtype_col:
        svtype = row.get(svtype_col, ".")
        if svtype and svtype != ".":
            svtype_upper = svtype.upper()
            if "INS" in svtype_upper:
                return "INS"
            if "DEL" in svtype_upper:
                return "DEL"
    ref_len = int(row.get(ref_col, "1"))
    alt_len = int(row.get(alt_col, "1"))
    return "INS" if alt_len > ref_len else "DEL"


def count_exclusion_variants(
    variant_table: Path,
    exclusion_tsvs: list[Path],
    output_csv: Path,
    bench_type: str,
    excl_name_mapping: Dict[str, str],
    log_path: Path,
) -> None:
    """Count variants per exclusion with size filtering.

    Args:
        variant_table: Path to variants.tsv
        exclusion_tsvs: Paths to per-exclusion BED metric TSV files
        output_csv: Path to output CSV
        bench_type: "smvar" or "stvar"
        excl_name_mapping: Maps EXCL_* IDs to canonical names
        log_path: Path to log file
    """
    is_smvar = bench_type == "smvar"
    # Counters: excl_id -> {total, snp, indel, sv_ins, sv_del}
    counts: Dict[str, Counter] = {}
    for excl_id in excl_name_mapping:
        counts[excl_id] = Counter()

    total_variants = 0
    filtered_variants = 0
    excluded_variants = 0

    with open(log_path, "w") as log:
        log.write(f"Counting exclusion variants for {bench_type}\n")
        log.write(f"Input: {variant_table}\n")
        log.write(f"Exclusion mapping: {excl_name_mapping}\n\n")

        with open(variant_table, "r") as f:
            reader = csv.DictReader(f, delimiter="\t")
            fieldnames = reader.fieldnames or []

            region_col = _find_col(fieldnames, {"INFO/REGION_IDS", "REGION_IDS"})
            type_col = _find_col(fieldnames, {"TYPE"})
            ref_col = _find_col(fieldnames, {"STRLEN(REF)"})
            alt_col = _find_col(fieldnames, {"STRLEN(ALT)"})
            svlen_col = _find_col(fieldnames, {"INFO/SVLEN", "SVLEN"})
            svtype_col = _find_col(fieldnames, {"INFO/SVTYPE", "SVTYPE"})

            if not region_col or not type_col or not ref_col or not alt_col:
                raise ValueError(f"Missing required columns. Found: {fieldnames}")

            log.write(f"Using columns: REGION_IDS={region_col}, TYPE={type_col}\n")
            log.write(f"SVLEN={svlen_col}, SVTYPE={svtype_col}\n\n")

            for row in reader:
                total_variants += 1
                region_ids_str = row.get(region_col, "")
                if not region_ids_str or region_ids_str == ".":
                    continue

                # Compute variant size and apply filter
                var_size = _compute_variant_size(row, ref_col, alt_col, svlen_col)
                if is_smvar and var_size >= 50:
                    continue
                if not is_smvar and var_size < 50:
                    continue
                filtered_variants += 1

                region_ids = [r.strip() for r in region_ids_str.split(",")]
                excl_ids = [r for r in region_ids if r.startswith("EXCL_")]

                if not excl_ids:
                    continue
                excluded_variants += 1

                var_type = row.get(type_col, "UNKNOWN")

                for excl_id in excl_ids:
                    if excl_id not in counts:
                        counts[excl_id] = Counter()
                    counts[excl_id]["total"] += 1
                    if is_smvar:
                        if var_type == "SNP":
                            counts[excl_id]["snp"] += 1
                        else:
                            counts[excl_id]["indel"] += 1
                    else:
                        sv_type = _classify_sv_type(row, svtype_col, ref_col, alt_col)
                        if sv_type == "INS":
                            counts[excl_id]["sv_ins"] += 1
                        else:
                            counts[excl_id]["sv_del"] += 1

                if total_variants % 100000 == 0:
                    log.write(f"Processed {total_variants:,} variants...\n")

        log.write(f"\nTotal variants: {total_variants:,}\n")
        log.write(f"Size-filtered variants: {filtered_variants:,}\n")
        log.write(f"Excluded variants: {excluded_variants:,}\n\n")

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
                        "snp_count",
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
                        "sv_ins_count",
                        "sv_del_count",
                    ]
                )

            for excl_id, canon_name in sorted(
                excl_name_mapping.items(), key=lambda x: x[1]
            ):
                metrics = bed_metrics.get(canon_name, {})
                c = counts.get(excl_id, Counter())
                if is_smvar:
                    writer.writerow(
                        [
                            canon_name,
                            excl_id,
                            metrics.get("exclusion_bp", 0),
                            metrics.get("dip_intersect_bp", 0),
                            metrics.get("pct_of_exclusion", 0),
                            metrics.get("pct_of_dip", 0),
                            c["total"],
                            c["snp"],
                            c["indel"],
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
                            c["total"],
                            c["sv_ins"],
                            c["sv_del"],
                        ]
                    )

        log.write(f"Output written to {output_csv}\n")


if __name__ == "__main__":
    # Parse bench_type from benchmark wildcard
    benchmark = snakemake.wildcards.benchmark
    bench_type = "stvar" if "stvar" in benchmark else "smvar"

    # Build excl_name_mapping from params
    excl_name_mapping = dict(snakemake.params.excl_name_mapping)

    count_exclusion_variants(
        variant_table=Path(snakemake.input.variant_table),
        exclusion_tsvs=[Path(p) for p in snakemake.input.exclusion_tsvs],
        output_csv=Path(snakemake.output.csv),
        bench_type=bench_type,
        excl_name_mapping=excl_name_mapping,
        log_path=Path(snakemake.log[0]),
    )
