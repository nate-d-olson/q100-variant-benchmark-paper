#!/usr/bin/env python3
"""
Count variants by stratification region from variant table.

Reads variant table TSV and counts variants in each stratification region.
Handles comma-separated STRAT_IDS field where variants can overlap multiple regions.
"""

from pathlib import Path
from typing import Dict, Counter
import csv


def count_variants_by_stratification(
    variant_table_path: Path, output_path: Path, log_path: Path
) -> None:
    """
    Count variants by stratification region and variant type.

    Args:
        variant_table_path: Path to variants.tsv file
        output_path: Path to output CSV file
        log_path: Path to log file

    Output format:
        strat_name,var_type,count
    """
    # Counters: {(strat_name, var_type): count}
    counts: Dict[tuple[str, str], int] = Counter()
    total_variants = 0
    variants_with_strat = 0
    line_num = 0

    with open(log_path, "w") as log:
        log.write(f"Counting variants by stratification\n")
        log.write(f"Input: {variant_table_path}\n")
        log.write(f"Output: {output_path}\n\n")

        try:
            with open(variant_table_path, "r") as f:
                reader = csv.DictReader(f, delimiter="\t")
                fieldnames = reader.fieldnames or []

                log.write(f"Found columns: {fieldnames}\n\n")

                # Helper to find column by normalized name
                def find_col(candidates):
                    for f in fieldnames:
                        # Handle bcftools [n] prefix (e.g., "[1]CHROM")
                        norm = f.split("]", 1)[1] if "]" in f else f
                        if norm in candidates:
                            return f
                    return None

                # Resolve columns
                chrom_col = find_col({"#CHROM", "CHROM"})
                type_col = find_col({"TYPE"})
                strat_col = find_col({"INFO/STRAT_IDS", "STRAT_IDS"})

                # Verify required columns exist
                missing = []
                if not chrom_col:
                    missing.append("CHROM")
                if not type_col:
                    missing.append("TYPE")
                if not strat_col:
                    missing.append("STRAT_IDS")

                if missing:
                    raise ValueError(f"Missing required columns: {missing}")

                log.write(
                    f"Using columns: TYPE='{type_col}', STRAT_IDS='{strat_col}'\n\n"
                )

                for row in reader:
                    line_num += 1
                    total_variants += 1

                    var_type = row.get(type_col, "UNKNOWN")
                    strat_ids = row.get(strat_col, "")

                    # Handle missing or empty STRAT_IDS
                    if not strat_ids or strat_ids == "." or strat_ids == "":
                        # Variant not in any stratification
                        continue

                    variants_with_strat += 1

                    # Split comma-separated stratification IDs
                    for strat_name in strat_ids.split(","):
                        strat_name = strat_name.strip()
                        if strat_name:
                            counts[(strat_name, var_type)] += 1

                    if line_num % 100000 == 0:
                        log.write(f"Processed {line_num:,} variants...\n")

        except Exception as e:
            log.write(f"ERROR: {str(e)}\n")
            raise

        log.write(f"\nProcessing complete:\n")
        log.write(f"  Total variants: {total_variants:,}\n")
        log.write(f"  Variants with stratification: {variants_with_strat:,}\n")
        log.write(f"  Unique (strat, var_type) combinations: {len(counts)}\n\n")

        # Write output CSV
        with open(output_path, "w", newline="") as out_f:
            writer = csv.writer(out_f)
            writer.writerow(["strat_name", "var_type", "count"])

            # Sort by stratification name, then variant type
            for (strat_name, var_type), count in sorted(counts.items()):
                writer.writerow([strat_name, var_type, count])

        log.write(f"Output written to {output_path}\n")
        log.write(f"Output lines: {len(counts) + 1} (including header)\n")


if __name__ == "__main__":
    # Snakemake provides these variables
    variant_table = Path(snakemake.input.tsv)
    output_csv = Path(snakemake.output.csv)
    log_file = Path(snakemake.log[0])

    count_variants_by_stratification(variant_table, output_csv, log_file)
