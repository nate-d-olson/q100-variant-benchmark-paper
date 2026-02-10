#!/usr/bin/env python3
"""
Count variants by genomic context from variant table.

Reads variant table TSV and counts variants in each genomic context.
Handles comma-separated CONTEXT_IDS field where variants can overlap multiple regions.
"""

import csv
from pathlib import Path
from typing import Counter, Dict


def count_variants_by_genomic_context(
    variant_table_path: Path, output_path: Path, log_path: Path
) -> None:
    """
    Count variants by genomic context region and variant type.

    Args:
        variant_table_path: Path to variants.tsv file
        output_path: Path to output CSV file
        log_path: Path to log file

    Output format:
        context_name,var_type,count
    """
    # Counters: {(context_name, var_type): count}
    counts: Dict[tuple[str, str], int] = Counter()
    total_variants = 0
    variants_with_context = 0
    line_num = 0

    with open(log_path, "w") as log:
        log.write("Counting variants by genomic context\n")
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
                context_col = find_col({"INFO/CONTEXT_IDS", "CONTEXT_IDS"})

                # Verify required columns exist
                missing = []
                if not chrom_col:
                    missing.append("CHROM")
                if not type_col:
                    missing.append("TYPE")
                if not context_col:
                    missing.append("CONTEXT_IDS")

                if missing:
                    raise ValueError(f"Missing required columns: {missing}")

                log.write(
                    f"Using columns: TYPE='{type_col}', CONTEXT_IDS='{context_col}'\n\n"
                )

                for row in reader:
                    line_num += 1
                    total_variants += 1

                    var_type = row.get(type_col, "UNKNOWN")
                    context_ids = row.get(context_col, "")

                    # Handle missing or empty CONTEXT_IDS
                    if not context_ids or context_ids == "." or context_ids == "":
                        # Variant not in any genomic context
                        continue

                    variants_with_context += 1

                    # Split comma-separated genomic context IDs
                    for context_name in context_ids.split(","):
                        context_name = context_name.strip()
                        if context_name:
                            counts[(context_name, var_type)] += 1

                    if line_num % 100000 == 0:
                        log.write(f"Processed {line_num:,} variants...\n")

        except Exception as e:
            log.write(f"ERROR: {str(e)}\n")
            raise

        log.write("\nProcessing complete:\n")
        log.write(f"  Total variants: {total_variants:,}\n")
        log.write(f"  Variants with genomic context: {variants_with_context:,}\n")
        log.write(f"  Unique (context, var_type) combinations: {len(counts)}\n\n")

        # Write output CSV
        with open(output_path, "w", newline="") as out_f:
            writer = csv.writer(out_f)
            writer.writerow(["context_name", "var_type", "count"])

            # Sort by genomic context name, then variant type
            for (context_name, var_type), count in sorted(counts.items()):
                writer.writerow([context_name, var_type, count])

        log.write(f"Output written to {output_path}\n")
        log.write(f"Output lines: {len(counts) + 1} (including header)\n")


if __name__ == "__main__":
    # Snakemake provides these variables
    variant_table = Path(snakemake.input.tsv)
    output_csv = Path(snakemake.output.csv)
    log_file = Path(snakemake.log[0])

    count_variants_by_genomic_context(variant_table, output_csv, log_file)
