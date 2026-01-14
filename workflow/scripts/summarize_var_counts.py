#!/usr/bin/env python3
"""
Summarize variant counts by stratification.

Aggregates variant type counts into summary statistics per stratification.
"""

import csv
from pathlib import Path
from collections import defaultdict
from typing import Dict


def summarize_variant_counts(
    input_csv: Path,
    output_csv: Path,
    log_path: Path
) -> None:
    """
    Summarize variant counts by stratification.

    Args:
        input_csv: Path to variants_by_stratification.csv
        output_csv: Path to output summary CSV
        log_path: Path to log file

    Output format:
        strat_name,total_variants,snp_count,indel_count,other_count
    """
    # {strat_name: {var_type: count}}
    strat_summary: Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))

    with open(log_path, 'w') as log:
        log.write(f"Summarizing variant counts\n")
        log.write(f"Input: {input_csv}\n")
        log.write(f"Output: {output_csv}\n\n")

        try:
            # Read input counts
            with open(input_csv, 'r') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    strat_name = row['strat_name']
                    var_type = row['var_type']
                    count = int(row['count'])

                    strat_summary[strat_name][var_type] += count

            log.write(f"Processed {len(strat_summary)} stratifications\n\n")

            # Write summary
            with open(output_csv, 'w', newline='') as out_f:
                writer = csv.writer(out_f)
                writer.writerow([
                    'strat_name',
                    'total_variants',
                    'snp_count',
                    'indel_count',
                    'del_count',
                    'ins_count',
                    'complex_count',
                    'other_count'
                ])

                for strat_name in sorted(strat_summary.keys()):
                    var_counts = strat_summary[strat_name]

                    # Calculate totals
                    total = sum(var_counts.values())
                    snp = var_counts.get('SNP', 0)
                    indel = var_counts.get('INDEL', 0)
                    dels = var_counts.get('DEL', 0)
                    ins = var_counts.get('INS', 0)
                    complex_var = var_counts.get('COMPLEX', 0)

                    # "Other" is everything not SNP/INDEL/DEL/INS/COMPLEX
                    known_types = {'SNP', 'INDEL', 'DEL', 'INS', 'COMPLEX'}
                    other = sum(count for vtype, count in var_counts.items()
                               if vtype not in known_types)

                    writer.writerow([
                        strat_name,
                        total,
                        snp,
                        indel,
                        dels,
                        ins,
                        complex_var,
                        other
                    ])

                    log.write(f"{strat_name}: {total:,} total variants\n")

            log.write(f"\nSummary written to {output_csv}\n")

        except Exception as e:
            log.write(f"ERROR: {str(e)}\n")
            raise


if __name__ == '__main__':
    # Snakemake provides these variables
    input_file = Path(snakemake.input.csv)
    output_file = Path(snakemake.output.summary)
    log_file = Path(snakemake.log[0])

    summarize_variant_counts(input_file, output_file, log_file)
