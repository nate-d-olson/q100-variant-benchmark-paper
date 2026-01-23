#!/usr/bin/env python3
"""
Aggregate truvari stratify results into summary tables.

Processes multiple TSV files from truvari stratify and generates two summaries:
1. summary_by_variant.tsv: Variant-level metrics (TP, FP, FN, Precision, Recall, F1)
2. summary_by_region.tsv: Region-level statistics (counts, bases, average size)

Input TSV format from truvari stratify:
    chrom  start  end  tpbase  tp  fn  fp

Output format (summary_by_variant.csv):
    stratification,region_type,mode,TP_base,TP,FN,FP,Precision,Recall,F1

Output format (summary_by_region.csv):
    stratification,region_type,mode,Total_Regions,Total_Bases,Avg_Region_Size
"""

from pathlib import Path
from typing import Tuple
import sys

import pandas as pd
import truvari


def parse_filename(filepath: Path) -> Tuple[str, str, str]:
    """
    Extract metadata from TSV filename.

    Args:
        filepath: Path to TSV file (e.g., TR_regions_overlap.tsv)

    Returns:
        Tuple of (stratification, region_type, mode)

    Example:
        >>> parse_filename(Path("TR_regions_overlap.tsv"))
        ('TR', 'regions', 'overlap')
    """
    stem = filepath.stem  # Remove .tsv extension
    parts = stem.split("_")

    if len(parts) != 3:
        raise ValueError(
            f"Unexpected filename format: {filepath.name}. "
            f"Expected: {{strat}}_{{region_type}}_{{mode}}.tsv"
        )

    stratification, region_type, mode = parts
    return stratification, region_type, mode


def calculate_region_metrics(df: pd.DataFrame) -> dict:
    """
    Calculate region-level statistics from stratify dataframe.

    Args:
        df: DataFrame with chrom, start, end columns

    Returns:
        Dict with Total_Regions, Total_Bases, Avg_Region_Size
    """
    total_regions = len(df)
    total_bases = (df["end"] - df["start"]).sum()
    avg_size = total_bases / total_regions if total_regions > 0 else 0.0

    return {
        "Total_Regions": total_regions,
        "Total_Bases": total_bases,
        "Avg_Region_Size": avg_size,
    }


def main():
    """Main aggregation logic."""
    # Snakemake provides input/output as special objects
    input_tsvs = snakemake.input.tsvs  # type: ignore
    variant_summary_path = snakemake.output.variant_summary  # type: ignore
    region_summary_path = snakemake.output.region_summary  # type: ignore
    log_path = snakemake.log[0]  # type: ignore

    variant_results = []
    region_results = []

    with open(log_path, "w") as log:
        log.write(f"Processing {len(input_tsvs)} TSV files\n\n")

        for tsv_path in input_tsvs:
            tsv_path = Path(tsv_path)
            log.write(f"Processing: {tsv_path.name}\n")

            try:
                # Extract metadata from filename
                strat, region_type, mode = parse_filename(tsv_path)

                # Read TSV with pandas
                df = pd.read_csv(tsv_path, sep="\t")
                df.columns = ["chrom", "start", "end", "tpbase", "tp", "fn", "fp"]

                # Calculate totals
                tp_base_total = df["tpbase"].sum()
                tp_total = df["tp"].sum()
                fn_total = df["fn"].sum()
                fp_total = df["fp"].sum()

                # Use truvari.performance_metrics to calculate precision, recall, F1
                precision, recall, f1 = truvari.performance_metrics(
                    tp_base_total, tp_total, fn_total, fp_total
                )

                # Calculate region metrics
                region_metrics = calculate_region_metrics(df)

                # Store variant summary
                variant_results.append(
                    {
                        "stratification": strat,
                        "region_type": region_type,
                        "mode": mode,
                        "TP_base": tp_base_total,
                        "TP": tp_total,
                        "FN": fn_total,
                        "FP": fp_total,
                        "Precision": precision,
                        "Recall": recall,
                        "F1": f1,
                    }
                )

                # Store region summary
                region_results.append(
                    {
                        "stratification": strat,
                        "region_type": region_type,
                        "mode": mode,
                        "Total_Regions": region_metrics["Total_Regions"],
                        "Total_Bases": region_metrics["Total_Bases"],
                        "Avg_Region_Size": region_metrics["Avg_Region_Size"],
                    }
                )

                log.write(
                    f"  Regions: {region_metrics['Total_Regions']}, "
                    f"TP: {tp_total}, FP: {fp_total}, FN: {fn_total}\n"
                )

            except Exception as e:
                log.write(f"  ERROR: {e}\n")
                sys.exit(1)

        log.write(f"\nWriting variant summary to {variant_summary_path}\n")
        log.write(f"Writing region summary to {region_summary_path}\n")

    # Write variant summary TSV using pandas
    variant_df = pd.DataFrame(variant_results)
    variant_df.to_csv(variant_summary_path, sep="\t", index=False)

    # Write region summary TSV using pandas
    region_df = pd.DataFrame(region_results)
    region_df.to_csv(region_summary_path, sep="\t", index=False)


if __name__ == "__main__":
    main()
