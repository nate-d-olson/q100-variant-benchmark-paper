#!/usr/bin/env python3
"""
Combine stratification coverage metrics with variant counts.

Joins stratification_coverage_table.csv with stratification_summary.csv
to create a unified table for analysis and plotting.
"""

import csv
from pathlib import Path
from typing import Any, Dict


def combine_metrics_and_counts(
    metrics_csv: Path, counts_csv: Path, output_csv: Path, log_path: Path
) -> None:
    """
    Combine genomic context coverage metrics with variant counts.

    Args:
        metrics_csv: Path to genomic_context_coverage_table.csv
        counts_csv: Path to genomic_context_summary.csv
        output_csv: Path to output combined CSV
        log_path: Path to log file

    Output format:
        context_name,context_bp,intersect_bp,pct_of_context,pct_of_bench,
        total_variants,snp_count,indel_count,del_count,ins_count,
        complex_count,other_count,variant_density_per_mb
    """
    metrics: Dict[str, Dict[str, Any]] = {}
    counts: Dict[str, Dict[str, int]] = {}

    with open(log_path, "w") as log:
        log.write(f"Combining metrics and counts\n")
        log.write(f"Metrics input: {metrics_csv}\n")
        log.write(f"Counts input: {counts_csv}\n")
        log.write(f"Output: {output_csv}\n\n")

        try:
            # Read metrics
            with open(metrics_csv, "r") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    context_name = row["context_name"]
                    metrics[context_name] = {
                        "context_bp": int(row["context_bp"]),
                        "intersect_bp": int(row["intersect_bp"]),
                        "pct_of_context": float(row["pct_of_context"]),
                        "pct_of_bench": float(row["pct_of_bench"]),
                    }

            log.write(f"Loaded metrics for {len(metrics)} genomic contexts\n")

            # Read counts
            with open(counts_csv, "r") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    context_name = row["context_name"]
                    counts[context_name] = {
                        "total_variants": int(row["total_variants"]),
                        "snp_count": int(row["snp_count"]),
                        "indel_count": int(row["indel_count"]),
                        "del_count": int(row["del_count"]),
                        "ins_count": int(row["ins_count"]),
                        "complex_count": int(row["complex_count"]),
                        "other_count": int(row["other_count"]),
                    }

            log.write(f"Loaded counts for {len(counts)} genomic contexts\n\n")

            # Combine data
            all_contexts = set(metrics.keys()) | set(counts.keys())
            log.write(f"Total unique genomic contexts: {len(all_contexts)}\n")

            # Check for mismatches
            metrics_only = set(metrics.keys()) - set(counts.keys())
            counts_only = set(counts.keys()) - set(metrics.keys())

            if metrics_only:
                log.write(
                    f"WARNING: Genomic contexts in metrics but not counts: {metrics_only}\n"
                )
            if counts_only:
                log.write(
                    f"WARNING: Genomic contexts in counts but not metrics: {counts_only}\n"
                )

            log.write("\n")

            # Write combined output
            with open(output_csv, "w", newline="") as out_f:
                writer = csv.writer(out_f)
                writer.writerow(
                    [
                        "context_name",
                        "context_bp",
                        "intersect_bp",
                        "pct_of_context",
                        "pct_of_bench",
                        "total_variants",
                        "snp_count",
                        "indel_count",
                        "del_count",
                        "ins_count",
                        "complex_count",
                        "other_count",
                        "variant_density_per_mb",
                    ]
                )

                for context_name in sorted(all_contexts):
                    # Get metrics (default to 0 if missing)
                    m = metrics.get(
                        context_name,
                        {
                            "context_bp": 0,
                            "intersect_bp": 0,
                            "pct_of_context": 0.0,
                            "pct_of_bench": 0.0,
                        },
                    )

                    # Get counts (default to 0 if missing)
                    c = counts.get(
                        context_name,
                        {
                            "total_variants": 0,
                            "snp_count": 0,
                            "indel_count": 0,
                            "del_count": 0,
                            "ins_count": 0,
                            "complex_count": 0,
                            "other_count": 0,
                        },
                    )

                    # Calculate variant density (variants per Mb of benchmark-covered genomic context)
                    intersect_mb = m["intersect_bp"] / 1_000_000
                    if intersect_mb > 0:
                        density = c["total_variants"] / intersect_mb
                    else:
                        density = 0.0

                    writer.writerow(
                        [
                            context_name,
                            m["context_bp"],
                            m["intersect_bp"],
                            f"{m['pct_of_context']:.6f}",
                            f"{m['pct_of_bench']:.6f}",
                            c["total_variants"],
                            c["snp_count"],
                            c["indel_count"],
                            c["del_count"],
                            c["ins_count"],
                            c["complex_count"],
                            c["other_count"],
                            f"{density:.2f}",
                        ]
                    )

            log.write(f"Combined output written to {output_csv}\n")
            log.write(f"Output rows: {len(all_contexts)}\n")

        except Exception as e:
            log.write(f"ERROR: {str(e)}\n")
            raise


if __name__ == "__main__":
    # Snakemake provides these variables
    metrics_file = Path(snakemake.input.metrics)
    counts_file = Path(snakemake.input.counts)
    output_file = Path(snakemake.output.combined)
    log_file = Path(snakemake.log[0])

    combine_metrics_and_counts(metrics_file, counts_file, output_file, log_file)
