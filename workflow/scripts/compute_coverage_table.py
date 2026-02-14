#!/usr/bin/env python3
# AI Disclosure: This script was developed with assistance from Claude (Anthropic).
"""
Compute genomic context coverage metrics for all contexts in a benchmark.

For each genomic context BED file, computes:
- context_name: Name of the genomic context region
- context_bp: Total bases in the context (after merging overlaps)
- intersect_bp: Bases overlapping between context and benchmark regions
- pct_of_context: Percent of context covered by benchmark
- pct_of_bench: Percent of benchmark covered by context

Output: Single CSV file with all context metrics.
"""

import csv

from compute_bed_metrics import compute_bed_size, compute_intersection
from logging_config import setup_logger

logger = setup_logger(__name__)


def main():
    """Main entry point for Snakemake script execution."""
    context_beds = snakemake.input.context_beds
    bench_bed = snakemake.input.bench_bed
    output_csv = snakemake.output.csv
    context_names = snakemake.params.context_names

    logger.info(f"Computing coverage table for {len(context_beds)} contexts")
    logger.info(f"Benchmark BED: {bench_bed}")

    bench_size = compute_bed_size(bench_bed, bench_bed.endswith(".gz"))
    logger.info(f"Benchmark size: {bench_size} bp")

    rows = []
    for context_name, context_bed in zip(context_names, context_beds):
        is_gzipped = context_bed.endswith(".gz")

        context_size = compute_bed_size(context_bed, is_gzipped)
        intersect = compute_intersection(
            context_bed,
            bench_bed,
            a_gzipped=is_gzipped,
            b_gzipped=bench_bed.endswith(".gz"),
        )

        pct_of_context = (intersect / context_size * 100) if context_size > 0 else 0.0
        pct_of_bench = (intersect / bench_size * 100) if bench_size > 0 else 0.0

        rows.append(
            {
                "context_name": context_name,
                "context_bp": context_size,
                "intersect_bp": intersect,
                "pct_of_context": f"{pct_of_context:.6f}",
                "pct_of_bench": f"{pct_of_bench:.6f}",
            }
        )

        logger.info(
            f"  {context_name}: size={context_size}, intersect={intersect}, "
            f"pct_of_context={pct_of_context:.2f}%, pct_of_bench={pct_of_bench:.2f}%"
        )

    with open(output_csv, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "context_name",
                "context_bp",
                "intersect_bp",
                "pct_of_context",
                "pct_of_bench",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    logger.info(f"Wrote {len(rows)} rows to {output_csv}")


if __name__ == "__main__":
    main()
