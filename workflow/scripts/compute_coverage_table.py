#!/usr/bin/env python3
# AI Disclosure: This script was developed with assistance from Claude (Anthropic).
"""
Summarize bedtools coverage output into a genomic context coverage table.

Reads _cov.bed files (from bedtools coverage) and a benchmark BED to produce
a CSV with per-context overlap metrics:
- context_name: Genomic context identifier (HP, TR, SD, MAP, etc.)
- context_bp: Total bases in the genomic context
- intersect_bp: Bases of context covered by benchmark regions
- pct_of_context: Percent of context covered by benchmark
- pct_of_bench: Percent of benchmark covered by context
"""

import csv

from logging_config import setup_logger

logger = setup_logger(__name__)


def compute_bed_total(bed_path: str) -> int:
    """Sum interval lengths from a BED file (end - start per row)."""
    total = 0
    with open(bed_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split("\t")
            total += int(fields[2]) - int(fields[1])
    return total


def summarize_coverage_bed(cov_bed_path: str) -> tuple[int, int]:
    """
    Summarize a bedtools coverage output file.

    Bedtools coverage columns: chrom, start, end, n_overlap, bases_cov, ivl_len, frac_cov

    Returns:
        Tuple of (context_bp, intersect_bp)
    """
    context_bp = 0
    intersect_bp = 0
    with open(cov_bed_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split("\t")
            intersect_bp += int(fields[4])
            context_bp += int(fields[5])
    return context_bp, intersect_bp


def main() -> None:
    """Main entry point for Snakemake script execution."""
    cov_beds = snakemake.input.cov_beds
    bench_bed = snakemake.input.bench_bed
    output_csv = snakemake.output.csv
    context_names = snakemake.params.context_names

    logger.info(f"Summarizing {len(cov_beds)} coverage files")

    bench_size = compute_bed_total(bench_bed)
    logger.info(f"Benchmark size: {bench_size} bp")

    rows = []
    for context_name, cov_bed in zip(context_names, cov_beds, strict=True):
        context_bp, intersect_bp = summarize_coverage_bed(cov_bed)

        pct_of_context = (intersect_bp / context_bp * 100) if context_bp > 0 else 0.0
        pct_of_bench = (intersect_bp / bench_size * 100) if bench_size > 0 else 0.0

        rows.append(
            {
                "context_name": context_name,
                "context_bp": context_bp,
                "intersect_bp": intersect_bp,
                "pct_of_context": f"{pct_of_context:.6f}",
                "pct_of_bench": f"{pct_of_bench:.6f}",
            }
        )

        logger.info(
            f"  {context_name}: context={context_bp}, intersect={intersect_bp}, "
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
