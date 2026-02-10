#!/usr/bin/env python3
# AI Disclosure: This script was modified with assistance from Claude (Anthropic)
# for bug fix and code review.
"""
Compute overlap metrics between two BED files.

Uses bedtools to calculate:
- size_a: Total bases in BED A (after merge)
- size_b: Total bases in BED B (after merge)
- intersect_bp: Bases overlapping between A and B
- pct_of_a: Percent of A covered by B
- pct_of_b: Percent of B covered by A
"""

import subprocess
from logging_config import setup_logger

logger = setup_logger(__name__)


def compute_bed_size(bed_path: str, is_gzipped: bool = False) -> int:
    """
    Compute total size of BED after merging overlaps.

    Args:
        bed_path: Path to BED file
        is_gzipped: Whether the BED file is gzipped

    Returns:
        Total base pairs in merged BED regions
    """
    if is_gzipped:
        cmd = (
            f"gzip -dc {bed_path} | "
            f"bedtools sort -i - | "
            f"bedtools merge -i - | "
            f"awk '{{sum+=$3-$2}} END {{print sum+0}}'"
        )
    else:
        cmd = (
            f"bedtools sort -i {bed_path} | "
            f"bedtools merge -i - | "
            f"awk '{{sum+=$3-$2}} END {{print sum+0}}'"
        )

    result = subprocess.run(
        cmd, shell=True, capture_output=True, text=True, executable="/bin/bash"
    )
    if result.returncode != 0:
        raise RuntimeError(f"Failed to compute BED size: {result.stderr}")

    return int(result.stdout.strip() or 0)


def compute_intersection(
    bed_a: str, bed_b: str, a_gzipped: bool = False, b_gzipped: bool = False
) -> int:
    """
    Compute overlap between two BED files.

    Args:
        bed_a: Path to first BED file
        bed_b: Path to second BED file
        a_gzipped: Whether bed_a is gzipped
        b_gzipped: Whether bed_b is gzipped

    Returns:
        Total base pairs of overlap
    """
    if a_gzipped:
        a_cmd = f"gzip -dc {bed_a} | bedtools sort -i -"
    else:
        a_cmd = f"bedtools sort -i {bed_a}"

    if b_gzipped:
        b_cmd = f"gzip -dc {bed_b} | bedtools sort -i -"
    else:
        b_cmd = f"bedtools sort -i {bed_b}"

    cmd = (
        f"bedtools intersect -a <({a_cmd}) -b <({b_cmd}) | "
        f"bedtools sort -i - | bedtools merge -i - | "
        f"awk '{{sum+=$3-$2}} END {{print sum+0}}'"
    )

    result = subprocess.run(
        cmd, shell=True, capture_output=True, text=True, executable="/bin/bash"
    )
    if result.returncode != 0:
        raise RuntimeError(f"Failed to compute intersection: {result.stderr}")

    return int(result.stdout.strip() or 0)


def main():
    """Main entry point for Snakemake script execution."""
    if len(snakemake.input) < 2:
        raise ValueError("Expected at least 2 input BED files")

    bed_a = snakemake.input[0]
    bed_b = snakemake.input[1]

    a_gzipped = bed_a.endswith(".gz")
    b_gzipped = bed_b.endswith(".gz")

    logger.info(f"Computing metrics for {bed_a} vs {bed_b}")

    size_a = compute_bed_size(bed_a, a_gzipped)
    size_b = compute_bed_size(bed_b, b_gzipped)
    intersect = compute_intersection(bed_a, bed_b, a_gzipped, b_gzipped)

    pct_of_a = (intersect / size_a * 100) if size_a > 0 else 0.0
    pct_of_b = (intersect / size_b * 100) if size_b > 0 else 0.0

    # Get region name from wildcards (different field names in different rules)
    region_name = snakemake.wildcards.get("exclusion") or snakemake.wildcards.get(
        "genomic_context"
    )

    with open(snakemake.output[0], "w") as f:
        f.write(
            f"{region_name}\t{size_a}\t{intersect}\t{pct_of_a:.6f}\t{pct_of_b:.6f}\n"
        )

    logger.info(
        f"Metrics: size_a={size_a}, intersect={intersect}, "
        f"pct_of_a={pct_of_a:.2f}%, pct_of_b={pct_of_b:.2f}%"
    )


if __name__ == "__main__":
    main()
