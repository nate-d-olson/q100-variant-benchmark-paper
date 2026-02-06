#!/usr/bin/env python3
"""
Generate VCF INFO header lines for bcftools annotate.

Usage:
    python generate_header_lines.py \
        --output header_lines.txt
"""

import argparse


def main():
    parser = argparse.ArgumentParser(description="Generate VCF header lines")
    parser.add_argument("--output", required=True, help="Output header file")
    args = parser.parse_args()

    lines = []

    # Genomic Context IDs field (comma-separated list of genomic context region IDs)
    lines.append(
        "##INFO=<ID=CONTEXT_IDS,Number=.,Type=String,"
        'Description="Comma-separated list of genomic context region IDs overlapping variant">'
    )

    # Region IDs field (includes Benchmark Regions and Exclusions)
    lines.append(
        "##INFO=<ID=REGION_IDS,Number=.,Type=String,"
        'Description="Comma-separated list of region IDs (Benchmark, Exclusions) overlapping variant">'
    )

    with open(args.output, "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
