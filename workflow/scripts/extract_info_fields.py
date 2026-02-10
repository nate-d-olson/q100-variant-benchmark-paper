#!/usr/bin/env python3
"""
Extract INFO field names from VCF header for bcftools query.

Usage:
    python extract_info_fields.py --vcf input.vcf.gz --output fields.txt
"""

import argparse
import gzip
import re
import sys


def main():
    parser = argparse.ArgumentParser(description="Extract INFO fields from VCF")
    parser.add_argument("--vcf", required=True, help="Input VCF file")
    parser.add_argument("--output", required=True, help="Output file with field names")
    parser.add_argument(
        "--exclude", nargs="*", default=[], help="INFO fields to exclude from output"
    )
    args = parser.parse_args()

    info_fields = []
    opener = gzip.open if args.vcf.endswith(".gz") else open

    try:
        with opener(args.vcf, "rt") as f:
            for line in f:
                if not line.startswith("##"):
                    break  # End of header
                if line.startswith("##INFO="):
                    match = re.search(r"ID=([^,>]+)", line)
                    if match:
                        field_id = match.group(1)
                        if field_id not in args.exclude:
                            info_fields.append(field_id)
    except FileNotFoundError:
        print(f"ERROR: File not found: {args.vcf}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: Failed to read VCF: {e}", file=sys.stderr)
        sys.exit(1)

    with open(args.output, "w") as f:
        f.write("\n".join(info_fields))


if __name__ == "__main__":
    main()
