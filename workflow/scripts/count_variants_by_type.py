#!/usr/bin/env python3
"""
Count variants by type from bcftools query output.

Reads variant types from stdin (one per line), counts occurrences,
and outputs long-format TSV with all variant types including zeros.
Auto-detects variant types present in the data for transparency.

Usage:
    bcftools query -f '%TYPE\n' input.vcf | python count_variants_by_type.py benchmark ref context
    bcftools query -f '%SVTYPE\n' input.vcf | python count_variants_by_type.py benchmark ref context
"""

import sys
from collections import Counter


def main():
    if len(sys.argv) != 4:
        print("Usage: count_variants_by_type.py <benchmark> <ref> <context>", file=sys.stderr)
        sys.exit(1)
    
    benchmark = sys.argv[1]
    ref = sys.argv[2]
    context = sys.argv[3]
    
    # Read variant types from stdin
    variant_types = []
    for line in sys.stdin:
        line = line.strip()
        if line:  # Skip empty lines
            variant_types.append(line)
    
    # Count variant types
    type_counts = Counter(variant_types)
    
    # Get all unique types (auto-detected)
    all_types = sorted(type_counts.keys()) if type_counts else []
    
    # If no variants found, output header only with note
    if not all_types:
        print("benchmark\tref\tcontext\tvariant_type\tcount")
        print(f"{benchmark}\t{ref}\t{context}\tNO_VARIANTS\t0")
        return
    
    # Output long-format TSV
    print("benchmark\tref\tcontext\tvariant_type\tcount")
    for variant_type in all_types:
        count = type_counts.get(variant_type, 0)
        print(f"{benchmark}\t{ref}\t{context}\t{variant_type}\t{count}")


if __name__ == "__main__":
    main()
