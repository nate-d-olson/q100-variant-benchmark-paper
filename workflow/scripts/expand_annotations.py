#!/usr/bin/env python3
"""
Expand STRAT_IDS and REGION_IDS columns into binary flags in TSV.
Reads from stdin, writes to stdout.

Usage:
    python expand_annotations.py --strat-ids ID1 ID2 ... --region-ids ID3 ID4 ...
"""
import argparse
import sys
import csv


def main():
    parser = argparse.ArgumentParser(description='Expand annotation IDs into binary columns')
    parser.add_argument('--strat-ids', nargs='*', default=[],
                        help='Stratification IDs to expand')
    parser.add_argument('--region-ids', nargs='*', default=[],
                        help='Region IDs to expand')
    args = parser.parse_args()

    reader = csv.DictReader(sys.stdin, delimiter='\t')

    # Prepare new fieldnames
    fieldnames = list(reader.fieldnames)

    # Remove container columns
    if 'STRAT_IDS' in fieldnames:
        fieldnames.remove('STRAT_IDS')
    if 'REGION_IDS' in fieldnames:
        fieldnames.remove('REGION_IDS')

    # Add expansion columns
    # Stratifications: STRAT_{ID}
    strat_cols = [f"STRAT_{sid}" for sid in args.strat_ids]
    fieldnames.extend(strat_cols)

    # Regions: {ID} (BMKREGIONS or EXCL_...)
    region_cols = args.region_ids
    fieldnames.extend(region_cols)

    writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames, delimiter='\t')
    writer.writeheader()

    for row in reader:
        # Parse STRAT_IDS
        strat_ids_str = row.pop('STRAT_IDS', '')
        strats = set(strat_ids_str.split(',')) if strat_ids_str else set()
        for sid, col in zip(args.strat_ids, strat_cols):
            row[col] = '1' if sid in strats else '0'

        # Parse REGION_IDS
        region_ids_str = row.pop('REGION_IDS', '')
        regions = set(region_ids_str.split(',')) if region_ids_str else set()
        for rid in args.region_ids:
            row[rid] = '1' if rid in regions else '0'

        writer.writerow(row)


if __name__ == '__main__':
    main()
