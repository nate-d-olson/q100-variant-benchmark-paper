#!/usr/bin/env python3
"""
Combine multiple BED files into a single sorted file with ID column.
Flattens overlapping intervals and merges IDs.

For use with bcftools annotate --columns CHROM,FROM,TO,<ID>

Usage:
    python combine_beds_with_id.py \
        --beds file1.bed:ID1 file2.bed:ID2 ... \
        --output combined.bed
"""
import argparse
import gzip
import sys
from collections import defaultdict


def open_bed(path):
    """Open BED file, handling gzip compression."""
    if str(path).endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'r')


def flatten_intervals(intervals):
    """
    Flatten overlapping intervals and merge IDs using sweep line algorithm.

    Args:
        intervals: list of (start, end, id) tuples

    Returns:
        list of (start, end, id_string) tuples with flattened intervals
    """
    if not intervals:
        return []

    # Create events for sweep line algorithm
    # BED is 0-based, half-open [start, end)
    events = []
    for start, end, annot_id in intervals:
        events.append((start, 1, annot_id))  # Start event
        events.append((end, -1, annot_id))   # End event

    # Sort events by position, with start events (1) before end events (-1) at same position
    events.sort(key=lambda x: (x[0], -x[1]))

    current_ids = set()
    result = []

    for i in range(len(events)):
        pos, event_type, annot_id = events[i]

        # Process event
        if event_type == 1:
            current_ids.add(annot_id)
        else:
            current_ids.discard(annot_id)

        # Get next position
        if i + 1 < len(events):
            next_pos = events[i + 1][0]
        else:
            next_pos = pos

        # Emit interval if:
        # 1. There's a gap to the next position
        # 2. There are active IDs
        if next_pos > pos and current_ids:
            ids_str = ",".join(sorted(list(current_ids)))
            # Merge with previous interval if IDs are identical and intervals are adjacent
            if result and result[-1][2] == ids_str and result[-1][1] == pos:
                result[-1] = (result[-1][0], next_pos, ids_str)
            else:
                result.append((pos, next_pos, ids_str))

    return result


def main():
    parser = argparse.ArgumentParser(description='Combine BED files with ID column')
    parser.add_argument('--beds', nargs='+', required=True,
                        help='BED files with IDs in format path:ID')
    parser.add_argument('--output', required=True, help='Output BED file')
    args = parser.parse_args()

    intervals_by_chrom = defaultdict(list)

    # Read all input BED files
    for bed_spec in args.beds:
        if ':' not in bed_spec:
            print(f"ERROR: Invalid bed spec '{bed_spec}'. Expected format: path:ID", file=sys.stderr)
            sys.exit(1)

        path, annot_id = bed_spec.rsplit(':', 1)

        try:
            with open_bed(path) as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    fields = line.strip().split('\t')
                    if len(fields) >= 3:
                        chrom, start, end = fields[0], int(fields[1]), int(fields[2])
                        intervals_by_chrom[chrom].append((start, end, annot_id))
        except FileNotFoundError:
            print(f"ERROR: File not found: {path}", file=sys.stderr)
            sys.exit(1)
        except Exception as e:
            print(f"ERROR: Failed to read {path}: {e}", file=sys.stderr)
            sys.exit(1)

    # Sort chromosomes
    def sort_key(chrom):
        """Sort chromosomes numerically if possible, otherwise alphabetically."""
        if chrom.startswith('chr'):
            chrom_num = chrom[3:]
        else:
            chrom_num = chrom
        try:
            return (0, int(chrom_num))
        except ValueError:
            return (1, chrom_num)

    sorted_chroms = sorted(intervals_by_chrom.keys(), key=sort_key)

    # Write output
    with open(args.output, 'w') as out:
        for chrom in sorted_chroms:
            flattened = flatten_intervals(intervals_by_chrom[chrom])
            for start, end, ids in flattened:
                out.write(f'{chrom}\t{start}\t{end}\t{ids}\n')


if __name__ == '__main__':
    main()
