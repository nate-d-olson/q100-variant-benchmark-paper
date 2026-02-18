#!/usr/bin/env python3
# AI Disclosure: This script was developed with assistance from Claude (Anthropic)
"""
find_chr8_inversion.py — detect the largest inversion from SyRI output.

Parses a SyRI output file to find the largest inversion (by query size) on the
specified chromosome, then writes reference/query coordinates together with
pre-computed zoom region boundaries to a JSON file.

The JSON file is consumed by make_chr8_figure.py and recorded as an explicit
pipeline output so the detected coordinates are traceable.

Usage
-----
  python find_chr8_inversion.py \\
      --syri-rp  results/chr8_synteny/syri/ref_patsyri.out \\
      --chrom    chr8 \\
      --min-inv-size 0 \\
      --zoom-padding 2000000 \\
      --output   results/chr8_synteny/inversion_coords.json
"""

import argparse
import json
import sys
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--syri-rp",
        required=True,
        metavar="FILE",
        help="SyRI output file with PAT as query (ref_patsyri.out)",
    )
    p.add_argument("--chrom", default="chr8", help="Chromosome name (default: chr8)")
    p.add_argument(
        "--min-inv-size",
        type=int,
        default=0,
        metavar="BP",
        help="Minimum inversion size to consider (default: 0)",
    )
    p.add_argument(
        "--zoom-padding",
        type=int,
        default=2_000_000,
        metavar="BP",
        help="Flanking bp to add around the inversion for the zoom region (default: 2000000)",
    )
    p.add_argument(
        "--output",
        required=True,
        metavar="FILE",
        help="Output JSON file path",
    )
    return p.parse_args()


def find_largest_inversion(syri_path: str, chrom: str, min_size: int = 0) -> dict:
    """
    Parse a SyRI output file and return coordinates of the largest INV call.

    SyRI column layout (1-based):
      1:  ref chromosome
      2-3: ref interval (start, end)
      7-8: query interval (start, end)
      11:  annotation type (INV, SYN, ...)

    Parameters
    ----------
    syri_path:
        Path to a SyRI .out file.
    chrom:
        Chromosome name to filter on (must match column 1).
    min_size:
        Minimum inversion size in bp (by query interval); smaller calls are skipped.

    Returns
    -------
    dict with keys ref_start, ref_end, qry_start, qry_end, size (all int).

    Raises
    ------
    RuntimeError if no qualifying inversion is found.
    """
    best = None
    best_score = None

    with open(syri_path) as fh:
        for line_no, line in enumerate(fh, start=1):
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11 or fields[10] != "INV" or fields[0] != chrom:
                continue

            try:
                ref_start, ref_end = sorted((int(fields[1]), int(fields[2])))
                qry_start, qry_end = sorted((int(fields[6]), int(fields[7])))
            except ValueError as exc:
                raise ValueError(
                    f"Invalid SyRI coordinates in {syri_path} at line {line_no}"
                ) from exc

            inv_size = qry_end - qry_start
            if inv_size < min_size:
                continue

            # Primary sort: query size; secondary: ref size (both descending)
            score = (inv_size, ref_end - ref_start)
            if best_score is None or score > best_score:
                best_score = score
                best = {
                    "ref_start": ref_start,
                    "ref_end": ref_end,
                    "qry_start": qry_start,
                    "qry_end": qry_end,
                    "size": inv_size,
                }

    if best is None:
        raise RuntimeError(
            f"No inversion found in {syri_path} for chrom '{chrom}' "
            f"with min size {min_size} bp"
        )

    return best


def main():
    args = parse_args()

    if not Path(args.syri_rp).exists():
        print(f"[ERROR] SyRI file not found: {args.syri_rp}", file=sys.stderr)
        sys.exit(1)

    print(f"[INFO] Parsing SyRI file: {args.syri_rp}")
    print(f"[INFO] Chromosome: {args.chrom}, min size: {args.min_inv_size} bp")

    inv = find_largest_inversion(args.syri_rp, args.chrom, args.min_inv_size)

    zoom_ref_start = max(0, inv["ref_start"] - args.zoom_padding)
    zoom_ref_end = inv["ref_end"] + args.zoom_padding

    coords = {
        "chrom": args.chrom,
        "ref_start": inv["ref_start"],
        "ref_end": inv["ref_end"],
        "qry_start": inv["qry_start"],
        "qry_end": inv["qry_end"],
        "size": inv["size"],
        "zoom_ref_start": zoom_ref_start,
        "zoom_ref_end": zoom_ref_end,
    }

    print(
        f"[INFO] Largest inversion — "
        f"REF {args.chrom}:{inv['ref_start']}-{inv['ref_end']}, "
        f"PAT {args.chrom}:{inv['qry_start']}-{inv['qry_end']}, "
        f"size {inv['size']:,} bp"
    )
    print(
        f"[INFO] Zoom region (REF) — "
        f"{args.chrom}:{zoom_ref_start}-{zoom_ref_end} "
        f"(padding ±{args.zoom_padding:,} bp)"
    )

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as fh:
        json.dump(coords, fh, indent=2)

    print(f"[INFO] Coordinates written to {args.output}")


if __name__ == "__main__":
    main()
