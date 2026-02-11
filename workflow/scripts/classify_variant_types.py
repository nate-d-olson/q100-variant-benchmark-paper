#!/usr/bin/env python3
"""
Reclassify variant types based on REF/ALT lengths and SVTYPE/SVLEN for SVs.

Bcftools %TYPE only provides basic classification (snp/indel/mnp) and doesn't
distinguish between insertions and deletions. This script reads the variant table,
classifies variants based on size and available fields, and writes an updated table.

Classification logic:
- For variants >50bp with SVTYPE: Use SVTYPE field (DEL, INS, DUP, INV, BND, CNV, etc.)
- For variants ≤50bp or missing SVTYPE:
  - REF_len == ALT_len == 1: SNP
  - REF_len == ALT_len > 1: MNP (multi-nucleotide polymorphism)
  - REF_len < ALT_len: INS (insertion)
  - REF_len > ALT_len: DEL (deletion)
  - Other cases: COMPLEX
"""

import csv
import sys
from pathlib import Path
from typing import Optional, TextIO


def classify_small_variant(ref_len: int, alt_len: int) -> str:
    """
    Classify small variant (<= 50bp) based on REF and ALT lengths.

    Args:
        ref_len: Length of REF allele
        alt_len: Length of ALT allele

    Returns:
        Variant type: SNP, MNP, INS, DEL, or COMPLEX
    """
    # Equal lengths
    if ref_len == alt_len:
        if ref_len == 1:
            return "SNP"
        else:
            return "MNP"  # Multi-nucleotide polymorphism

    # Insertion (ALT longer than REF)
    elif alt_len > ref_len:
        return "INS"

    # Deletion (REF longer than ALT)
    elif ref_len > alt_len:
        return "DEL"

    # Shouldn't reach here, but fallback to COMPLEX
    return "COMPLEX"


def classify_variant_type(
    ref_len: int,
    alt_len: int,
    svtype: Optional[str] = None,
    svlen: Optional[int] = None
) -> str:
    """
    Classify variant type using size-dependent logic.

    For variants >50bp with SVTYPE annotation, use SVTYPE directly.
    For smaller variants or those without SVTYPE, classify by REF/ALT lengths.

    Args:
        ref_len: Length of REF allele
        alt_len: Length of ALT allele
        svtype: SVTYPE INFO field value (if present)
        svlen: SVLEN INFO field value (if present)

    Returns:
        Variant type: SNP, MNP, INS, DEL, COMPLEX, or SV type (DUP, INV, BND, CNV)
    """
    # Determine variant size
    if svlen is not None:
        var_size = abs(svlen)
    else:
        var_size = abs(alt_len - ref_len)

    # For structural variants (>50bp) with SVTYPE, use it directly
    if var_size > 50 and svtype:
        # Normalize SVTYPE to uppercase
        svtype_upper = svtype.upper()
        # Common SV types: DEL, INS, DUP, INV, BND, CNV
        return svtype_upper

    # For small variants or those without SVTYPE, use length-based classification
    return classify_small_variant(ref_len, alt_len)


def reclassify_variant_table(
    input_tsv: Path,
    output_tsv: Path,
    log_file: TextIO
) -> None:
    """
    Read variant table, reclassify types, and write updated table.

    Uses SVTYPE/SVLEN for structural variants (>50bp) when available,
    and REF/ALT length-based classification for small variants.

    Args:
        input_tsv: Path to input variant table
        output_tsv: Path to output variant table
        log_file: Open log file handle
    """
    log_file.write("Reclassifying variant types\n")
    log_file.write(f"Input: {input_tsv}\n")
    log_file.write(f"Output: {output_tsv}\n\n")

    # Track statistics (use defaultdict for dynamic SV types)
    from collections import defaultdict
    stats = defaultdict(int)
    stats["total"] = 0
    sv_count = 0  # Track variants classified using SVTYPE

    try:
        with open(input_tsv, "r") as infile, open(output_tsv, "w") as outfile:
            reader = csv.DictReader(infile, delimiter="\t")
            fieldnames = reader.fieldnames or []

            log_file.write(f"Found {len(fieldnames)} columns\n")

            # Find required columns
            def find_col(candidates):
                for f in fieldnames:
                    # Handle bcftools [n] prefix
                    norm = f.split("]", 1)[1] if "]" in f else f
                    if norm in candidates:
                        return f
                return None

            type_col = find_col({"TYPE"})
            ref_len_col = find_col({"STRLEN(REF)"})
            alt_len_col = find_col({"STRLEN(ALT)"})
            svtype_col = find_col({"SVTYPE", "INFO/SVTYPE"})
            svlen_col = find_col({"SVLEN", "INFO/SVLEN"})

            if not type_col or not ref_len_col or not alt_len_col:
                missing = []
                if not type_col:
                    missing.append("TYPE")
                if not ref_len_col:
                    missing.append("STRLEN(REF)")
                if not alt_len_col:
                    missing.append("STRLEN(ALT)")
                raise ValueError(f"Missing required columns: {missing}")

            log_msg = (
                f"Using columns: TYPE='{type_col}', REF_LEN='{ref_len_col}', "
                f"ALT_LEN='{alt_len_col}'"
            )
            if svtype_col:
                log_msg += f", SVTYPE='{svtype_col}'"
            if svlen_col:
                log_msg += f", SVLEN='{svlen_col}'"
            log_file.write(log_msg + "\n\n")

            # Write header
            writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()

            # Process rows
            for row in reader:
                stats["total"] += 1

                try:
                    ref_len = int(row[ref_len_col])
                    alt_len = int(row[alt_len_col])

                    # Extract SVTYPE and SVLEN if present
                    svtype = None
                    svlen = None

                    if svtype_col and row.get(svtype_col):
                        svtype_val = row[svtype_col].strip()
                        # Handle missing values (., empty string)
                        if svtype_val and svtype_val != ".":
                            svtype = svtype_val

                    if svlen_col and row.get(svlen_col):
                        svlen_val = row[svlen_col].strip()
                        # Handle missing values (., empty string)
                        if svlen_val and svlen_val != ".":
                            try:
                                svlen = int(svlen_val)
                            except ValueError:
                                pass  # Keep svlen as None

                    # Determine variant size for tracking
                    if svlen is not None:
                        var_size = abs(svlen)
                    else:
                        var_size = abs(alt_len - ref_len)

                    # Reclassify
                    new_type = classify_variant_type(ref_len, alt_len, svtype, svlen)
                    row[type_col] = new_type
                    stats[new_type] += 1

                    # Track if classified as SV using SVTYPE
                    if var_size > 50 and svtype:
                        sv_count += 1

                except (ValueError, KeyError) as e:
                    log_file.write(f"WARNING: Skipping row due to error: {e}\n")
                    continue

                writer.writerow(row)

                if stats["total"] % 100000 == 0:
                    log_file.write(f"Processed {stats['total']:,} variants...\n")

        # Write summary
        log_file.write("\nReclassification complete:\n")
        log_file.write(f"  Total variants: {stats['total']:,}\n")
        log_file.write(f"  Classified using SVTYPE: {sv_count:,}\n\n")

        # Sort variant types for consistent output
        type_order = ["SNP", "MNP", "INS", "DEL", "COMPLEX"]
        sv_types = sorted([k for k in stats.keys() if k not in type_order and k != "total"])

        for vtype in type_order:
            if vtype in stats:
                log_file.write(f"  {vtype}: {stats[vtype]:,}\n")

        if sv_types:
            log_file.write("\n  Structural variant types:\n")
            for vtype in sv_types:
                log_file.write(f"  {vtype}: {stats[vtype]:,}\n")

    except Exception as e:
        log_file.write(f"ERROR: {str(e)}\n")
        raise


if __name__ == "__main__":
    # Snakemake provides these variables
    input_table = Path(snakemake.input.tsv)
    output_table = Path(snakemake.output.tsv)

    with open(snakemake.log[0], "w") as log:
        reclassify_variant_table(input_table, output_table, log)
