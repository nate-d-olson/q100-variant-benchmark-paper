#!/usr/bin/env python3
# AI Disclosure: This script was developed with assistance from Claude (Anthropic).
"""
Generate variant Parquet table from annotated VCF using Truvari VariantRecord API.

Reads a fully annotated VCF (with INFO/CONTEXT_IDS and INFO/REGION_IDS from
bcftools annotate) and produces a clean Parquet table with:
- Correct variant type classification using VariantRecord.var_type()
- Variant sizes using VariantRecord.var_size()
- Size filtering (smvar <50bp, stvar >=50bp)
- Truvari size bins (SZBINS) as strings for R compatibility
- All columns required by R schema for caching and validation

This refactored version uses Truvari's VariantFile and VariantRecord classes
instead of vcf_to_df(), which provides:
- More robust variant classification (no manual enum mapping)
- Correct handling of all variant types (SNV, INDEL, INS, DEL)
- Simpler, more maintainable code
- Direct access to variant properties via tested methods

Output columns match R/schemas.R variant_table schema:
- bench_version, ref, bench_type (metadata)
- chrom, pos, end, gt (genomic coordinates and genotype)
- var_type, var_size, szbin (variant classification)
- ref_len, alt_len, qual, filter, is_pass (VCF quality info)
- context_ids, region_ids (annotations)
"""

import logging
import re
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import truvari


def parse_benchmark_id(benchmark_id: str) -> Dict[str, str]:
    """Parse benchmark ID into components.

    Args:
        benchmark_id: e.g. "v5.0q_GRCh38_smvar"

    Returns:
        Dict with bench_version, ref, bench_type
    """
    pattern = r"^([A-Za-z0-9.-]+)_([A-Za-z0-9.-]+)_(smvar|stvar)$"
    match = re.match(pattern, benchmark_id)
    if not match:
        raise ValueError(f"Invalid benchmark ID: {benchmark_id}")
    return {
        "bench_version": match.group(1),
        "ref": match.group(2),
        "bench_type": match.group(3),
    }


def get_size_bin(var_size: int) -> int:
    """Get Truvari size bin for variant size.

    Args:
        var_size: Variant size (can be negative for deletions)

    Returns:
        Size bin index (0-based)
    """
    abs_size = abs(var_size)
    for i, cutoff in enumerate(truvari.SZBINS):
        if abs_size < cutoff:
            return i
    return len(truvari.SZBINS)


def normalize_annotation(value) -> Optional[str]:
    """Normalize VCF INFO annotation to string.

    Handles Truvari's tuple/list returns for multi-valued INFO fields.

    Args:
        value: Raw INFO field value (can be None, str, tuple, or list)

    Returns:
        Normalized string or None
    """
    if value is None or value == ".":
        return None
    if isinstance(value, (tuple, list)):
        # Convert tuple/list to comma-separated string
        return ",".join(str(item) for item in value)
    return str(value)


def generate_variant_parquet(
    vcf_path: str,
    output_path: str,
    bench_type: str,
    benchmark_id: str,
    log_path: str,
) -> None:
    """Generate variant Parquet table from annotated VCF using VariantRecord API.

    Args:
        vcf_path: Path to fully annotated VCF (with CONTEXT_IDS, REGION_IDS)
        output_path: Path to output Parquet file
        bench_type: "smvar" or "stvar"
        benchmark_id: Full benchmark ID (e.g., "v5.0q_GRCh38_smvar")
        log_path: Path to log file
    """
    logging.basicConfig(
        filename=log_path,
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )

    logging.info(
        "Generating variant Parquet for %s (bench_type=%s)", benchmark_id, bench_type
    )
    logging.info("Input VCF: %s", vcf_path)
    logging.info("Output Parquet: %s", output_path)

    # Size threshold for smvar vs stvar
    size_threshold = 50

    # Open VCF with Truvari
    logging.info("Opening VCF with Truvari VariantFile...")
    vcf = truvari.VariantFile(vcf_path)

    # Get sample name (assuming single-sample VCF)
    samples = list(vcf.header.samples)
    if not samples:
        raise ValueError(f"VCF has no samples: {vcf_path}")
    sample = samples[0]
    logging.info("Processing sample: %s", sample)

    # Iterate through variants and collect data
    variants: List[Dict] = []
    filtered_size = 0
    filtered_non = 0

    logging.info("Processing variants...")
    for record in vcf:
        # Create VariantRecord wrapper
        vr = truvari.VariantRecord(record)

        # Get variant type and size using Truvari's methods
        var_type = vr.var_type()
        var_size = vr.var_size()

        # Filter out NON (non-variant) records
        if var_type == "NON":
            filtered_non += 1
            continue

        # Size filtering based on bench_type
        abs_size = abs(var_size)
        if bench_type == "smvar" and abs_size >= size_threshold:
            filtered_size += 1
            continue
        if bench_type == "stvar" and abs_size < size_threshold:
            filtered_size += 1
            continue

        # Get size bin (convert to string for R schema compatibility)
        szbin = str(get_size_bin(var_size))

        # Extract INFO annotations (handle tuple/list → string)
        context_ids = normalize_annotation(record.info.get("CONTEXT_IDS"))
        region_ids = normalize_annotation(record.info.get("REGION_IDS"))

        # Extract genotype using VariantRecord.gt()
        gt = vr.gt()

        # Extract allele lengths for R schema
        ref_len = len(record.ref)
        alt_len = len(record.alts[0]) if record.alts else None

        # Extract quality and filter for R schema
        qual = float(record.qual) if record.qual is not None else None
        filter_val = ";".join(record.filter.keys()) if record.filter else "PASS"
        is_pass = len(record.filter) == 0 or "PASS" in record.filter

        # Collect variant data (includes all columns from R schema)
        variants.append(
            {
                "chrom": record.chrom,
                "pos": record.pos,
                "end": vr.end,
                "gt": gt,
                "var_type": var_type,
                "var_size": var_size,
                "szbin": szbin,
                "ref_len": ref_len,
                "alt_len": alt_len,
                "qual": qual,
                "filter": filter_val,
                "is_pass": is_pass,
                "context_ids": context_ids,
                "region_ids": region_ids,
            }
        )

    vcf.close()

    logging.info("Collected %d variants", len(variants))
    if filtered_non > 0:
        logging.info("Filtered %d NON (non-variant) records", filtered_non)
    if filtered_size > 0:
        logging.info("Filtered %d variants by size (%s threshold)", filtered_size, bench_type)

    if not variants:
        logging.warning("No variants after filtering — creating empty DataFrame")
        # Create empty DataFrame with correct schema (matching R schema)
        df = pd.DataFrame(
            columns=[
                "chrom",
                "pos",
                "end",
                "gt",
                "var_type",
                "var_size",
                "szbin",
                "ref_len",
                "alt_len",
                "qual",
                "filter",
                "is_pass",
                "context_ids",
                "region_ids",
            ]
        )
    else:
        # Create DataFrame from collected variants
        df = pd.DataFrame(variants)

    # Add benchmark metadata columns at the beginning
    meta = parse_benchmark_id(benchmark_id)
    df.insert(0, "bench_version", meta["bench_version"])
    df.insert(1, "ref", meta["ref"])
    df.insert(2, "bench_type", meta["bench_type"])

    # Set proper dtypes for string columns
    for col in ["context_ids", "region_ids"]:
        if col in df.columns:
            df[col] = df[col].astype("string")

    # Log variant type distribution
    if not df.empty:
        logging.info(
            "Variant type distribution:\n%s", df["var_type"].value_counts().to_string()
        )
        logging.info("Size bin distribution:\n%s", df["szbin"].value_counts().to_string())
    else:
        logging.info("No variants to report statistics")

    # Optimize memory before writing
    if not df.empty:
        pre, post = truvari.optimize_df_memory(df)
        logging.info("Memory optimization: %d -> %d bytes", pre, post)

    # Write Parquet
    logging.info("Writing Parquet with %d variants...", len(df))
    table = pa.Table.from_pandas(df, preserve_index=False)

    pq.write_table(
        table,
        output_path,
        compression="zstd",
        compression_level=3,
        use_dictionary=True,
        write_statistics=True,
    )

    file_size_mb = Path(output_path).stat().st_size / (1024 * 1024)
    logging.info("Wrote %s (%.1f MB)", output_path, file_size_mb)
    logging.info("Complete!")


if __name__ == "__main__":
    benchmark_id = snakemake.wildcards.benchmark

    generate_variant_parquet(
        vcf_path=snakemake.input.vcf,
        output_path=snakemake.output.parquet,
        bench_type=snakemake.params.bench_type,
        benchmark_id=benchmark_id,
        log_path=snakemake.log[0],
    )
