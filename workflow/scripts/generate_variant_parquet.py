#!/usr/bin/env python3
# AI Disclosure: This script was developed with assistance from Claude (Anthropic).
"""
Generate variant Parquet table from annotated VCF using Truvari VariantRecord API.

Reads a fully annotated VCF (with INFO/CONTEXT_IDS and INFO/REGION_IDS from
bcftools annotate) and produces a clean Parquet table with:
- Correct variant type classification using VariantRecord.var_type()
- Variant sizes using VariantRecord.var_size()
- Genotype classification using get_gt()
- Size filtering (smvar <50bp, stvar >=50bp)
- Truvari size bins (SZBINS) as strings for R compatibility
- All columns required by R schema for caching and validation

This refactored version uses Truvari's VariantFile and VariantRecord classes
instead of vcf_to_df(), which provides:
- More robust variant classification (no manual enum mapping)
- Correct handling of all variant types (SNV, INDEL, INS, DEL)
- Simpler, more maintainable code
- Direct access to variant properties via tested methods

Truvari Enum Conversions:
- var_type: SV enum → .name extracts string ("DEL", "INS", "SNP")
- gt: tuple → get_gt() → GT enum → .name extracts string ("HET", "HOM", "REF")
- szbin: get_sizebin() returns string directly ("SNP", "[50,100)", etc.)

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

    # Open VCF with Truvari using context manager for proper cleanup
    logging.info("Opening VCF with Truvari VariantFile...")
    with truvari.VariantFile(vcf_path) as vcf:
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

            # Get variant type (as string from enum)
            var_type = vr.var_type().name

            # Get variant size (Truvari returns unsigned/absolute value)
            abs_size = vr.var_size()

            # Get reference and alternate allele lengths
            ref_len = len(vr.get_ref())
            alt_allele = vr.get_alt()
            if alt_allele is None:
                logging.warning(
                    "Skipping %s variant with missing ALT at %s:%d",
                    var_type,
                    record.chrom,
                    record.pos,
                )
                continue
            alt_len = len(alt_allele)
            
            # Apply sign convention: positive for INS, negative for DEL
            if var_type == "DEL":
                var_size = -abs_size
            elif var_type == "INS":
                var_size = abs_size
            elif var_type in ["UNK", "INV"]:
                # Recalculate size with sign for ambiguous types
                logging.info(
                    "Reclassifying %s to INS/DEL at %s:%d",
                    var_type,
                    record.chrom,
                    record.pos,
                )
                var_size = alt_len - ref_len  # Positive for INS, negative for DEL
                var_type = "INS" if var_size > 0 else "DEL"
            else:
                # Other types (SNP, DUP, BND) keep unsigned size
                var_size = abs_size

            # Convert INS/DEL to INDEL for small variants
            if bench_type == "smvar" and var_type in ["INS", "DEL"]:
                var_type = "INDEL"

                # Size filtering based on bench_type
                abs_size = abs(var_size)
                if bench_type == "smvar" and abs_size >= size_threshold:
                    filtered_size += 1
                    continue
                if bench_type == "stvar" and abs_size < size_threshold:
                    filtered_size += 1
                    continue

            # Get size bin (convert to string for R schema compatibility)
            szbin = truvari.get_sizebin(abs_size)

            # Extract INFO annotations (handle tuple/list → string)
            context_ids = normalize_annotation(record.info.get("CONTEXT_IDS"))
            region_ids = normalize_annotation(record.info.get("REGION_IDS"))

            # Extract genotype: gt() returns tuple, convert to GT enum, then extract name
            gt = truvari.get_gt(vr.gt()).name  # "HET", "HOM", "REF", "NON", "UNK"

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

    logging.info("Collected %d variants", len(variants))
    if filtered_non > 0:
        logging.info("Filtered %d NON (non-variant) records", filtered_non)
    if filtered_size > 0:
        logging.info(
            "Filtered %d variants by size (%s threshold)", filtered_size, bench_type
        )

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
        logging.info(
            "Size bin distribution:\n%s", df["szbin"].value_counts().to_string()
        )
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
