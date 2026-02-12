#!/usr/bin/env python3
# AI Disclosure: This script was developed with assistance from Claude (Anthropic).
"""
Generate variant Parquet table from annotated VCF using Truvari.

Reads a fully annotated VCF (with INFO/CONTEXT_IDS and INFO/REGION_IDS from
bcftools annotate) and produces a clean Parquet table with:
- Correct variant type classification (SNV vs INDEL for smvar; INS/DEL for stvar)
- Size filtering (smvar <50bp, stvar >=50bp)
- Truvari size bins (SZBINS)
- All INFO fields included automatically (no version-specific column lists)

Replaces: generate_var_table (bcftools query), extract_info_fields,
          and the R-side tidy_smvar/tidy_stvar/get_bench_var_cols functions.
"""

import logging
import re
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import truvari


def parse_benchmark_id(benchmark_id: str) -> dict:
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


def classify_smvar(svtype_series: pd.Series) -> pd.Series:
    """Classify small variant types.

    SNP (enum 0) -> SNV, everything else -> INDEL.
    """
    return svtype_series.map(lambda x: "SNV" if x == 0 else "INDEL")


def classify_stvar(
    svtype_series: pd.Series, ref_len_series: pd.Series, alt_len_series: pd.Series
) -> pd.Series:
    """Classify structural variant types.

    DEL (1) -> DEL, INS (2) -> INS.
    DUP (3), INV (4), UNK (6) are reclassified as INS or DEL based on
    whether ALT is longer (INS) or REF is longer (DEL).
    """
    # Start with direct mapping for DEL and INS
    result = svtype_series.map({0: "SNV", 1: "DEL", 2: "INS", 5: "NON"})

    # Reclassify DUP, INV, UNK based on allele lengths
    needs_reclassify = svtype_series.isin([3, 4, 6])
    if needs_reclassify.any():
        n_reclassified = needs_reclassify.sum()
        logging.info(
            "Reclassifying %d variants (DUP=%d, INV=%d, UNK=%d) as INS/DEL by allele size",
            n_reclassified,
            (svtype_series[needs_reclassify] == 3).sum(),
            (svtype_series[needs_reclassify] == 4).sum(),
            (svtype_series[needs_reclassify] == 6).sum(),
        )
        result[needs_reclassify] = pd.Series(
            [
                "INS" if alt > ref else "DEL"
                for alt, ref in zip(
                    alt_len_series[needs_reclassify], ref_len_series[needs_reclassify]
                )
            ],
            index=result[needs_reclassify].index,
        )

    # Fill any remaining NaN (shouldn't happen, but defensive)
    result = result.fillna("DEL")
    return result


def generate_variant_parquet(
    vcf_path: str,
    output_path: str,
    bench_type: str,
    benchmark_id: str,
    log_path: str,
) -> None:
    """Generate variant Parquet table from annotated VCF.

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

    # Read VCF with Truvari — gets all INFO fields, FORMAT fields, and alleles
    logging.info("Reading VCF with truvari.vcf_to_df()...")
    df = truvari.vcf_to_df(
        str(vcf_path), with_info=True, with_format=True, no_prefix=True, alleles=True
    )
    logging.info("Loaded %d raw variants with columns: %s", len(df), list(df.columns))

    # Filter to variants in benchmark regions (REGION_IDS contains BMKREGIONS)
    if "REGION_IDS" in df.columns:
        before = len(df)
        df = df[df["REGION_IDS"].str.contains("BMKREGIONS", na=False)].copy()
        logging.info("Filtered to benchmark regions: %d -> %d", before, len(df))
    else:
        logging.warning("No REGION_IDS column found — keeping all variants")

    # Filter out non-variant genotypes (0/0, ./.)
    # Truvari's is_pass and svtype handle most of this, but explicit filter
    # for NON (monomorphic ref) variants
    before = len(df)
    df = df[df["svtype"] != 5].copy()  # SV.NON = 5
    if len(df) < before:
        logging.info("Removed %d non-variant (NON) records", before - len(df))

    # Size filtering
    logging.info("Applying size filter for %s...", bench_type)
    abs_svlen = df["svlen"].abs()
    before = len(df)
    if bench_type == "smvar":
        df = df[abs_svlen < 50].copy()
    else:
        df = df[abs_svlen >= 50].copy()
    logging.info("Size filtering: %d -> %d variants", before, len(df))

    # Compute ref_len and alt_len from allele strings before classification
    # (stvar classification needs allele lengths to reclassify DUP/INV/UNK)
    # Note: Truvari's "ref" column is the REF allele, not the reference genome
    if "ref" in df.columns and "alt" in df.columns:
        df["ref_len"] = df["ref"].str.len().astype("Int32")
        df["alt_len"] = df["alt"].str.len().astype("Int32")
        df.drop(columns=["ref", "alt"], inplace=True)
    else:
        logging.warning("No ref/alt allele columns — ref_len/alt_len will be missing")

    # Classify variant types
    logging.info("Classifying variant types...")
    if bench_type == "smvar":
        df["var_type"] = classify_smvar(df["svtype"])
    else:
        df["var_type"] = classify_stvar(df["svtype"], df["ref_len"], df["alt_len"])

    # Rename Truvari columns to project conventions
    rename_map = {
        "svlen": "var_size",  # project uses var_size
    }
    # Truvari uses 'start' for 0-based position; project uses 'pos'
    if "start" in df.columns:
        rename_map["start"] = "pos"
    # Normalize INFO field names to lowercase
    if "CONTEXT_IDS" in df.columns:
        rename_map["CONTEXT_IDS"] = "context_ids"
    if "REGION_IDS" in df.columns:
        rename_map["REGION_IDS"] = "region_ids"
    # Rename GT format field to lowercase
    if "GT" in df.columns:
        rename_map["GT"] = "gt"
    df.rename(columns=rename_map, inplace=True)

    # Drop the raw svtype enum — var_type is the project classification
    if "svtype" in df.columns:
        df.drop(columns=["svtype"], inplace=True)

    # Clean up context_ids and region_ids: replace '.' with None
    for col in ["context_ids", "region_ids"]:
        if col in df.columns:
            df[col] = df[col].replace({".": pd.NA, "": pd.NA}).astype("string")

    # Add benchmark metadata columns
    meta = parse_benchmark_id(benchmark_id)
    df.insert(0, "bench_version", meta["bench_version"])
    df.insert(1, "ref", meta["ref"])
    df.insert(2, "bench_type", meta["bench_type"])

    # Log variant type distribution
    logging.info(
        "Variant type distribution:\n%s", df["var_type"].value_counts().to_string()
    )
    logging.info("Size bin distribution:\n%s", df["szbin"].value_counts().to_string())

    # Optimize memory before writing
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
