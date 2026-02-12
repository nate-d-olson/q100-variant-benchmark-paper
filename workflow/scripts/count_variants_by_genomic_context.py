#!/usr/bin/env python3
# AI Disclosure: This script was developed with assistance from Claude (Anthropic).
"""
Count variants by genomic context from Parquet variant table.

Reads the Parquet variant table (already size-filtered and type-classified)
and produces two count views in a single output Parquet:
1. Counts by (context_name, var_type) — for total/per-type counts
2. Counts by (context_name, var_type, szbin) — for size distribution analysis
"""

import logging

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


def count_variants_by_genomic_context(
    parquet_path: str, output_path: str, log_path: str
) -> None:
    """Count variants by genomic context and variant type from Parquet.

    Args:
        parquet_path: Path to variants.parquet
        output_path: Path to output Parquet file
        log_path: Path to log file
    """
    logging.basicConfig(
        filename=log_path,
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )

    logging.info("Counting variants by genomic context")
    logging.info("Input: %s", parquet_path)
    logging.info("Output: %s", output_path)

    # Read only needed columns from Parquet
    df = pd.read_parquet(parquet_path, columns=["context_ids", "var_type", "szbin"])
    logging.info("Loaded %d variants", len(df))

    # Filter to variants with genomic context annotations
    # Empty strings are already converted to NA in generate_variant_parquet.py
    df = df[df["context_ids"].notna()].copy()
    logging.info("%d variants have genomic context annotations", len(df))

    # Explode comma-separated context_ids into separate rows
    df["context_ids"] = df["context_ids"].str.split(",")
    df_exploded = df.explode("context_ids")
    df_exploded["context_ids"] = df_exploded["context_ids"].str.strip()
    # Remove empty strings from split
    df_exploded = df_exploded[df_exploded["context_ids"] != ""]
    logging.info("Exploded to %d (context, variant) pairs", len(df_exploded))

    # Count by (context_name, var_type, szbin)
    counts = (
        df_exploded.groupby(["context_ids", "var_type", "szbin"], observed=True)
        .size()
        .reset_index(name="count")
        .rename(columns={"context_ids": "context_name"})
        .sort_values(["context_name", "var_type", "szbin"])
    )

    logging.info("Result: %d rows", len(counts))
    logging.info("Contexts: %s", sorted(counts["context_name"].unique()))
    logging.info("Variant types: %s", sorted(counts["var_type"].unique()))

    # Write Parquet
    table = pa.Table.from_pandas(counts, preserve_index=False)
    pq.write_table(
        table,
        output_path,
        compression="zstd",
        compression_level=3,
    )
    logging.info("Wrote %d rows to %s", len(counts), output_path)


if __name__ == "__main__":
    count_variants_by_genomic_context(
        parquet_path=snakemake.input.parquet,
        output_path=snakemake.output.parquet,
        log_path=snakemake.log[0],
    )
