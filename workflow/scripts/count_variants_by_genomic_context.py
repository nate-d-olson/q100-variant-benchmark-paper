#!/usr/bin/env python3
# AI Disclosure: This script was developed with assistance from Claude (Anthropic).
"""
Count variants by genomic context from Parquet variant table.

Reads the Parquet variant table (with boolean genomic context columns) and
produces variant counts by (context_name, var_type, szbin).

The Parquet table has boolean columns for each genomic context (e.g., HP, TR,
SD, MAP). This script melts those columns into long format for counting.
"""

import logging
from typing import List

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


def count_variants_by_genomic_context(
    parquet_path: str,
    output_path: str,
    context_columns: List[str],
    log_path: str,
) -> None:
    """Count variants by genomic context and variant type from Parquet.

    Args:
        parquet_path: Path to variants.parquet (with boolean context columns)
        output_path: Path to output Parquet file
        context_columns: List of boolean context column names (e.g., ["HP", "TR", "SD"])
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
    logging.info("Context columns: %s", context_columns)

    # Read only needed columns from Parquet
    columns_to_read = ["var_type", "szbin"] + context_columns
    df = pd.read_parquet(parquet_path, columns=columns_to_read)
    logging.info("Loaded %d variants", len(df))

    # Verify context columns exist
    available_contexts = [c for c in context_columns if c in df.columns]
    missing_contexts = [c for c in context_columns if c not in df.columns]
    if missing_contexts:
        logging.warning("Context columns not found in Parquet: %s", missing_contexts)
    logging.info("Available context columns: %s", available_contexts)

    if not available_contexts:
        logging.warning("No context columns found — creating empty output")
        empty_df = pd.DataFrame(
            columns=["context_name", "var_type", "szbin", "count"]
        )
        table = pa.Table.from_pandas(empty_df, preserve_index=False)
        pq.write_table(table, output_path, compression="zstd", compression_level=3)
        return

    # Melt boolean context columns into long format:
    # Each row becomes one (variant, context) pair where the context column is True
    melted = df.melt(
        id_vars=["var_type", "szbin"],
        value_vars=available_contexts,
        var_name="context_name",
        value_name="in_context",
    )
    # Keep only rows where variant overlaps the context
    melted = melted[melted["in_context"]].drop(columns=["in_context"])
    logging.info("Melted to %d (context, variant) pairs", len(melted))

    # Count by (context_name, var_type, szbin)
    counts = (
        melted.groupby(["context_name", "var_type", "szbin"], observed=True)
        .size()
        .reset_index(name="count")
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
        context_columns=list(snakemake.params.context_columns),
        log_path=snakemake.log[0],
    )
