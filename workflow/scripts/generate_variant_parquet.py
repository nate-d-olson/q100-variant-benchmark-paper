#!/usr/bin/env python3
# AI Disclosure: This script was developed with assistance from Claude (Anthropic).
"""
Generate variant Parquet table from VCF with column-wise region annotations.

Reads a VCF (pre-split, pre-svinfo) and annotates each variant against
genomic context and region BED files using interval overlap queries in Python.
Produces a Parquet table with:
- Boolean columns for each genomic context (e.g., HP, TR, SD, MAP)
- Boolean column `in_benchmark` for benchmark region membership
- Boolean columns for each exclusion (e.g., excl_flanks, excl_satellites)
- Correct variant type classification using VariantRecord.var_type()
- Size filtering (smvar <50bp, stvar >=50bp)
- Truvari size bins as strings for R compatibility

This replaces the previous bcftools annotate approach which used comma-separated
INFO fields (CONTEXT_IDS, REGION_IDS). Column-wise boolean annotations simplify
downstream genomic context and exclusion analysis.

Output columns match R/schemas.R variant_table schema:
- bench_version, ref, bench_type (metadata)
- chrom, pos, end, gt (genomic coordinates and genotype)
- var_type, var_size, szbin (variant classification)
- ref_len, alt_len, qual, filter, is_pass (VCF quality info)
- <context_name> columns (boolean, one per genomic context)
- in_benchmark (boolean)
- excl_<name> columns (boolean, one per exclusion)
"""

import bisect
import gzip
import logging
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import truvari


# ============================================================================
# Interval Index for BED overlap queries
# ============================================================================


class IntervalIndex:
    """Efficient interval index using sorted intervals with augmented max-end.

    Supports fast point-in-interval queries using binary search. After adding
    all intervals via add(), call build() once, then query with contains().

    The index uses a prefix-max array of end positions to answer "does any
    interval with start <= pos have end > pos?" in O(log n) time.
    """

    def __init__(self, name: str = ""):
        self.name = name
        self._data: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
        self._max_end: Dict[str, List[int]] = {}
        self._built = False
        self._interval_count = 0

    def add(self, chrom: str, start: int, end: int) -> None:
        """Add an interval [start, end) on the given chromosome."""
        self._data[chrom].append((start, end))
        self._interval_count += 1

    def build(self) -> None:
        """Sort intervals and build prefix-max end arrays for fast queries."""
        for chrom in self._data:
            self._data[chrom].sort()
            # Build prefix-max of end values
            max_ends = []
            current_max = 0
            for _start, end in self._data[chrom]:
                current_max = max(current_max, end)
                max_ends.append(current_max)
            self._max_end[chrom] = max_ends
        self._built = True

    def contains(self, chrom: str, pos: int) -> bool:
        """Check if position (0-based) falls within any stored interval.

        An interval [start, end) contains pos if start <= pos < end.
        Uses binary search + prefix-max for O(log n) query time.
        """
        if not self._built:
            raise RuntimeError("Must call build() before querying")
        if chrom not in self._data:
            return False

        ivls = self._data[chrom]
        max_ends = self._max_end[chrom]

        # Find first index where start > pos using binary search.
        # bisect_right with (pos, inf) finds insertion point after all (pos, *) entries,
        # so all entries at indices [0, hi) have start <= pos.
        hi = bisect.bisect_right(ivls, (pos, float("inf")))

        # Check if any interval in [0, hi) has end > pos (i.e., contains pos).
        # The prefix-max array gives us the maximum end value in [0, hi).
        if hi > 0 and max_ends[hi - 1] > pos:
            return True
        return False

    @property
    def interval_count(self) -> int:
        return self._interval_count

    @property
    def chrom_count(self) -> int:
        return len(self._data)


# ============================================================================
# BED file loading
# ============================================================================


def _open_bed(path: str):
    """Open BED file, handling gzip compression."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def load_bed_to_index(bed_paths: List[str], name: str) -> IntervalIndex:
    """Load one or more BED files into a single IntervalIndex.

    Multiple files for the same region (e.g., pair-type exclusions with
    start/end BEDs) are combined into a single index.

    Args:
        bed_paths: Paths to BED files (plain or gzipped)
        name: Name for this index (used in logging)

    Returns:
        Built IntervalIndex ready for queries
    """
    index = IntervalIndex(name=name)
    for bed_path in bed_paths:
        with _open_bed(bed_path) as f:
            for line in f:
                if line.startswith("#") or line.startswith("track"):
                    continue
                fields = line.strip().split("\t")
                if len(fields) >= 3:
                    chrom = fields[0]
                    start = int(fields[1])
                    end = int(fields[2])
                    index.add(chrom, start, end)
    index.build()
    return index


def parse_bed_specs(bed_specs: List[str]) -> List[Tuple[str, str]]:
    """Parse BED spec strings into (path, id) tuples.

    Args:
        bed_specs: List of "path:ID" strings

    Returns:
        List of (path, id) tuples
    """
    result = []
    for spec in bed_specs:
        if ":" not in spec:
            raise ValueError(f"Invalid bed spec '{spec}'. Expected format: path:ID")
        path, bed_id = spec.rsplit(":", 1)
        result.append((path, bed_id))
    return result


def region_id_to_column(region_id: str) -> str:
    """Convert a region ID to a Parquet column name.

    Args:
        region_id: Region ID from BED spec (e.g., "BMKREGIONS", "EXCL_FLANKS")

    Returns:
        Column name (e.g., "in_benchmark", "excl_flanks")
    """
    if region_id == "BMKREGIONS":
        return "in_benchmark"
    if region_id.startswith("EXCL_"):
        return region_id.lower()
    return region_id.lower()


# ============================================================================
# Benchmark ID parsing
# ============================================================================


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


# ============================================================================
# Main variant table generation
# ============================================================================


def generate_variant_parquet(
    vcf_path: str,
    output_path: str,
    bench_type: str,
    benchmark_id: str,
    context_bed_specs: List[str],
    region_bed_specs: List[str],
    log_path: str,
) -> None:
    """Generate variant Parquet table from VCF with column-wise region annotations.

    Reads a VCF, builds interval indices from BED files, and annotates each
    variant with boolean columns indicating overlap with each genomic context
    and region.

    Args:
        vcf_path: Path to input VCF (pre-annotation, post-split/svinfo)
        output_path: Path to output Parquet file
        bench_type: "smvar" or "stvar"
        benchmark_id: Full benchmark ID (e.g., "v5.0q_GRCh38_smvar")
        context_bed_specs: List of "path:CONTEXT_ID" specs for genomic contexts
        region_bed_specs: List of "path:REGION_ID" specs for regions
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

    # ---- Load BED files into interval indices ----
    logging.info("Loading genomic context BED files...")
    context_specs = parse_bed_specs(context_bed_specs)
    context_indices: Dict[str, IntervalIndex] = {}
    context_columns: List[str] = []  # Ordered list of context column names

    # Group context specs by ID (each context typically has one BED)
    context_beds_by_id: Dict[str, List[str]] = defaultdict(list)
    for bed_path, context_id in context_specs:
        context_beds_by_id[context_id].append(bed_path)

    for context_id, bed_paths in context_beds_by_id.items():
        logging.info(
            "  Loading context BED(s) for %s: %s", context_id, bed_paths
        )
        idx = load_bed_to_index(bed_paths, context_id)
        context_indices[context_id] = idx
        context_columns.append(context_id)
        logging.info(
            "    Loaded %d intervals across %d chromosomes",
            idx.interval_count,
            idx.chrom_count,
        )

    logging.info("Loading region BED files...")
    region_specs = parse_bed_specs(region_bed_specs)
    region_indices: Dict[str, IntervalIndex] = {}
    region_columns: List[str] = []  # Ordered list of region column names

    # Group region specs by column name (pair-type exclusions have multiple BEDs)
    region_beds_by_col: Dict[str, List[str]] = defaultdict(list)
    region_col_order: List[str] = []  # Preserve insertion order
    for bed_path, region_id in region_specs:
        col_name = region_id_to_column(region_id)
        region_beds_by_col[col_name].append(bed_path)
        if col_name not in region_col_order:
            region_col_order.append(col_name)

    for col_name in region_col_order:
        bed_paths = region_beds_by_col[col_name]
        logging.info(
            "  Loading region BED(s) for %s: %s", col_name, bed_paths
        )
        idx = load_bed_to_index(bed_paths, col_name)
        region_indices[col_name] = idx
        region_columns.append(col_name)
        logging.info(
            "    Loaded %d intervals across %d chromosomes",
            idx.interval_count,
            idx.chrom_count,
        )

    logging.info(
        "Annotation columns: contexts=%s, regions=%s", context_columns, region_columns
    )

    # ---- Process VCF ----
    size_threshold = 50

    logging.info("Opening VCF with Truvari VariantFile...")
    with truvari.VariantFile(vcf_path) as vcf:
        samples = list(vcf.header.samples)
        if not samples:
            raise ValueError(f"VCF has no samples: {vcf_path}")
        sample = samples[0]
        logging.info("Processing sample: %s", sample)

        variants: List[Dict] = []
        filtered_size = 0
        filtered_non = 0
        annotation_hits: Dict[str, int] = defaultdict(int)

        logging.info("Processing variants...")
        record_count = 0
        for record in vcf:
            record_count += 1
            if record_count % 500_000 == 0:
                logging.info("  Processed %d records...", record_count)

            vr = truvari.VariantRecord(record)

            # Get variant type (as string from SV enum)
            var_type = vr.var_type().name

            # Skip NON (non-variant) records
            if var_type == "NON":
                filtered_non += 1
                continue

            # Get variant size with sign convention
            abs_size = vr.var_size()
            ref_len = len(vr.get_ref())
            alt = vr.get_alt()
            alt_len = len(alt) if alt is not None else 0

            if var_type == "DEL":
                var_size = -abs_size
            elif var_type == "INS":
                var_size = abs_size
            elif var_type in ["UNK", "INV"]:
                logging.debug(
                    "Reclassifying %s to INS/DEL at %s:%d",
                    var_type,
                    record.chrom,
                    record.pos,
                )
                if alt is None:
                    logging.warning(
                        "Skipping %s variant with missing ALT at %s:%d",
                        var_type,
                        record.chrom,
                        record.pos,
                    )
                    continue
                var_size = alt_len - ref_len
                var_type = "INS" if var_size > 0 else "DEL"
            else:
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

            # Get size bin
            szbin = str(truvari.get_sizebin(var_size))

            # Extract genotype
            gt = truvari.get_gt(vr.gt()).name

            # Extract quality and filter
            qual = float(record.qual) if record.qual is not None else None
            filter_val = (
                ";".join(record.filter.keys()) if record.filter else "PASS"
            )
            is_pass = len(record.filter) == 0 or "PASS" in record.filter

            # ---- Annotate with region overlaps ----
            # Use variant start position (0-based, from pysam) for overlap check
            var_pos = record.pos  # pysam provides 0-based position

            variant_data = {
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
            }

            # Add genomic context boolean columns
            for ctx_name in context_columns:
                hit = context_indices[ctx_name].contains(record.chrom, var_pos)
                variant_data[ctx_name] = hit
                if hit:
                    annotation_hits[ctx_name] += 1

            # Add region boolean columns
            for rgn_col in region_columns:
                hit = region_indices[rgn_col].contains(record.chrom, var_pos)
                variant_data[rgn_col] = hit
                if hit:
                    annotation_hits[rgn_col] += 1

            variants.append(variant_data)

        logging.info("Finished processing %d VCF records", record_count)

    logging.info("Collected %d variants after filtering", len(variants))
    if filtered_non > 0:
        logging.info("Filtered %d NON (non-variant) records", filtered_non)
    if filtered_size > 0:
        logging.info(
            "Filtered %d variants by size (%s threshold=%d)",
            filtered_size,
            bench_type,
            size_threshold,
        )

    # Log annotation hit counts
    logging.info("Annotation overlap counts:")
    for col_name in context_columns + region_columns:
        count = annotation_hits.get(col_name, 0)
        pct = (count / len(variants) * 100) if variants else 0
        logging.info("  %s: %d variants (%.1f%%)", col_name, count, pct)

    # ---- Build DataFrame ----
    all_columns = (
        ["chrom", "pos", "end", "gt", "var_type", "var_size", "szbin",
         "ref_len", "alt_len", "qual", "filter", "is_pass"]
        + context_columns
        + region_columns
    )

    if not variants:
        logging.warning("No variants after filtering - creating empty DataFrame")
        df = pd.DataFrame(columns=all_columns)
    else:
        df = pd.DataFrame(variants)

    # Add benchmark metadata columns at the beginning
    meta = parse_benchmark_id(benchmark_id)
    df.insert(0, "bench_version", meta["bench_version"])
    df.insert(1, "ref", meta["ref"])
    df.insert(2, "bench_type", meta["bench_type"])

    # Ensure boolean columns have correct dtype
    for col in context_columns + region_columns:
        if col in df.columns:
            df[col] = df[col].astype(bool)

    # Log variant type distribution
    if not df.empty:
        logging.info(
            "Variant type distribution:\n%s",
            df["var_type"].value_counts().to_string(),
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
    logging.info("Writing Parquet with %d variants, %d columns...", len(df), len(df.columns))
    logging.info("Columns: %s", list(df.columns))
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
        context_bed_specs=snakemake.params.context_bed_specs,
        region_bed_specs=snakemake.params.region_bed_specs,
        log_path=snakemake.log[0],
    )
