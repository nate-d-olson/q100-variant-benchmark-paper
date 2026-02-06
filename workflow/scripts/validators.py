#!/usr/bin/env python3
"""
Data validation utilities for Q100 variant benchmark pipeline.

Provides validation functions for VCF, BED, and TSV file formats to detect
corruption, format violations, and data quality issues early in the pipeline.
"""

from pathlib import Path
from typing import Optional, List, Set
import gzip
from logging_config import setup_logger
from exceptions import ValidationError, DataFormatError

logger = setup_logger(__name__)


def validate_file_exists(file_path: Path, file_type: str = "file") -> None:
    """
    Validate that a file exists and is readable.

    Args:
        file_path: Path to file to validate
        file_type: Description of file type for error messages

    Raises:
        ValidationError: If file doesn't exist or isn't readable
    """
    if not file_path.exists():
        raise ValidationError(
            f"{file_type} not found",
            file_path=file_path,
            expected="Existing file",
            actual="File does not exist",
        )

    if not file_path.is_file():
        raise ValidationError(
            f"{file_type} is not a regular file",
            file_path=file_path,
            expected="Regular file",
            actual=f"{'Directory' if file_path.is_dir() else 'Special file'}",
        )

    if not file_path.stat().st_size > 0:
        raise ValidationError(
            f"{file_type} is empty",
            file_path=file_path,
            expected="Non-empty file",
            actual="0 bytes",
        )


def open_maybe_gzip(file_path: Path, mode: str = "rt"):
    """
    Open a file, automatically handling gzip compression.

    Args:
        file_path: Path to file (can be .gz or uncompressed)
        mode: File open mode (default: 'rt' for text read)

    Returns:
        File handle (gzip or regular)
    """
    if str(file_path).endswith(".gz"):
        return gzip.open(file_path, mode)
    return open(file_path, mode)


def validate_bed_format(
    file_path: Path, check_sorted: bool = True, max_lines_to_check: int = 1000
) -> dict:
    """
    Validate BED file format and return basic statistics.

    Args:
        file_path: Path to BED file
        check_sorted: Whether to verify file is sorted by position
        max_lines_to_check: Maximum number of lines to validate

    Returns:
        Dictionary with validation statistics

    Raises:
        DataFormatError: If BED format is invalid
        ValidationError: If file structure is wrong
    """
    logger.info(f"Validating BED file: {file_path}")

    validate_file_exists(file_path, "BED file")

    line_count = 0
    interval_count = 0
    chromosomes: Set[str] = set()
    prev_chrom: Optional[str] = None
    prev_start: Optional[int] = None
    errors: List[str] = []

    try:
        with open_maybe_gzip(file_path) as f:
            for line_num, line in enumerate(f, start=1):
                line = line.strip()

                # Skip empty lines and comments
                if not line or line.startswith("#") or line.startswith("track"):
                    continue

                line_count += 1
                fields = line.split("\t")

                # BED requires at least 3 columns (chrom, start, end)
                if len(fields) < 3:
                    errors.append(
                        f"Line {line_num}: Expected at least 3 columns, found {len(fields)}"
                    )
                    if len(errors) >= 5:  # Stop after 5 errors
                        break
                    continue

                chrom, start_str, end_str = fields[0:3]
                chromosomes.add(chrom)

                # Validate coordinates are integers
                try:
                    start = int(start_str)
                    end = int(end_str)
                except ValueError:
                    errors.append(
                        f"Line {line_num}: Start ({start_str}) or end ({end_str}) is not an integer"
                    )
                    if len(errors) >= 5:
                        break
                    continue

                # Validate start < end
                if start >= end:
                    errors.append(f"Line {line_num}: Start ({start}) >= end ({end})")
                    if len(errors) >= 5:
                        break
                    continue

                # Check sorted order if requested
                if check_sorted and prev_chrom is not None:
                    if chrom == prev_chrom and start < prev_start:
                        errors.append(
                            f"Line {line_num}: File not sorted - position {start} comes after {prev_start} on {chrom}"
                        )

                prev_chrom = chrom
                prev_start = start
                interval_count += 1

                if line_count >= max_lines_to_check:
                    logger.info(
                        f"Checked {max_lines_to_check} lines (validation limit)"
                    )
                    break

    except Exception as e:
        raise DataFormatError(
            f"Failed to parse BED file: {e}", file_path=file_path, format_type="BED"
        ) from e

    if errors:
        error_summary = "\n  ".join(errors)
        raise DataFormatError(
            f"BED format validation failed:\n  {error_summary}",
            file_path=file_path,
            format_type="BED",
        )

    stats = {
        "file": str(file_path),
        "lines_checked": line_count,
        "intervals": interval_count,
        "chromosomes": len(chromosomes),
        "chrom_names": sorted(chromosomes),
    }

    logger.info(
        f"BED validation passed: {interval_count} intervals, "
        f"{len(chromosomes)} chromosomes"
    )

    return stats


def validate_vcf_header(file_path: Path) -> dict:
    """
    Validate VCF file has proper header structure.

    Args:
        file_path: Path to VCF file

    Returns:
        Dictionary with header statistics

    Raises:
        DataFormatError: If VCF header is invalid
        ValidationError: If file structure is wrong
    """
    logger.info(f"Validating VCF header: {file_path}")

    validate_file_exists(file_path, "VCF file")

    has_fileformat = False
    has_column_header = False
    header_lines = 0
    required_columns = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]

    try:
        with open_maybe_gzip(file_path) as f:
            for line_num, line in enumerate(f, start=1):
                line = line.strip()

                if line.startswith("##"):
                    header_lines += 1
                    if line.startswith("##fileformat=VCF"):
                        has_fileformat = True

                elif line.startswith("#CHROM"):
                    has_column_header = True
                    columns = line.split("\t")

                    # Validate required columns
                    missing = []
                    for req_col in required_columns:
                        if req_col not in columns:
                            missing.append(req_col)

                    if missing:
                        raise DataFormatError(
                            "VCF column header missing required columns",
                            file_path=file_path,
                            format_type="VCF",
                            required_columns=required_columns,
                            missing_columns=missing,
                        )
                    break

                elif line and not line.startswith("#"):
                    # Reached data without finding header
                    break

                if line_num > 10000:  # Safety limit
                    break

    except DataFormatError:
        raise
    except Exception as e:
        raise DataFormatError(
            f"Failed to parse VCF header: {e}", file_path=file_path, format_type="VCF"
        ) from e

    if not has_fileformat:
        raise DataFormatError(
            "VCF header missing ##fileformat line",
            file_path=file_path,
            format_type="VCF",
        )

    if not has_column_header:
        raise DataFormatError(
            "VCF header missing #CHROM column header line",
            file_path=file_path,
            format_type="VCF",
        )

    stats = {
        "file": str(file_path),
        "header_lines": header_lines,
        "has_fileformat": has_fileformat,
        "has_column_header": has_column_header,
    }

    logger.info(f"VCF header validation passed: {header_lines} header lines")

    return stats


def validate_tsv_columns(
    file_path: Path,
    required_columns: List[str],
    optional_columns: Optional[List[str]] = None,
) -> dict:
    """
    Validate TSV file has required column headers.

    Args:
        file_path: Path to TSV file
        required_columns: List of required column names
        optional_columns: List of optional column names

    Returns:
        Dictionary with column statistics

    Raises:
        DataFormatError: If required columns are missing
        ValidationError: If file structure is wrong
    """
    logger.info(f"Validating TSV columns: {file_path}")

    validate_file_exists(file_path, "TSV file")

    try:
        with open(file_path) as f:
            header_line = f.readline().strip()

            if not header_line:
                raise DataFormatError(
                    "TSV file has no header line",
                    file_path=file_path,
                    format_type="TSV",
                    required_columns=required_columns,
                )

            columns = header_line.split("\t")

            # Check for required columns
            missing = [col for col in required_columns if col not in columns]

            if missing:
                raise DataFormatError(
                    "TSV file missing required columns",
                    file_path=file_path,
                    format_type="TSV",
                    required_columns=required_columns,
                    missing_columns=missing,
                )

            # Identify which optional columns are present
            present_optional = []
            if optional_columns:
                present_optional = [col for col in optional_columns if col in columns]

    except DataFormatError:
        raise
    except Exception as e:
        raise DataFormatError(
            f"Failed to parse TSV header: {e}", file_path=file_path, format_type="TSV"
        ) from e

    stats = {
        "file": str(file_path),
        "total_columns": len(columns),
        "required_present": len(required_columns),
        "optional_present": len(present_optional) if optional_columns else 0,
        "column_names": columns,
    }

    logger.info(
        f"TSV validation passed: {len(columns)} columns, "
        f"{len(required_columns)} required found"
    )

    return stats
