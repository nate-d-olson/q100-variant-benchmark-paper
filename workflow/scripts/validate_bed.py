#!/usr/bin/env python3
"""
Validate BED file format and write validation report.

Snakemake script to validate BED files and generate validation reports.
"""

from pathlib import Path
from validators import validate_file_exists, validate_bed_format
from logging_config import setup_logger

logger = setup_logger(__name__)


def main(snakemake):
    """
    Main validation function called by Snakemake.

    Args:
        snakemake: Snakemake object with input, output, log attributes
    """
    input_bed = Path(snakemake.input.bed)
    output_report = Path(snakemake.output.report)
    output_report.parent.mkdir(parents=True, exist_ok=True)

    logger.info(f"Validating BED: {input_bed}")

    try:
        # Validate file exists
        validate_file_exists(input_bed, "BED file")

        # Validate BED format
        stats = validate_bed_format(
            input_bed,
            check_sorted=True,
            max_lines_to_check=10000,  # Check first 10k lines
        )

        # Write validation report
        with open(output_report, "w") as f:
            f.write("BED Validation Report\n")
            f.write("=" * 60 + "\n")
            f.write(f"File: {input_bed}\n")
            f.write("Status: PASS\n")
            f.write("\n")
            f.write("Statistics:\n")
            for key, value in stats.items():
                if key == "chrom_names":
                    f.write(f"  {key}: {', '.join(value[:10])}")
                    if len(value) > 10:
                        f.write(f" ... ({len(value)} total)")
                    f.write("\n")
                else:
                    f.write(f"  {key}: {value}\n")
            f.write("\n")
            f.write("All validation checks passed.\n")

        logger.info(f"Validation passed: {input_bed}")

    except Exception as e:
        logger.error(f"Validation failed: {e}")

        # Write failure report
        with open(output_report, "w") as f:
            f.write("BED Validation Report\n")
            f.write("=" * 60 + "\n")
            f.write(f"File: {input_bed}\n")
            f.write("Status: FAIL\n")
            f.write("\n")
            f.write(f"Error:\n{e}\n")

        raise


if __name__ == "__main__":
    main(snakemake)
