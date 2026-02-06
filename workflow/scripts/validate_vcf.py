#!/usr/bin/env python3
"""
Validate VCF file format and write validation report.

Snakemake script to validate VCF files and generate validation reports.
"""

from pathlib import Path
from validators import validate_file_exists, validate_vcf_header
from logging_config import setup_logger


logger = setup_logger(__name__)


def main(snakemake):
    """
    Main validation function called by Snakemake.

    Args:
        snakemake: Snakemake object with input, output, log attributes
    """
    input_vcf = Path(snakemake.input.vcf)
    output_report = Path(snakemake.output.report)
    output_report.parent.mkdir(parents=True, exist_ok=True)

    logger.info(f"Validating VCF: {input_vcf}")

    try:
        # Validate file exists
        validate_file_exists(input_vcf, "VCF file")

        # Validate VCF header structure
        stats = validate_vcf_header(input_vcf)

        # Write validation report
        with open(output_report, "w") as f:
            f.write("VCF Validation Report\n")
            f.write("=" * 60 + "\n")
            f.write(f"File: {input_vcf}\n")
            f.write(f"Status: PASS\n")
            f.write("\n")
            f.write("Statistics:\n")
            for key, value in stats.items():
                f.write(f"  {key}: {value}\n")
            f.write("\n")
            f.write("All validation checks passed.\n")

        logger.info(f"Validation passed: {input_vcf}")

    except Exception as e:
        logger.error(f"Validation failed: {e}")

        # Write failure report
        with open(output_report, "w") as f:
            f.write("VCF Validation Report\n")
            f.write("=" * 60 + "\n")
            f.write(f"File: {input_vcf}\n")
            f.write(f"Status: FAIL\n")
            f.write("\n")
            f.write(f"Error:\n{e}\n")

        raise


if __name__ == "__main__":
    main(snakemake)
