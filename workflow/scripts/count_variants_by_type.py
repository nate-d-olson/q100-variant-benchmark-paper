#!/usr/bin/env python3
"""
Count variants by type from bcftools query output.

Reads variant types from stdin (one per line), counts occurrences,
and outputs long-format TSV with all variant types including zeros.
Auto-detects variant types present in the data for transparency.

Usage:
    bcftools query -f '%TYPE\n' input.vcf | python count_variants_by_type.py benchmark ref context
    bcftools query -f '%SVTYPE\n' input.vcf | python count_variants_by_type.py benchmark ref context
"""

import sys
from collections import Counter
from typing import List, Dict
from logging_config import setup_logger, log_context
from exceptions import ValidationError, ProcessingError

# Initialize logger
logger = setup_logger(__name__)


def main() -> None:
    """Main entry point for variant counting."""
    try:
        if len(sys.argv) != 4:
            raise ValidationError(
                "Incorrect number of arguments",
                expected="3 arguments: benchmark, ref, context",
                actual=f"{len(sys.argv) - 1} arguments provided"
            )

        benchmark = sys.argv[1]
        ref = sys.argv[2]
        context = sys.argv[3]

        logger.info(
            f"Starting variant counting: {log_context(benchmark=benchmark, ref=ref, context=context)}"
        )

        # Read variant types from stdin
        variant_types: List[str] = []
        line_count = 0

        for line in sys.stdin:
            line_count += 1
            line = line.strip()
            if line:  # Skip empty lines
                variant_types.append(line)

        logger.info(
            f"Read variants from stdin: {log_context(total_lines=line_count, variants=len(variant_types))}"
        )

        # Count variant types
        type_counts: Dict[str, int] = Counter(variant_types)

        # Get all unique types (auto-detected)
        all_types = sorted(type_counts.keys()) if type_counts else []

        logger.info(
            f"Detected variant types: {log_context(unique_types=len(all_types), types=','.join(all_types) if all_types else 'none')}"
        )

        # If no variants found, output header only with note
        if not all_types:
            logger.warning("No variants found in input data")
            print("benchmark\tref\tcontext\tvariant_type\tcount")
            print(f"{benchmark}\t{ref}\t{context}\tNO_VARIANTS\t0")
            return

        # Output long-format TSV
        print("benchmark\tref\tcontext\tvariant_type\tcount")
        for variant_type in all_types:
            count = type_counts.get(variant_type, 0)
            print(f"{benchmark}\t{ref}\t{context}\t{variant_type}\t{count}")

        logger.info(
            f"Variant counting complete: {log_context(output_rows=len(all_types), total_variants=sum(type_counts.values()))}"
        )

    except ValidationError as e:
        logger.error(f"Validation error: {e}")
        print(
            "Usage: count_variants_by_type.py <benchmark> <ref> <context>",
            file=sys.stderr,
        )
        sys.exit(1)
    except Exception as e:
        logger.error(f"Unexpected error during variant counting: {e}", exc_info=True)
        raise ProcessingError(
            "Failed to count variants",
            operation="variant type counting",
            context={"benchmark": benchmark if 'benchmark' in locals() else "unknown"}
        ) from e


if __name__ == "__main__":
    main()
