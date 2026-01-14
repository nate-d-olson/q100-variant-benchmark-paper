#!/usr/bin/env python3
"""
Expand STRAT_IDS and REGION_IDS columns into binary flags in TSV.
Reads from stdin, writes to stdout.

Usage:
    python expand_annotations.py --strat-ids ID1 ID2 ... --region-ids ID3 ID4 ...
"""
import argparse
import sys
import csv
from typing import List, Set
from logging_config import setup_logger, log_context
from exceptions import DataFormatError, ProcessingError

# Initialize logger
logger = setup_logger(__name__)


def main() -> None:
    """Main entry point for annotation expansion."""
    try:
        parser = argparse.ArgumentParser(description='Expand annotation IDs into binary columns')
        parser.add_argument('--strat-ids', nargs='*', default=[],
                            help='Stratification IDs to expand')
        parser.add_argument('--region-ids', nargs='*', default=[],
                            help='Region IDs to expand')
        args = parser.parse_args()

        logger.info(
            f"Starting annotation expansion: {log_context(strat_ids=len(args.strat_ids), region_ids=len(args.region_ids))}"
        )

        reader = csv.DictReader(sys.stdin, delimiter='\t')

        # Validate that reader has fieldnames
        if not reader.fieldnames:
            raise DataFormatError(
                "Input TSV has no header row",
                format_type="TSV",
                required_columns=["At least one column header"]
            )

        # Prepare new fieldnames
        fieldnames: List[str] = list(reader.fieldnames)

        logger.info(f"Input columns: {log_context(columns=len(fieldnames), names=','.join(fieldnames[:5]))}")

        # Remove container columns
        removed_cols: List[str] = []
        if 'STRAT_IDS' in fieldnames:
            fieldnames.remove('STRAT_IDS')
            removed_cols.append('STRAT_IDS')
        if 'REGION_IDS' in fieldnames:
            fieldnames.remove('REGION_IDS')
            removed_cols.append('REGION_IDS')

        if removed_cols:
            logger.info(f"Removed container columns: {', '.join(removed_cols)}")

        # Add expansion columns
        # Stratifications: STRAT_{ID}
        strat_cols: List[str] = [f"STRAT_{sid}" for sid in args.strat_ids]
        fieldnames.extend(strat_cols)

        # Regions: {ID} (BMKREGIONS or EXCL_...)
        region_cols: List[str] = args.region_ids
        fieldnames.extend(region_cols)

        logger.info(
            f"Added expansion columns: {log_context(strat_cols=len(strat_cols), region_cols=len(region_cols), total_output_cols=len(fieldnames))}"
        )

        writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()

        row_count = 0
        for row in reader:
            row_count += 1

            # Parse STRAT_IDS
            strat_ids_str = row.pop('STRAT_IDS', '')
            strats: Set[str] = set(strat_ids_str.split(',')) if strat_ids_str else set()
            for sid, col in zip(args.strat_ids, strat_cols):
                row[col] = '1' if sid in strats else '0'

            # Parse REGION_IDS
            region_ids_str = row.pop('REGION_IDS', '')
            regions: Set[str] = set(region_ids_str.split(',')) if region_ids_str else set()
            for rid in args.region_ids:
                row[rid] = '1' if rid in regions else '0'

            writer.writerow(row)

        logger.info(
            f"Annotation expansion complete: {log_context(rows_processed=row_count, output_columns=len(fieldnames))}"
        )

    except DataFormatError as e:
        logger.error(f"Data format error: {e}")
        sys.exit(1)
    except csv.Error as e:
        logger.error(f"CSV parsing error: {e}", exc_info=True)
        raise ProcessingError(
            "Failed to parse TSV input",
            operation="annotation expansion",
            context={"row": row_count if 'row_count' in locals() else 0}
        ) from e
    except Exception as e:
        logger.error(f"Unexpected error during annotation expansion: {e}", exc_info=True)
        raise ProcessingError(
            "Failed to expand annotations",
            operation="annotation expansion"
        ) from e


if __name__ == '__main__':
    main()
