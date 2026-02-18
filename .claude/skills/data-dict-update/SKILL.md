---
name: data-dict-update
description: Check and update docs/data-dictionary.md after changes to R loading/schema files or pipeline output formats
user-invocable: false
---

# Data Dictionary Update

The data dictionary at `docs/data-dictionary.md` must be kept in sync with the codebase.
This skill is invoked automatically by Claude whenever relevant files are modified.

## Trigger Files

Invoke this skill whenever any of these files are modified:
- `R/data_loading.R` — loading functions, column names, return values
- `R/schemas.R` — Arrow schemas, factor levels, validation rules
- `R/cache.R` — caching behavior changes
- `workflow/scripts/generate_variant_parquet.py` — variant table output columns
- `workflow/scripts/count_variants_by_genomic_context.py` — variant count output columns
- `workflow/scripts/compute_bed_metrics.py` — metrics output columns
- Any `workflow/rules/*.smk` file that changes output file formats

## What to Check

1. **New columns**: Added to a loading function or pipeline script → document in relevant section
2. **Removed columns**: Dropped from a function → remove from data dictionary
3. **Renamed columns**: Renamed anywhere → update both the data dictionary and verify no other code uses the old name
4. **Factor levels changed**: Added/removed values in `schemas.R` → update documented valid values
5. **New cached datasets**: New `load_*` function added → add full dataset documentation section
6. **New pipeline outputs**: New rule output format → add to pipeline outputs section

## Process

1. Read `docs/data-dictionary.md` to understand current state
2. Read the modified file(s) to understand what changed
3. Identify which sections of the data dictionary are affected
4. Propose specific updates — show exact text changes
5. Apply updates after user confirmation (or directly if changes are unambiguous)
