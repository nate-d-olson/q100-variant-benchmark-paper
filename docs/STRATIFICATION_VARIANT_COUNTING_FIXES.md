# Stratification Variant Counting and Coverage Metrics Fixes

## Summary

Fixed critical issues in the variant counting and benchmark region coverage analysis for difficult stratifications (TR, HP, SD, MAP, etc.). The rules now properly execute for all benchmark sets with correct data aggregation and metric computation.

## Issues Fixed

### 1. Lambda Function Calls in download_stratification Rule

**File:** `workflow/rules/downloads.smk` (lines 233-260)

**Problem:**

- `params.url` was using `lambda w: get_stratification_url` (passing function reference instead of calling it)
- `sha256` constraint was using `get_stratification_sha256` directly instead of wrapped in lambda

**Impact:**

- Snakemake couldn't resolve stratification download URLs dynamically
- SHA256 validation couldn't access wildcard values

**Fix:**

```python
# Before (broken):
params:
    url=lambda w: get_stratification_url,
output:
    bed=ensure(
        "resources/stratifications/{ref}_{strat_name}.bed.gz",
        sha256=get_stratification_sha256,
        ...
    ),

# After (fixed):
params:
    url=lambda w: get_stratification_url(w),
output:
    bed=ensure(
        "resources/stratifications/{ref}_{strat_name}.bed.gz",
        sha256=lambda w: get_stratification_sha256(w),
        ...
    ),
```

### 2. Wildcard Constraints for download_stratification Rule

**File:** `workflow/rules/downloads.smk` (lines 233-260)

**Problem:**

- No constraints on `{ref}` and `{strat_name}` wildcards
- Snakemake couldn't unambiguously parse filenames like `CHM13v2.0_TR10kb.bed.gz` (due to dots in ref names)
- Could match as `ref=CHM13v2` + `strat_name=0_TR10kb` or `ref=CHM13v2.0` + `strat_name=TR10kb`

**Impact:**

- Wildcard resolution was ambiguous
- Snakemake might fail to match dependencies correctly
- Rules could end up with incorrect wildcard values

**Fix:**

Added explicit wildcard_constraints to clarify the parsing:

```python
wildcard_constraints:
    ref="GRCh37|GRCh38|CHM13v2.0",
    strat_name="TR|TR10kb|HP|SD|SD10kb|MAP",
```

### 3. MacOS Compatibility (zcat Replacement)

**File:** `workflow/rules/strat_metrics.smk`

**Problem:**

- The pipeline used `zcat` to decompress `.bed.gz` files.
- On MacOS, `zcat` expects files to have `.Z` extension or similar, causing output errors like `gzip: ... .bed.gz.Z: No such file or directory`.

**Impact:**

- Pipeline execution failed on MacOS environments.

**Fix:**

Replaced `zcat` with `gzip -dc` which is portable across Linux and MacOS.

```python
shell:
    """
    gzip -dc {input.bed} | \
    awk '{{sum += $3-$2}} END {{print sum}}' > {output.size}
    """
```

### 4. Robust Variant Table Header Parsing

**File:** `workflow/scripts/count_variants_by_strat.py`

**Problem:**

- `bcftools query` outputs headers with indices (e.g., `[1]CHROM`) when multiple files are queried.
- Some VCF fields like `INFO/STRAT_IDS` might appear with or without the `INFO/` prefix depending on the tool chain.
- The script strictly looked for `INFO/STRAT_IDS`, causing `ValueError: Missing required columns`.

**Impact:**

- Variant counting failed for all stratifications with `ValueError`.

**Fix:**

Implemented a robust column finding function that normalizes headers.

```python
def find_col(candidates):
    # candidates is a list of acceptable names e.g., ['STRAT_IDS', 'INFO/STRAT_IDS']
    for f in fieldnames:
        norm = f.split(']', 1)[1] if ']' in f else f
        if norm in candidates:
            return f
    return None
```

## Affected Rules

The fixes enable the following rules to work correctly:

1. **Stratification Metrics Rules:**
   - `download_stratification` - Downloads GIAB stratification BED files
   - `materialize_stratification` - Links downloaded files for processing
   - `compute_stratification_size` - Calculates total stratification bp
   - `compute_stratification_metrics` - Computes overlap with benchmark regions
   - `aggregate_stratification_metrics` - Combines all metrics into summary table

2. **Variant Counting Rules:**
   - `count_variants_by_stratification` - Counts variants in each stratification from variant table
   - `summarize_variant_counts` - Aggregates variant counts by type
   - `combine_metrics_and_counts` - Joins stratification coverage with variant counts

## Output Files Generated

After these fixes, the following output files are now properly generated:

### Stratification Coverage Metrics

- `results/strat_metrics/{benchmark}/stratification_coverage_table.csv`
  - Columns: strat_name, strat_bp, intersect_bp, pct_of_strat, pct_of_dip

### Variant Counts

- `results/var_counts/{benchmark}/variants_by_stratification.csv`
  - Columns: strat_name, var_type, count
- `results/var_counts/{benchmark}/stratification_summary.csv`
  - Columns: strat_name, total_variants, snp_count, indel_count, etc.

### Combined Metrics (for Analysis)

- `results/var_counts/{benchmark}/stratification_combined_metrics.csv`
  - Combines coverage metrics with variant counts and density calculations

## Testing

Run the following to verify the fixes:

```bash
# Dry-run to verify workflow builds
snakemake -n

# Run the workflow
snakemake --cores 8

# Test specific benchmark
snakemake results/var_counts/v5q_chm13_smvar/stratification_combined_metrics.csv
```

## References

- Development instructions: `.github/instructions/development.instructions.md`
- Snakemake standards: `.github/instructions/snakemake.instructions.md`
- Analysis guide: `docs/DIFFICULT_REGIONS_ANALYSIS.md`
