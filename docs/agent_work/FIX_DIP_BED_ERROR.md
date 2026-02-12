# Fix: KeyError for dip_bed in Older Benchmarks

## Problem

The pipeline was failing with the following error when processing older benchmarks (v421, v06):

```
InputFunctionException in rule download_benchmark_dip_bed
KeyError: 'dip_bed'
Wildcards: benchmark=v421_grch38_smvar
```

## Root Cause

**Two related issues:**

1. The stratification metrics and variant counting modules tried to process ALL benchmarks, but only v5q benchmarks have `dip.bed` files configured
2. The `download_benchmark_dip_bed` rule had no wildcard constraints, so Snakemake considered it a match for ANY benchmark pattern, causing lambda evaluation errors

**Dependency chain:**

```
var_counts (combine_metrics_and_counts)
    ↓ depends on
strat_metrics (aggregate_stratification_metrics)
    ↓ depends on
dip.bed file (compute_stratification_metrics needs it)
    ↓ would trigger
download_benchmark_dip_bed (for ANY benchmark, causing KeyError)
```

**Configuration reality:**

- ✅ v5q benchmarks have `dip_bed` configured in config.yaml
- ❌ v421 and v06 benchmarks do NOT have `dip_bed` configured

## Solution

Two fixes were required:

### Fix 1: Filter Helper Functions

Updated helper functions in [workflow/rules/common.smk](/Users/nolson/active/q100-papers/q100-variant-benchmark/workflow/rules/common.smk) to only request stratification metrics for benchmarks with `dip_bed`:

**Before (Broken):**

```python
def get_strat_metrics_inputs(wildcards):
    """Generate list for all benchmarks."""
    return [
        f"results/strat_metrics/{benchmark}/stratification_coverage_table.csv"
        for benchmark in config["benchmarksets"]
    ]

def get_var_counts_inputs(wildcards):
    """Generate list for all benchmarks."""
    return [
        f"results/var_counts/{benchmark}/stratification_combined_metrics.csv"
        for benchmark in config["benchmarksets"]
    ]
```

**After (Fixed):**

```python
def get_strat_metrics_inputs(wildcards):
    """Generate list only for benchmarks with dip_bed."""
    return [
        f"results/strat_metrics/{benchmark}/stratification_coverage_table.csv"
        for benchmark, conf in config["benchmarksets"].items()
        if "dip_bed" in conf  # ← Filter condition added
    ]

def get_var_counts_inputs(wildcards):
    """Generate list only for benchmarks with dip_bed."""
    return [
        f"results/var_counts/{benchmark}/stratification_combined_metrics.csv"
        for benchmark, conf in config["benchmarksets"].items()
        if "dip_bed" in conf  # ← Filter condition added
    ]
```

### Fix 2: Make Download Rule Defensive + Add Wildcard Constraint

Modified `download_benchmark_dip_bed` in [workflow/rules/downloads.smk](/Users/nolson/active/q100-papers/q100-variant-benchmark/workflow/rules/downloads.smk) with two changes:

**Before (Broken):**

```python
rule download_benchmark_dip_bed:
    output:
        bed=ensure(
            "resources/benchmarksets/{benchmark}_dip.bed",
            non_empty=True,
            sha256=lambda w: config["benchmarksets"][w.benchmark]["dip_bed"]["sha256"],
        ),
    params:
        url=lambda w: config["benchmarksets"][w.benchmark]["dip_bed"]["url"],
```

**After (Fixed):**

```python
rule download_benchmark_dip_bed:
    """Download dipcall BED file - only for benchmarks with dip_bed configured."""
    output:
        bed=ensure(
            "resources/benchmarksets/{benchmark}_dip.bed",
            non_empty=True,
            # Defensive: use .get() to avoid KeyError during DAG construction
            sha256=lambda w: config["benchmarksets"][w.benchmark].get("dip_bed", {}).get("sha256", ""),
        ),
    params:
        # Defensive: use .get() to avoid KeyError during DAG construction
        url=lambda w: config["benchmarksets"][w.benchmark].get("dip_bed", {}).get("url", ""),
    wildcard_constraints:
        # Only allow benchmarks that have dip_bed configured
        benchmark="|".join([
            b for b, conf in config["benchmarksets"].items()
            if "dip_bed" in conf
        ]),
```

**Why both changes are necessary:**

1. **Defensive lambdas** (`.get()` with defaults): Prevents KeyError when Snakemake evaluates lambdas during DAG construction to understand rule structure
2. **Wildcard constraint**: Ensures the rule only matches v5q benchmarks, so the defensive defaults are never actually used in practice

Without defensive lambdas, Snakemake throws KeyError during parsing even before wildcard constraints are applied. Without the wildcard constraint, the rule could theoretically match non-v5q benchmarks (though it wouldn't be requested due to Fix 1).

## Impact

### Benchmarks Processed (After Fix)

**Stratification metrics & variant counts:**

- ✅ v5q_chm13_smvar
- ✅ v5q_chm13_stvar
- ✅ v5q_grch37_smvar
- ✅ v5q_grch37_stvar
- ✅ v5q_grch38_smvar
- ✅ v5q_grch38_stvar

**Excluded from stratification analysis:**

- ❌ v421_grch38_smvar (no dip_bed)
- ❌ v06_grch37_stvar (no dip_bed)

**Still processed for other outputs:**

- ✅ ALL benchmarks still get variant tables generated (unchanged)
- ✅ ALL benchmarks with exclusions still get exclusion metrics (unchanged)

## Why Both Fixes are Needed

1. **Fix 1 (filtering)** - Prevents `rule all` from requesting stratification metrics for v421/v06
2. **Fix 2 (constraint)** - Prevents Snakemake from even considering the download rule as a match for v421/v06

Without Fix 2, even though stratification metrics aren't requested, Snakemake might still evaluate the download rule during DAG construction for other purposes, causing the KeyError.

## Pattern for Future Rules

When creating rules that only apply to certain benchmarks based on optional config keys:

```python
rule some_conditional_rule:
    output:
        "path/{benchmark}_output.txt"
    params:
        # Use .get() with defaults to prevent KeyError during DAG construction
        config_val=lambda w: config["benchmarksets"][w.benchmark].get("some_key", {}).get("value", "")
    wildcard_constraints:
        # Limit which benchmarks can match this rule
        benchmark="|".join([
            b for b, conf in config["benchmarksets"].items()
            if "some_key" in conf
        ])
    shell:
        "process {output}"
```

**Key principles:**

1. **Use defensive lambdas** with `.get()` and default values to prevent errors during Snakemake parsing
2. **Add wildcard constraints** to explicitly limit which wildcard values can match the rule
3. The combination ensures both parse-time safety and runtime correctness

## Verification

The fix ensures:

1. ✅ No KeyError for benchmarks without dip_bed
2. ✅ v5q benchmarks still get full stratification analysis
3. ✅ Older benchmarks still get variant tables and other analyses
4. ✅ Pipeline completes successfully for all configured benchmarks

## Related Files Modified

- [workflow/rules/common.smk:158-189](/Users/nolson/active/q100-papers/q100-variant-benchmark/workflow/rules/common.smk#L158-L189) - Filter helper functions
- [workflow/rules/downloads.smk:119-120](/Users/nolson/active/q100-papers/q100-variant-benchmark/workflow/rules/downloads.smk#L119-L120) - Add wildcard constraint

## Date Fixed

2026-01-13
