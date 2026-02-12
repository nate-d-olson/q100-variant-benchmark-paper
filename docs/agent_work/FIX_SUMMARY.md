# Complete Fix Summary: dip_bed KeyError

## Error Fixed

```
InputFunctionException in rule download_benchmark_dip_bed
KeyError: 'dip_bed'
Wildcards: benchmark=v421_grch38_smvar
```

## Root Cause

Snakemake evaluates lambda functions during DAG construction/parsing to understand rule structure. The `download_benchmark_dip_bed` rule's lambdas tried to access `config["benchmarksets"]["v421_grch38_smvar"]["dip_bed"]`, which doesn't exist.

## Complete Solution (3 Parts)

### 1. Filter Input Helper Functions

**File:** [workflow/rules/common.smk:158-189](/Users/nolson/active/q100-papers/q100-variant-benchmark/workflow/rules/common.smk#L158-L189)

```python
def get_strat_metrics_inputs(wildcards):
    return [
        f"results/strat_metrics/{benchmark}/stratification_coverage_table.csv"
        for benchmark, conf in config["benchmarksets"].items()
        if "dip_bed" in conf  # Only v5q benchmarks
    ]
```

**Purpose:** Prevents requesting stratification metrics for benchmarks without dip_bed.

### 2. Defensive Lambda Functions

**File:** [workflow/rules/downloads.smk:115-118](/Users/nolson/active/q100-papers/q100-variant-benchmark/workflow/rules/downloads.smk#L115-L118)

```python
output:
    bed=ensure(
        "resources/benchmarksets/{benchmark}_dip.bed",
        non_empty=True,
        sha256=lambda w: config["benchmarksets"][w.benchmark].get("dip_bed", {}).get("sha256", ""),
    ),
params:
    url=lambda w: config["benchmarksets"][w.benchmark].get("dip_bed", {}).get("url", ""),
```

**Purpose:** Prevents KeyError when Snakemake evaluates lambdas during DAG construction.

### 3. Wildcard Constraint

**File:** [workflow/rules/downloads.smk:119-120](/Users/nolson/active/q100-papers/q100-variant-benchmark/workflow/rules/downloads.smk#L119-L120)

```python
wildcard_constraints:
    benchmark="|".join([b for b, conf in config["benchmarksets"].items() if "dip_bed" in conf]),
```

**Purpose:** Explicitly limits which benchmarks can match this rule.

## Why All 3 Parts Are Needed

| Part | Prevents | When Applied |
|------|----------|--------------|
| 1. Filter helpers | Requesting stratification outputs for non-v5q | During `rule all` input evaluation |
| 2. Defensive lambdas | KeyError during lambda evaluation | During rule parsing/DAG construction |
| 3. Wildcard constraint | Rule matching incorrect benchmarks | During wildcard resolution |

Without part 2, the error occurs during Snakemake parsing before parts 1 and 3 can take effect.

## Affected Benchmarks

**Processed (with dip_bed):**

- v5q_chm13_smvar ✅
- v5q_chm13_stvar ✅
- v5q_grch37_smvar ✅
- v5q_grch37_stvar ✅
- v5q_grch38_smvar ✅
- v5q_grch38_stvar ✅

**Skipped (no dip_bed):**

- v421_grch38_smvar ⏭️
- v06_grch37_stvar ⏭️

## Testing

```bash
# Should now run without errors
snakemake --cores 8

# Or dry-run to verify
snakemake --dry-run
```

## Files Modified

1. [workflow/rules/common.smk](/Users/nolson/active/q100-papers/q100-variant-benchmark/workflow/rules/common.smk) - Lines 158-189 (filter helpers)
2. [workflow/rules/downloads.smk](/Users/nolson/active/q100-papers/q100-variant-benchmark/workflow/rules/downloads.smk) - Lines 115-120 (defensive lambdas + constraint)

## Documentation

- [docs/FIX_DIP_BED_ERROR.md](/Users/nolson/active/q100-papers/q100-variant-benchmark/docs/FIX_DIP_BED_ERROR.md) - Detailed explanation with examples
- [docs/FIX_SUMMARY.md](/Users/nolson/active/q100-papers/q100-variant-benchmark/docs/FIX_SUMMARY.md) - This file

## Date

2026-01-13
