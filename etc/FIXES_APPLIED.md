# Fixes Applied Checklist

## Error Being Fixed

```
InputFunctionException in rule download_benchmark_dip_bed
KeyError: 'dip_bed'
Wildcards: benchmark=v421_grch38_smvar
```

## ✅ All Fixes Applied

### ✅ Fix 1: Filter Helper Functions
**File:** `workflow/rules/common.smk` lines 158-189

**Verified:**
```python
def get_strat_metrics_inputs(wildcards):
    return [
        f"results/strat_metrics/{benchmark}/stratification_coverage_table.csv"
        for benchmark, conf in config["benchmarksets"].items()
        if "dip_bed" in conf  # ← Filter present
    ]
```

**Status:** ✅ APPLIED

---

### ✅ Fix 2: Defensive Lambda Functions
**File:** `workflow/rules/downloads.smk` lines 115-118

**Verified:**
```python
sha256=lambda w: config["benchmarksets"][w.benchmark].get("dip_bed", {}).get("sha256", ""),
...
url=lambda w: config["benchmarksets"][w.benchmark].get("dip_bed", {}).get("url", ""),
```

**Status:** ✅ APPLIED

---

### ✅ Fix 3: Wildcard Constraint
**File:** `workflow/rules/downloads.smk` lines 119-120

**Verified:**
```python
wildcard_constraints:
    benchmark="|".join([b for b, conf in config["benchmarksets"].items() if "dip_bed" in conf]),
```

**Status:** ✅ APPLIED

---

## Summary

| Fix | File | Lines | Status |
|-----|------|-------|--------|
| 1. Filter helpers | common.smk | 158-189 | ✅ |
| 2. Defensive lambdas | downloads.smk | 115-118 | ✅ |
| 3. Wildcard constraint | downloads.smk | 119-120 | ✅ |

## Expected Behavior

**Benchmarks with stratification analysis (6):**
- v5q_chm13_smvar
- v5q_chm13_stvar
- v5q_grch37_smvar
- v5q_grch37_stvar
- v5q_grch38_smvar
- v5q_grch38_stvar

**Benchmarks without stratification (2):**
- v421_grch38_smvar
- v06_grch37_stvar

All benchmarks still get variant tables and other standard outputs.

## Testing

```bash
# The pipeline should now run without errors
snakemake --cores 8

# Or test with dry-run first
snakemake --dry-run --quiet
```

## Documentation Created

1. ✅ [docs/FIX_DIP_BED_ERROR.md](/Users/nolson/active/q100-papers/q100-variant-benchmark/docs/FIX_DIP_BED_ERROR.md) - Detailed explanation
2. ✅ [docs/FIX_SUMMARY.md](/Users/nolson/active/q100-papers/q100-variant-benchmark/docs/FIX_SUMMARY.md) - Quick reference
3. ✅ [FIXES_APPLIED.md](/Users/nolson/active/q100-papers/q100-variant-benchmark/FIXES_APPLIED.md) - This checklist

## Date Applied

2026-01-13

## Next Steps for User

1. Run the pipeline: `snakemake --cores 8`
2. Verify outputs are generated for all 6 v5q benchmarks
3. Confirm no errors for v421 or v06 benchmarks

---

**All fixes verified and applied!** ✅
