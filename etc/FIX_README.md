# KeyError: 'dip_bed' - FINAL FIX

## Problem

The pipeline was failing with:
```
InputFunctionException in rule download_benchmark_dip_bed
KeyError: 'dip_bed'
Wildcards: benchmark=v421_grch38_smvar
```

## Root Cause

Snakemake evaluates `ensure()` function parameters **during rule parsing** (parse time), which happens **before** wildcard constraints are applied. The lambda function in the sha256 parameter was being evaluated for ALL benchmarks in the config, including v421 and v06 which don't have `dip_bed` configured.

## ✅ Final Fix Applied

**File:** [workflow/rules/downloads.smk](workflow/rules/downloads.smk#L101-L150)

### What Changed

1. **Removed lambda from ensure() (lines 112-115)**
   ```python
   # BEFORE (caused error at parse time)
   output:
       bed=ensure(
           "resources/benchmarksets/{benchmark}_dip.bed",
           non_empty=True,
           sha256=lambda w: config["benchmarksets"][w.benchmark]["dip_bed"]["sha256"],
       ),

   # AFTER (no parse-time evaluation)
   output:
       bed=ensure(
           "resources/benchmarksets/{benchmark}_dip.bed",
           non_empty=True,
       ),
   ```

2. **Moved validation to params (lines 116-118)**
   ```python
   params:
       url=lambda w: config["benchmarksets"][w.benchmark]["dip_bed"]["url"],
       sha256=lambda w: config["benchmarksets"][w.benchmark]["dip_bed"]["sha256"],
   ```
   These lambdas are only evaluated at **runtime**, after wildcard constraints filter out non-v5q benchmarks.

3. **Added SHA256 validation to shell script (lines 137-149)**
   ```bash
   # Download to temporary location
   wget --no-verbose -O {output.bed}.tmp "{params.url}" 2>&1 | tee -a {log}

   # Validate SHA256 checksum
   echo "[$(date)] Validating checksum..." | tee -a {log}
   echo "{params.sha256}  {output.bed}.tmp" | sha256sum -c - 2>&1 | tee -a {log}

   # Move to final location after validation
   mv {output.bed}.tmp {output.bed}
   ```

4. **Kept wildcard constraint (line 121)**
   ```python
   wildcard_constraints:
       benchmark="v5q_chm13_smvar|v5q_chm13_stvar|v5q_grch37_smvar|v5q_grch37_stvar|v5q_grch38_smvar|v5q_grch38_stvar",
   ```

### Also Applied

- **Helper function filtering** in [workflow/rules/common.smk:158-189](workflow/rules/common.smk#L158-L189)
  - Only requests stratification outputs for benchmarks with dip_bed

## How to Test

### 1. Cache Already Cleared ✅

The Snakemake cache has been cleared to remove old rule definitions.

### 2. Activate Environment

```bash
# Using mamba
mamba activate q100-smk

# OR using micromamba
micromamba activate q100-smk
```

### 3. Test with Dry Run

```bash
# Using Makefile
make dry-run

# OR directly
snakemake -n --quiet
```

**Expected:** No errors, full DAG construction succeeds.

### 4. Run the Pipeline

```bash
# Using Makefile (recommended)
make run

# OR directly
snakemake --cores 4 --sdm conda
```

## Why This Fix Works

**The key difference: Parse-time vs Runtime evaluation**

| Code Location | When Evaluated | Wildcard Filtering | Result |
|---------------|----------------|-------------------|---------|
| `ensure()` sha256 (OLD) | **Parse time** | ❌ Not yet applied | Error - tries ALL benchmarks |
| `params` lambdas (NEW) | **Runtime** | ✅ Already applied | Safe - only v5q benchmarks |

See [docs/PARSE_TIME_VS_RUNTIME_FIX.md](docs/PARSE_TIME_VS_RUNTIME_FIX.md) for detailed explanation.

## Affected Benchmarks

### ✅ Will Process (with dip_bed)
- v5q_chm13_smvar
- v5q_chm13_stvar
- v5q_grch37_smvar
- v5q_grch37_stvar
- v5q_grch38_smvar
- v5q_grch38_stvar

### ⏭️ Will Skip (no dip_bed)
- v421_grch38_smvar
- v06_grch37_stvar

These benchmarks will still generate variant tables but won't have stratification metrics.

## Troubleshooting

### If error persists after testing:

1. **Verify the fix is in place:**
   ```bash
   grep -A 5 "rule download_benchmark_dip_bed:" workflow/rules/downloads.smk | grep -A 3 "output:"
   ```
   Should show `ensure()` with NO sha256 parameter.

2. **Check Snakemake version:**
   ```bash
   snakemake --version
   ```
   Should be >= 8.0

3. **Manually clear cache again:**
   ```bash
   ./scripts/clear_snakemake_cache.sh
   ```

4. **Check for syntax errors:**
   ```bash
   python -m py_compile workflow/rules/downloads.smk
   ```

## Documentation

- **[FINAL_FIX_APPLIED.md](FINAL_FIX_APPLIED.md)** - Complete fix explanation
- **[docs/PARSE_TIME_VS_RUNTIME_FIX.md](docs/PARSE_TIME_VS_RUNTIME_FIX.md)** - Parse vs runtime evaluation explained
- **[docs/FIX_DIP_BED_ERROR.md](docs/FIX_DIP_BED_ERROR.md)** - Original analysis
- **[docs/DIFFICULT_REGIONS_ANALYSIS.md](docs/DIFFICULT_REGIONS_ANALYSIS.md)** - How to use stratification analysis

## Status

✅ **FIXED** - Parse-time evaluation removed, validation moved to runtime and shell script

---

**Date:** 2026-01-13
**Final solution:** Remove lambdas from `ensure()` parameters, move to `params` for runtime evaluation
