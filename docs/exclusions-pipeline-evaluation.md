# Exclusions Pipeline Streamlining Evaluation

**Date:** 2026-02-14  
**Related Issue:** Evaluate streamlining exclusions pipeline metrics (compute_bed_metrics.py)  
**Context:** PR #37 streamlined genomic context metrics; evaluate if exclusions should follow

## Executive Summary

**Recommendation:** **Do NOT streamline the exclusions pipeline** to match genomic context refactor.

**Key Finding:** The exclusions pipeline has a fundamentally different architecture than genomic context metrics. The current implementation is already optimized and any changes would add complexity without clear benefits.

## 1. Current Pipeline Architecture Comparison

### Genomic Context Pipeline (After PR #37 Refactor)

```
Input: resources/stratifications/{ref}_{context}.bed.gz
  ↓ (genomic_context_coverage rule - bedtools coverage)
Intermediate: results/genomic_context/{benchmark}/coverage/{context}_cov.bed
  ↓ (compute_genomic_context_coverage_table rule - compute_coverage_table.py)
Output: results/genomic_context/{benchmark}/genomic_context_coverage_table.csv
```

**Characteristics:**
- Two-step process: (1) bedtools coverage → (2) summarize coverage files
- Single aggregation rule processes ALL contexts in one pass
- Uses bedtools coverage output format (7 columns with per-interval statistics)
- Input function `get_genomic_context_cov_beds()` returns list of coverage BED paths
- Script `compute_coverage_table.py` reads pre-computed coverage data

### Exclusions Pipeline (Current)

```
Input: resources/exclusions/{benchmark}/{exclusion}_{idx}.bed
  ↓ (materialize_exclusion rule - bedtools sort/merge)
Intermediate: results/exclusions/{benchmark}/{exclusion}.bed
  ↓ (compute_exclusion_metrics rule - compute_bed_metrics.py)
Output: results/exclusions/{benchmark}/coverage/{exclusion}.tsv
  ↓ (compute_exclusion_impact rule - count_exclusion_variants.py)
Final: results/exclusions/{benchmark}/exclusion_impact.csv
```

**Characteristics:**
- Three-step process: (1) materialize → (2) compute metrics → (3) combine with variant counts
- `compute_bed_metrics.py` computes metrics directly from BED files (not pre-computed coverage)
- Per-exclusion TSV files are intermediate outputs consumed by `compute_exclusion_impact`
- `compute_exclusion_impact` also reads variant table and merges BED + variant metrics

## 2. Key Architectural Differences

### Why Genomic Context Was Refactored

1. **Pre-computed coverage files existed:** The `genomic_context_coverage` rule already generated `_cov.bed` files
2. **No variant integration:** Genomic context metrics are purely BED overlap statistics
3. **Eliminated intermediate TSVs:** Replaced 3-rule chain with single aggregation rule
4. **Removed symlinks:** Eliminated `materialize_genomic_context` symlink rule

### Why Exclusions Should NOT Be Refactored

1. **No pre-computed coverage files:** Exclusions compute metrics directly from source BEDs
2. **Variant integration required:** `compute_exclusion_impact` merges BED metrics with variant counts
3. **Intermediate TSVs are functional:** They bridge `compute_exclusion_metrics` and `compute_exclusion_impact`
4. **Materialize is NOT symlinks:** `materialize_exclusion` performs meaningful work (sort/merge pairs)

## 3. Analysis of `compute_bed_metrics.py`

### Current Usage

**Only used by:** `workflow/rules/exclusions.smk` (rule `compute_exclusion_metrics`)

**Function:**
- Computes overlap metrics between two BED files
- Handles gzipped inputs
- Uses bedtools sort, merge, intersect under the hood
- Outputs: `{region_name}\t{size_a}\t{intersect}\t{pct_of_a}\t{pct_of_b}`

### Dead Code Analysis

**Lines 116-118:**
```python
region_name = snakemake.wildcards.get("exclusion") or snakemake.wildcards.get(
    "genomic_context"
)
```

**Status:** This is indeed dead code after PR #37 refactor
- Genomic context metrics no longer use `compute_bed_metrics.py`
- Only `exclusion` wildcard is ever used now
- The fallback `genomic_context` is never triggered

**Action:** Clean up the dead wildcard fallback

### Comparison with `compute_coverage_table.py`

`compute_coverage_table.py` is **NOT** a replacement for `compute_bed_metrics.py`:

| Feature | compute_bed_metrics.py | compute_coverage_table.py |
|---------|------------------------|---------------------------|
| Input | Two BED files (raw) | Pre-computed _cov.bed files |
| Processing | Runs bedtools commands | Summarizes existing coverage |
| Output | Single TSV row | Multi-row CSV table |
| Scope | One region pair | All contexts in one pass |
| Use case | Exclusions | Genomic contexts |

**Verdict:** These scripts serve different purposes and are not interchangeable.

## 4. Evaluation of `materialize_exclusion`

### Current Implementation

```snakemake
rule materialize_exclusion:
    input:
        files=get_exclusion_inputs,
    output:
        bed=ensure("results/exclusions/{benchmark}/{exclusion}.bed", non_empty=True),
    params:
        exclusion_type=get_exclusion_type,
    shell:
        """
        if [ "{params.exclusion_type}" == "single" ]; then
            bedtools sort -i {input.files} > {output.bed}
        else
            cat {input.files} | bedtools sort -i - | bedtools merge -i - > {output.bed}
        fi
        """
```

### Analysis

**Is this a symlink pattern?** **NO**

- This rule performs meaningful data transformation (sort and/or merge)
- It handles two exclusion types: `single` (one file) and `pair` (two files)
- For `pair` type, it concatenates, sorts, and merges the files
- Input function `get_exclusion_inputs()` returns actual file paths from `resources/exclusions/`

**Comparison with removed genomic context symlinks:**

The genomic context refactor removed `materialize_genomic_context` which was purely creating symlinks:
```snakemake
# OLD (removed in PR #37)
rule materialize_genomic_context:
    input: lambda wc: f"resources/stratifications/{wc.ref}_{wc.context}.bed.gz"
    output: "results/genomic_context/{benchmark}/{context}.bed.gz"
    shell: "ln -sf $(realpath {input}) {output}"
```

**Verdict:** `materialize_exclusion` performs legitimate data processing and should NOT be removed.

## 5. Streamlining Opportunity Assessment

### Could compute_exclusion_metrics be merged with compute_exclusion_impact?

**Analysis:**

Current flow:
1. `compute_exclusion_metrics`: BED overlap metrics (one TSV per exclusion)
2. `compute_exclusion_impact`: Reads TSVs + variant table, merges data

Potential consolidated approach:
- Single rule that computes BED metrics AND counts variants in one pass
- Would eliminate intermediate TSV files
- Script would need to:
  - Take all exclusion BEDs as input
  - Compute metrics for each exclusion (loop)
  - Read variant table
  - Merge and output final CSV

**Pros:**
- Eliminates intermediate TSV files (similar to genomic context refactor)
- Single rule execution (potentially faster)

**Cons:**
- **More complex script:** Mixing BED overlap computation with variant counting
- **Reduced modularity:** BED metrics and variant counting are conceptually separate
- **Harder to debug:** Failures would be harder to isolate
- **No performance gain:** BED metrics are fast; variant counting dominates runtime
- **Breaks separation of concerns:** One script doing two different tasks
- **Different conda environments:** BED metrics use bedtools.yaml, variant counting uses truvari.yaml

**Verdict:** Do NOT merge these rules. The current separation is cleaner and more maintainable.

## 6. Consistency Analysis

### Are genomic context and exclusions inconsistent?

**Short answer:** No, they are appropriately different.

**Key differences that justify separate patterns:**

| Aspect | Genomic Context | Exclusions | Justification |
|--------|-----------------|------------|---------------|
| Pre-computed coverage | Yes (_cov.bed files) | No (direct from BED) | Genomic contexts reuse coverage for multiple analyses |
| Variant integration | No (separate pipeline) | Yes (merged in impact) | Exclusion analysis specifically asks "which variants excluded?" |
| Intermediate files | Coverage BEDs | Metric TSVs | Both serve as functional checkpoints |
| Aggregation | Single table rule | Impact rule merges | Reflects different data flows |
| Input sources | Direct to stratifications/ | Materialized exclusions | Exclusions need sort/merge preprocessing |

### What patterns SHOULD be consistent?

✅ **Both pipelines correctly:**
- Use `ensure()` for output validation
- Include proper logging
- Specify resources (mem_mb, threads)
- Use conda environments
- Follow naming conventions

## 7. Recommendations

### 1. Clean Up Dead Code in `compute_bed_metrics.py`

**Action:** Remove the `genomic_context` wildcard fallback (lines 116-118)

**Change:**
```python
# OLD
region_name = snakemake.wildcards.get("exclusion") or snakemake.wildcards.get(
    "genomic_context"
)

# NEW
region_name = snakemake.wildcards.exclusion
```

**Rationale:** After PR #37, only exclusions use this script. Simplify the code.

### 2. Keep `compute_bed_metrics.py`

**Action:** Retain the script as-is (after cleanup)

**Rationale:**
- It's the correct tool for exclusions pipeline
- Not replaceable by `compute_coverage_table.py`
- Well-tested and functional
- Clear purpose and scope

### 3. Do NOT Refactor Exclusions Pipeline

**Action:** Keep the current 3-rule structure

**Rationale:**
- Current architecture is appropriate for its use case
- Merging rules would reduce maintainability
- No performance or complexity benefit
- Different from genomic context in meaningful ways

### 4. Update Documentation

**Action:** Document why exclusions differ from genomic contexts

**Files to update:**
- `workflow/README.md`: Clarify exclusions pipeline purpose
- `CLAUDE.md`: Document the architectural difference
- `docs/architecture.md` (if exists): Add section on pipeline patterns

### 5. Consider Future Optimization (Optional)

If exclusion metrics computation becomes a bottleneck, consider:
- Pre-computing bedtools coverage files for exclusions (similar to genomic contexts)
- Then could use `compute_coverage_table.py` pattern
- **However:** Current performance is likely acceptable; don't optimize prematurely

## 8. Implementation Plan

### Immediate Changes (This PR)

- [ ] Remove dead `genomic_context` wildcard from `compute_bed_metrics.py`
- [ ] Add docstring clarifying exclusions-only usage
- [ ] Update workflow/README.md to document the difference
- [ ] Add this evaluation document to `docs/`

### No Follow-up Required

No additional PRs needed. The exclusions pipeline is working as designed.

## 9. Conclusion

The genomic context refactor (PR #37) was appropriate for that specific pipeline, but the same changes do not apply to the exclusions pipeline. The architectural differences are intentional and justified:

1. **`compute_bed_metrics.py`** is the correct tool for exclusions and should be retained (with dead code cleanup)
2. **`materialize_exclusion`** performs meaningful data processing (not just symlinks) and should be retained
3. **No streamlining needed** - the current exclusions pipeline is well-designed for its purpose

The two pipelines are **appropriately different**, not inconsistent. Attempting to force them into the same pattern would reduce code quality without tangible benefits.
