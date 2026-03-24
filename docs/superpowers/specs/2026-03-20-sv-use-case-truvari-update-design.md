# SV Use-Case Evaluation: Truvari Metrics Update

**Date**: 2026-03-20
**Status**: Approved
**Implementation**: Approach 1 (minimal script updates for quick manuscript completion)

## Goal

Update the structural variant use-case evaluation to:
1. Compare ONT-Sniffles (mapping-based) vs ONT-Verkko (assembly-based) SV callsets
2. Use proper Truvari metrics extraction (refine summary + stratify utility)
3. Maintain existing Quarto notebook compatibility
4. Enable same-day manuscript draft completion

Future work: Full Snakemake integration (TODO.md line 60)

## Requirements

### Callsets
- **ONT-Sniffles**: Mapping-based SV calling (Sniffles2 v2.5.3)
  - Source: `/Users/nolson/Google Drive/.../ont-sniffles/v5hac40x_snifflesPhased_v2.5.3tr.pass.vcf.gz`
- **ONT-Verkko**: Assembly-based SV calling (Verkko + dipcall)
  - Source: `/Users/nolson/Google Drive/.../ont-verkko/twoFCverkko.dipcall.dip.svs.splitmulti.vcf.gz`

### Benchmark
- v5.0q GRCh38 structural variant benchmark
- GIAB v5.0q recommended Truvari parameters: `--pick ac --passonly -r 2000 -C 5000`

### Metrics
1. **Overall performance**: Precision, recall, F1 from `truvari refine` output
2. **Stratified by genomic context**: HP, MAP, SD, TR (using `truvari stratify`)
3. **Stratified by SV type**: DEL, INS, DUP, INV, etc.
4. **Stratified by size bin**: Using `truvari.get_sizebin()`

### Constraints
- Small variant analysis unchanged (rtg vcfeval for Roche NeuSomatic)
- Minimal code changes to existing scripts
- Output CSVs must match current format for Quarto notebook compatibility

## Architecture

### File Structure
```
scripts/
├── run_sv_use_case.sh           # Updated: new callsets, add stratify step
└── extract_sv_metrics.py        # Updated: Truvari API for metrics

results/use_case/
├── callsets/                    # Copied from Google Drive
│   ├── ont-sniffles.vcf.gz
│   └── ont-verkko.vcf.gz
└── stvar/                       # Benchmarking outputs
    ├── ont-sniffles_v5.0q/
    │   ├── summary.json         # truvari refine metrics
    │   ├── tp-base.vcf.gz       # post-refine TPs
    │   ├── fp.vcf.gz            # post-refine FPs
    │   ├── fn.vcf.gz            # post-refine FNs
    │   ├── stratify_HP.txt      # per-context counts
    │   ├── stratify_TR.txt
    │   ├── stratify_SD.txt
    │   └── stratify_MAP.txt
    ├── ont-verkko_v5.0q/
    │   └── (same structure)
    ├── stratified_metrics.csv   # Overall + context metrics
    ├── svtype_metrics.csv       # Metrics by SV type
    └── svtype_size_counts.csv   # TP/FP/FN counts by type + size
```

## Component Design

### 1. Bash Orchestration (`run_sv_use_case.sh`)

**Updated variables**:
```bash
CALLSET_NAMES="ont-sniffles ont-verkko"
SNIFFLES_VCF="${CALLSET_DIR}/ont-sniffles.vcf.gz"
VERKKO_VCF="${CALLSET_DIR}/ont-verkko.vcf.gz"
```

**Workflow per callset**:
1. Copy VCF from Google Drive to `results/use_case/callsets/` (preserve original names, then symlink)
2. Filter ALT="*" variants per GIAB README
3. Run `truvari bench` with GIAB v5.0q params
4. Run `truvari refine` with `--use-original-vcfs --threads 16 --debug`
5. **NEW**: Run `truvari stratify` for each genomic context:
   ```bash
   for context in HP TR SD MAP; do
     truvari stratify \
       --base ${BENCH_VCF} \
       --comp ${filtered_vcf} \
       --bench-dir ${run_dir} \
       --regions ${STRAT_DIR}/GRCh38_${context}.bed.gz \
       --output ${run_dir}/stratify_${context}.txt
   done
   ```

**Key changes**:
- Add Google Drive → local copy step
- Add stratify loop after refine
- Update callset names and file paths
- Remove old Baylor ONT / DRAGEN results

### 2. Metrics Extraction (`extract_sv_metrics.py`)

**Dependencies**:
```python
import json
import pandas as pd
import truvari
from pathlib import Path
```

**Output 1: Stratified Metrics** (`stratified_metrics.csv`)

Columns: `callset, context, tp, fp, fn, recall, precision`

**Overall metrics** (from `summary.json`):
```python
with open(run_dir / "summary.json") as f:
    summary = json.load(f)

precision, recall, f1 = truvari.performance_metrics(
    summary["TP-base"],
    summary["TP-call"],
    summary["FN"],
    summary["FP"]
)
```

**Per-context metrics** (from `stratify_*.txt`):
```python
for context in ["HP", "TR", "SD", "MAP"]:
    strat_df = pd.read_csv(
        run_dir / f"stratify_{context}.txt",
        sep='\t',
        names=['chrom', 'start', 'end', 'tpbase', 'tp', 'fn', 'fp']
    )

    # Sum across all regions in this context
    totals = strat_df[["tpbase", "tp", "fn", "fp"]].sum()

    # Calculate metrics using Truvari
    precision, recall, f1 = truvari.performance_metrics(
        totals["tpbase"],
        totals["tp"],
        totals["fn"],
        totals["fp"]
    )
```

**Output 2: SV Type Metrics** (`svtype_metrics.csv`)

Columns: `callset, svtype, tp, fp, fn, recall, precision`

```python
# Load post-refine VCFs
tp_df = truvari.vcf_to_df(str(run_dir / "tp-base.vcf.gz"))
fp_df = truvari.vcf_to_df(str(run_dir / "fp.vcf.gz"))
fn_df = truvari.vcf_to_df(str(run_dir / "fn.vcf.gz"))

# Group by svtype and calculate metrics
for svtype in set(tp_df["svtype"]) | set(fp_df["svtype"]) | set(fn_df["svtype"]):
    tp_count = len(tp_df[tp_df["svtype"] == svtype])
    fp_count = len(fp_df[fp_df["svtype"] == svtype])
    fn_count = len(fn_df[fn_df["svtype"] == svtype])

    precision, recall, f1 = truvari.performance_metrics(
        tp_count,  # tpbase
        tp_count,  # tp (same for post-refine)
        fn_count,
        fp_count
    )
```

**Output 3: Size-Binned Counts** (`svtype_size_counts.csv`)

Columns: `callset, category, svtype, size_bin, count`

```python
# Add size_bin column using Truvari's get_sizebin()
for df in [tp_df, fp_df, fn_df]:
    df["size_bin"] = df["svlen"].apply(
        lambda x: truvari.get_sizebin(abs(x) if x else 0)
    )

# Build output rows
size_rows = []
for category, df, label in [("TP", tp_df), ("FP", fp_df), ("FN", fn_df)]:
    counts = df.groupby(["svtype", "size_bin"]).size().reset_index(name="count")
    counts["category"] = label
    counts["callset"] = callset
    size_rows.extend(counts.to_dict("records"))
```

**Note**: `truvari.get_sizebin()` returns bracket notation (`"[50,100)"`) rather than the current format (`"50-99bp"`). Quarto notebook may need minor label updates.

### 3. Quarto Notebook Updates (`analysis/use_case_evaluation.qmd`)

**Required changes**:
1. Update callset labels:
   - `"ont-sniffles"` → `"ONT-Sniffles (Sniffles2)"`
   - `"ont-verkko"` → `"ONT-Verkko (dipcall)"`
2. Update color mapping:
   ```r
   callset_colors <- c(
     "ONT-Sniffles (Sniffles2)" = "#2ca02c",
     "ONT-Verkko (dipcall)" = "#1f77b4"
   )
   ```
3. Update size bin levels if needed (depending on `truvari.get_sizebin()` output format)
4. Update results text to reflect assembly-based vs mapping-based comparison

**Optional changes**:
- Update figure titles and axis labels
- Adjust phred scale breaks if performance ranges differ

## Key Design Decisions

### 1. Why Truvari stratify instead of bedtools?

**Current approach** (bedtools intersect):
```bash
bedtools intersect -u -a tp-base.vcf.gz -b GRCh38_HP.bed.gz | wc -l
```

**New approach** (truvari stratify):
```bash
truvari stratify --bench-dir run_dir --regions GRCh38_HP.bed.gz --output stratify_HP.txt
```

**Advantages**:
- Integrated with Truvari workflow (uses same bench directory)
- Outputs counts directly (no VCF parsing needed)
- Consistent with Truvari's variant matching logic
- Follows GIAB documentation examples

### 2. Why truvari.performance_metrics() everywhere?

**Rationale**:
- Handles edge cases (zero denominators, missing values)
- Consistent precision/recall calculation across all stratifications
- Matches metrics reported in `summary.json`
- Used by Truvari bench internally

### 3. Why keep bash + Python structure?

**For quick implementation** (Approach 1):
- Minimal changes to working code
- Familiar structure for debugging
- Can iterate quickly during manuscript revision

**Future Snakemake integration** (Approach 3):
- Will consolidate into proper pipeline rules
- Add dependency tracking and caching
- Conda environment management
- Tracked as TODO item

## Testing Strategy

1. **Dry run**: Check file paths and argument syntax
2. **Single callset**: Run ont-sniffles only, verify outputs
3. **Full run**: Both callsets, check CSV format matches expected columns
4. **Quarto render**: Verify notebook renders without errors
5. **Visual inspection**: Check figures for reasonable performance values

## Rollback Plan

If issues arise:
1. Git checkout original scripts from HEAD
2. Keep existing DRAGEN/Baylor results in separate directory
3. Fall back to current approach with manual metric calculation

## Future Work (TODO)

- [ ] Integrate SV use-case evaluation into Snakemake pipeline (Approach 3)
- [ ] Add config entries for callsets, parameters, stratifications
- [ ] Create `workflow/rules/use_case_evaluation.smk`
- [ ] Add Snakemake rules for download, bench, refine, stratify, extract
- [ ] Update documentation in `workflow/README.md`
