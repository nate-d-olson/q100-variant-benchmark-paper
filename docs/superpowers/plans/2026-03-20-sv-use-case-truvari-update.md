# SV Use-Case Truvari Update Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Update SV use-case evaluation to compare ONT-Sniffles vs ONT-Verkko using proper Truvari metrics extraction

**Architecture:** Two-script workflow (bash orchestration + Python metrics extraction) with Truvari stratify for genomic context metrics and Truvari API for size binning

**Tech Stack:** Bash, Python 3.11, Truvari 5.4.0, pandas, bcftools, Google Drive file access

**Design Spec:** `docs/superpowers/specs/2026-03-20-sv-use-case-truvari-update-design.md`

---

## File Structure

**Modified files:**
- `scripts/run_sv_use_case.sh` - Add callset copying, stratify step, update callset names
- `scripts/extract_sv_metrics.py` - Replace bedtools/bcftools with Truvari API
- `analysis/use_case_evaluation.qmd` - Update callset labels and colors

**New files:**
- None (minimal changes approach)

**Deleted outputs:**
- `results/use_case/stvar/baylor_ont_v5.0q/` - Old Baylor ONT results
- `results/use_case/stvar/dragen_v5.0q/` - Old DRAGEN results
- `results/use_case/callsets/hg002_ONT_BCM_GRCh38.vcf.gz*` - Old callset
- `results/use_case/callsets/DRAGEN_4.4.x_NIST_HG002.sv.vcf.gz*` - Old callset

---

## Task 1: Clean Up Old Results

**Files:**
- Delete: `results/use_case/stvar/baylor_ont_v5.0q/`
- Delete: `results/use_case/stvar/dragen_v5.0q/`
- Delete: `results/use_case/callsets/hg002_ONT_BCM_GRCh38.vcf.gz*`
- Delete: `results/use_case/callsets/DRAGEN_4.4.x_NIST_HG002.sv.vcf.gz*`

- [ ] **Step 1: Remove old benchmark results**

```bash
cd /Users/nolson/Desktop/active/q100-papers/q100-variant-benchmark
rm -rf results/use_case/stvar/baylor_ont_v5.0q results/use_case/stvar/dragen_v5.0q
```

Expected: Directories removed

- [ ] **Step 2: Remove old callsets**

```bash
rm -f results/use_case/callsets/hg002_ONT_BCM_GRCh38.vcf.gz*
rm -f results/use_case/callsets/DRAGEN_4.4.x_NIST_HG002.sv.vcf.gz*
```

Expected: VCF files and indices removed

- [ ] **Step 3: Verify cleanup**

```bash
ls -la results/use_case/stvar/
ls -la results/use_case/callsets/
```

Expected: Only RN_HG002_Illumina_PacBio_Oxford.vcf.gz* remains in callsets

- [ ] **Step 4: Commit cleanup**

```bash
git add -u results/use_case/
git commit -m "chore: remove old SV use-case results for ONT callset update"
```

---

## Task 2: Update Bash Orchestration Script

**Files:**
- Modify: `scripts/run_sv_use_case.sh`

- [ ] **Step 1: Verify PROJ_DIR variable exists (line 7)**

Check that this variable is defined:
```bash
PROJ_DIR="$(cd "$(dirname "$0")/.." && pwd)"
```

Expected: `PROJ_DIR=/Users/nolson/Desktop/active/q100-papers/q100-variant-benchmark`

- [ ] **Step 2: Verify stratification directory exists**

```bash
ls -lh resources/stratifications/GRCh38_{HP,TR,SD,MAP}.bed.gz
```

Expected: Four BED files with .tbi indices

- [ ] **Step 3: Update callset variables (lines 15-17)**

Replace:
```bash
CALLSET_NAMES="baylor_ont dragen"
BAYLOR_ONT_VCF="${CALLSET_DIR}/hg002_ONT_BCM_GRCh38.vcf.gz"
DRAGEN_VCF="${CALLSET_DIR}/DRAGEN_4.4.x_NIST_HG002.sv.vcf.gz"
```

With:
```bash
CALLSET_NAMES="ont-sniffles ont-verkko"
SNIFFLES_SRC="/Users/nolson/Google Drive/Shared drives/BBD_Human_Genomics/HG002-Q100v1.1-Variant-Benchmark-Evaluations/evaluations/stvar/callsets/ont-sniffles/v5hac40x_snifflesPhased_v2.5.3tr.pass.vcf.gz"
VERKKO_SRC="/Users/nolson/Google Drive/Shared drives/BBD_Human_Genomics/HG002-Q100v1.1-Variant-Benchmark-Evaluations/evaluations/stvar/callsets/ont-verkko/twoFCverkko.dipcall.dip.svs.splitmulti.vcf.gz"
SNIFFLES_VCF="${CALLSET_DIR}/ont-sniffles.vcf.gz"
VERKKO_VCF="${CALLSET_DIR}/ont-verkko.vcf.gz"
```

- [ ] **Step 4: Update callset selection case statement (lines 19-22)**

Replace:
```bash
for callset in ${CALLSET_NAMES}; do
    case "${callset}" in
        baylor_ont) vcf="${BAYLOR_ONT_VCF}" ;;
        dragen)     vcf="${DRAGEN_VCF}" ;;
    esac
```

With:
```bash
for callset in ${CALLSET_NAMES}; do
    case "${callset}" in
        ont-sniffles) src_vcf="${SNIFFLES_SRC}"; vcf="${SNIFFLES_VCF}" ;;
        ont-verkko)   src_vcf="${VERKKO_SRC}"; vcf="${VERKKO_VCF}" ;;
    esac
```

- [ ] **Step 5: Add callset copy step (after line 27, before filtering)**

Add after `echo "=== Processing ${callset} ==="`:
```bash
    # Step 0: Copy callset from Google Drive if not present
    if [ ! -f "${vcf}" ]; then
        echo "Copying ${callset} from Google Drive..."
        cp "${src_vcf}" "${vcf}"
        cp "${src_vcf}.tbi" "${vcf}.tbi"
    else
        echo "Callset ${vcf} already exists, skipping copy"
    fi
```

- [ ] **Step 6: Add STRAT_DIR variable (after line 12, after REF definition)**

Add:
```bash
STRAT_DIR="${PROJ_DIR}/resources/stratifications"
```

- [ ] **Step 7: Add stratify step (after truvari refine, before cleanup)**

Add after the `truvari refine` block (around line 57):
```bash
    # Step 5: Run truvari stratify for each genomic context
    echo "Running truvari stratify..."
    for context in HP TR SD MAP; do
        echo "  Stratifying ${context}..."
        truvari stratify \
            --base "${BENCH_VCF}" \
            --comp "${filtered_vcf}" \
            --bench-dir "${run_dir}" \
            --regions "${STRAT_DIR}/GRCh38_${context}.bed.gz" \
            --output "${run_dir}/stratify_${context}.txt"
    done
```

- [ ] **Step 8: Test script syntax**

```bash
bash -n scripts/run_sv_use_case.sh
```

Expected: No syntax errors

- [ ] **Step 9: Commit bash script updates**

```bash
git add scripts/run_sv_use_case.sh
git commit -m "feat: update SV use-case script for ONT callsets and Truvari stratify"
```

---

## Task 3: Update Python Metrics Extraction Script

**Files:**
- Modify: `scripts/extract_sv_metrics.py`

- [ ] **Step 1: Update imports (lines 12-15)**

Replace:
```python
import csv
import subprocess
import sys
from pathlib import Path
```

With:
```python
import csv
import json
import sys
from pathlib import Path
import pandas as pd
import truvari
```

- [ ] **Step 2: Update callset names constant (line 21)**

Replace:
```python
CALLSETS = ["baylor_ont", "dragen"]
```

With:
```python
CALLSETS = ["ont-sniffles", "ont-verkko"]
```

- [ ] **Step 3: Replace get_refine_vcfs function with simpler version (lines 32-42)**

Replace entire function with:
```python
def get_refine_vcfs(run_dir: Path) -> dict[str, Path]:
    """Return paths to post-refine TP/FP/FN VCFs."""
    return {
        "TP": run_dir / "tp-base.vcf.gz",
        "FP": run_dir / "fp.vcf.gz",
        "FN": run_dir / "fn.vcf.gz",
    }
```

- [ ] **Step 4: Remove count_vcf_records and count_intersecting_variants functions (lines 45-67)**

Delete both functions completely - they're replaced by Truvari API usage

- [ ] **Step 5: Keep get_svtype_counts function but mark deprecated (lines 69-82)**

Add comment at top:
```python
# NOTE: Deprecated - kept for reference, not used with Truvari API approach
def get_svtype_counts(vcf_path: Path) -> dict[str, int]:
```

- [ ] **Step 6: Remove get_svtype_size_counts function (lines 85-116)**

Delete entire function - replaced by Truvari API with `get_sizebin()`

- [ ] **Step 7: Remove compute_metrics function (lines 119-123)**

Delete function - replaced by `truvari.performance_metrics()`

- [ ] **Step 8: Rewrite stratified metrics section (lines 134-163)**

Replace entire section with:
```python
    # === 1. Stratified metrics (Overall + per-context) ===
    strat_rows: list[dict[str, str | int | float]] = []

    for callset in CALLSETS:
        run_dir = STVAR_DIR / f"{callset}_v5.0q"

        # Overall metrics from summary.json
        with open(run_dir / "summary.json") as f:
            summary = json.load(f)

        precision, recall, f1 = truvari.performance_metrics(
            summary["TP-base"],
            summary["TP-call"],
            summary["FN"],
            summary["FP"]
        )

        strat_rows.append({
            "callset": callset,
            "context": "Overall",
            "tp": summary["TP-call"],
            "fp": summary["FP"],
            "fn": summary["FN"],
            "recall": round(recall, 4),
            "precision": round(precision, 4),
        })

        # Per-context metrics from stratify text files
        for context in CONTEXTS:
            strat_file = run_dir / f"stratify_{context}.txt"
            strat_df = pd.read_csv(
                strat_file,
                sep='\t',
                names=['chrom', 'start', 'end', 'tpbase', 'tp', 'fn', 'fp']
            )

            # Sum across all regions in this context
            totals = strat_df[["tpbase", "tp", "fn", "fp"]].sum()

            precision, recall, f1 = truvari.performance_metrics(
                int(totals["tpbase"]),
                int(totals["tp"]),
                int(totals["fn"]),
                int(totals["fp"])
            )

            strat_rows.append({
                "callset": callset,
                "context": context,
                "tp": int(totals["tp"]),
                "fp": int(totals["fp"]),
                "fn": int(totals["fn"]),
                "recall": round(recall, 4),
                "precision": round(precision, 4),
            })
```

- [ ] **Step 9: Rewrite SV type metrics section (lines 172-192)**

Replace entire section with:
```python
    # === 2. SV type metrics ===
    svtype_rows: list[dict[str, str | int | float]] = []

    for callset in CALLSETS:
        run_dir = STVAR_DIR / f"{callset}_v5.0q"
        vcfs = get_refine_vcfs(run_dir)

        # Load VCFs via Truvari API
        tp_df = truvari.vcf_to_df(str(vcfs["TP"]))
        fp_df = truvari.vcf_to_df(str(vcfs["FP"]))
        fn_df = truvari.vcf_to_df(str(vcfs["FN"]))

        # Get all SV types present in any category
        all_types = sorted(set(tp_df["svtype"]) | set(fp_df["svtype"]) | set(fn_df["svtype"]))

        for svtype in all_types:
            tp = len(tp_df[tp_df["svtype"] == svtype])
            fp = len(fp_df[fp_df["svtype"] == svtype])
            fn = len(fn_df[fn_df["svtype"] == svtype])

            precision, recall, f1 = truvari.performance_metrics(tp, tp, fn, fp)

            svtype_rows.append({
                "callset": callset,
                "svtype": svtype,
                "tp": tp,
                "fp": fp,
                "fn": fn,
                "recall": round(recall, 4),
                "precision": round(precision, 4),
            })
```

- [ ] **Step 10: Rewrite size-binned counts section (lines 201-220)**

Replace entire section with:
```python
    # === 3. SV type x size bin counts ===
    size_rows: list[dict[str, str | int]] = []

    for callset in CALLSETS:
        run_dir = STVAR_DIR / f"{callset}_v5.0q"
        vcfs = get_refine_vcfs(run_dir)

        for category, vcf_path in vcfs.items():
            # Load VCF and add size_bin column
            df = truvari.vcf_to_df(str(vcf_path))
            df["size_bin"] = df["svlen"].apply(
                lambda x: truvari.get_sizebin(abs(x) if x else 0)
            )

            # Count by (svtype, size_bin)
            counts = df.groupby(["svtype", "size_bin"]).size().reset_index(name="count")
            counts["callset"] = callset
            counts["category"] = category

            size_rows.extend(counts.to_dict("records"))
```

- [ ] **Step 11: Update main execution block (lines 126-end)**

The main execution block should look like this after all updates:
```python
def main() -> None:
    # Verify all run directories exist
    for callset in CALLSETS:
        run_dir = STVAR_DIR / f"{callset}_v5.0q"
        if not run_dir.exists():
            print(f"ERROR: {run_dir} does not exist. Run scripts/run_sv_use_case.sh first.", file=sys.stderr)
            sys.exit(1)

    # === 1. Stratified metrics (Overall + per-context) ===
    strat_rows: list[dict[str, str | int | float]] = []
    [... code from Step 8 ...]

    # Write stratified metrics CSV
    strat_path = STVAR_DIR / "stratified_metrics.csv"
    with open(strat_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["callset", "context", "tp", "fp", "fn", "recall", "precision"])
        writer.writeheader()
        writer.writerows(strat_rows)
    print(f"Wrote {strat_path}")

    # === 2. SV type metrics ===
    svtype_rows: list[dict[str, str | int | float]] = []
    [... code from Step 9 ...]

    # Write SV type metrics CSV
    svtype_path = STVAR_DIR / "svtype_metrics.csv"
    with open(svtype_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["callset", "svtype", "tp", "fp", "fn", "recall", "precision"])
        writer.writeheader()
        writer.writerows(svtype_rows)
    print(f"Wrote {svtype_path}")

    # === 3. SV type x size bin counts ===
    size_rows: list[dict[str, str | int]] = []
    [... code from Step 10 ...]

    # Write size-binned counts CSV
    size_path = STVAR_DIR / "svtype_size_counts.csv"
    with open(size_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["callset", "category", "svtype", "size_bin", "count"])
        writer.writeheader()
        writer.writerows(size_rows)
    print(f"Wrote {size_path}")


if __name__ == "__main__":
    main()
```

**Note**: The `[... code from Step X ...]` placeholders represent the code blocks from Steps 8, 9, and 10.

- [ ] **Step 12: Test Python script syntax**

```bash
python3 -m py_compile scripts/extract_sv_metrics.py
```

Expected: No syntax errors

- [ ] **Step 13: Commit Python script updates**

```bash
git add scripts/extract_sv_metrics.py
git commit -m "feat: use Truvari API for SV metrics extraction"
```

---

## Task 4: Run Updated Pipeline

**Files:**
- Execute: `scripts/run_sv_use_case.sh`
- Execute: `scripts/extract_sv_metrics.py`

- [ ] **Step 1: Activate Truvari conda environment**

```bash
conda activate q100-smk
```

Expected: Environment activated (check with `which truvari`)

- [ ] **Step 2: Run benchmarking script (dry-run check)**

```bash
# Quick path check - don't run full benchmark yet
head -30 scripts/run_sv_use_case.sh | grep -E "CALLSET|VCF"
```

Expected: See new callset paths and names

- [ ] **Step 3: Run full benchmarking pipeline**

```bash
bash scripts/run_sv_use_case.sh 2>&1 | tee results/use_case/stvar/run.log
```

Expected output structure:
```
=== Processing ont-sniffles ===
Copying ont-sniffles from Google Drive...
Filtering ALT=* variants...
Running truvari bench...
Running truvari refine...
Running truvari stratify...
  Stratifying HP...
  Stratifying TR...
  Stratifying SD...
  Stratifying MAP...
=== Done ont-sniffles ===
=== Processing ont-verkko ===
...
All SV use-case benchmarks complete.
Run scripts/extract_sv_metrics.py to generate CSV metrics.
```

**Time estimate:** 60-120 minutes for both callsets (Truvari refine MSA step can be slow)

- [ ] **Step 4: Verify benchmark outputs exist**

```bash
ls -lh results/use_case/stvar/ont-sniffles_v5.0q/
ls -lh results/use_case/stvar/ont-verkko_v5.0q/
```

Expected files per callset:
- `summary.json`
- `tp-base.vcf.gz`, `fp.vcf.gz`, `fn.vcf.gz`
- `stratify_HP.txt`, `stratify_TR.txt`, `stratify_SD.txt`, `stratify_MAP.txt`

- [ ] **Step 5: Run metrics extraction**

```bash
python3 scripts/extract_sv_metrics.py
```

Expected output:
```
Wrote results/use_case/stvar/stratified_metrics.csv
Wrote results/use_case/stvar/svtype_metrics.csv
Wrote results/use_case/stvar/svtype_size_counts.csv
```

- [ ] **Step 6: Verify CSV output structure**

```bash
head -5 results/use_case/stvar/stratified_metrics.csv
head -5 results/use_case/stvar/svtype_metrics.csv
head -5 results/use_case/stvar/svtype_size_counts.csv
```

Expected columns:
- `stratified_metrics.csv`: callset, context, tp, fp, fn, recall, precision
- `svtype_metrics.csv`: callset, svtype, tp, fp, fn, recall, precision
- `svtype_size_counts.csv`: callset, category, svtype, size_bin, count

- [ ] **Step 7: Check for size bin label format**

```bash
cut -d, -f4 results/use_case/stvar/svtype_size_counts.csv | sort -u
```

Expected: `[50,100)`, `[100,300)`, etc. (Truvari bracket notation)

- [ ] **Step 8: Commit generated results**

```bash
git add results/use_case/stvar/*.csv
git commit -m "feat: generate SV use-case metrics with ONT callsets"
```

---

## Task 5: Update Quarto Notebook

**Files:**
- Modify: `analysis/use_case_evaluation.qmd`

- [ ] **Step 1: Update callset color mapping (lines 24-27)**

Replace:
```r
callset_colors <- c(
  "Baylor ONT (Sniffles2)" = "#2ca02c",
  "Illumina DRAGEN" = "#1f77b4"
)
```

With:
```r
callset_colors <- c(
  "ONT-Sniffles (Sniffles2)" = "#2ca02c",
  "ONT-Verkko (dipcall)" = "#1f77b4"
)
```

- [ ] **Step 2: Update callset label mapping (lines 444-447)**

Replace:
```r
    callset_label = case_when(
      callset == "baylor_ont" ~ "Baylor ONT (Sniffles2)",
      callset == "dragen" ~ "Illumina DRAGEN"
    ),
```

With:
```r
    callset_label = case_when(
      callset == "ont-sniffles" ~ "ONT-Sniffles (Sniffles2)",
      callset == "ont-verkko" ~ "ONT-Verkko (dipcall)"
    ),
```

- [ ] **Step 3: Update size bin levels if needed (check extraction output first)**

If `svtype_size_counts.csv` uses bracket notation, update lines 30-31:
```r
# Check actual size bins from data
sv_size_levels <- unique(sv_size_df$size_bin)
# Or manually define:
# sv_size_levels <- c("[50,100)", "[100,300)", "[300,1000)", "[1,10000)", ">=10000")
```

- [ ] **Step 4: Update figure axis labels (lines 510-515, 582-585)**

Replace `"Baylor ONT (Sniffles2)"` with `"ONT-Sniffles"` and `"Illumina DRAGEN"` with `"ONT-Verkko"` in:
- Line 512: x-axis label in scatter plot
- Line 513: y-axis label in scatter plot
- Line 582-583: x/y axis labels in size scatter plot

- [ ] **Step 5: Test Quarto render**

```bash
quarto render analysis/use_case_evaluation.qmd
```

Expected: HTML output generated without errors

- [ ] **Step 6: Visual inspection of key plots**

Open `analysis/use_case_evaluation.html` and check:
- Callset labels show "ONT-Sniffles" vs "ONT-Verkko"
- Colors match expected mapping
- Size bins display correctly (with or without brackets)
- Data points appear reasonable (recall/precision in 0-1 range)

- [ ] **Step 7: Update results text (lines 617-627)**

Replace section starting "To illustrate how the v5.0q SV benchmark enables..." with:

```markdown
To demonstrate the practical impact of SV calling methodology on benchmarking performance, we compared mapping-based (ONT-Sniffles2) and assembly-based (ONT-Verkko dipcall) structural variant callsets against the v5.0q GRCh38 benchmark using `truvari bench` with GIAB-recommended parameters (`--pick ac --passonly -r 2000 -C 5000`) followed by `truvari refine` for MSA-based re-evaluation. ALT="*" spanning deletion records were filtered prior to benchmarking per GIAB v5.0q README guidance.

[Add performance comparison results after reviewing actual metrics]
```

- [ ] **Step 8: Commit Quarto notebook updates**

```bash
git add analysis/use_case_evaluation.qmd
git commit -m "feat: update SV use-case notebook for ONT callset comparison"
```

---

## Task 6: Verification and Documentation

**Files:**
- Update: `TODO.md`
- Create: `CHANGELOG.md` entry (if exists)

- [ ] **Step 1: Verify all CSV outputs have expected row counts**

```bash
wc -l results/use_case/stvar/*.csv
```

Expected:
- `stratified_metrics.csv`: ~11 rows (1 header + 5 contexts × 2 callsets)
- `svtype_metrics.csv`: ~10-20 rows (varies by SV types present)
- `svtype_size_counts.csv`: ~100+ rows (TP/FP/FN × callsets × types × bins)

- [ ] **Step 2: Verify figure generation**

```bash
ls -lh manuscript/figs/use_case_stvar.png
```

Expected: PNG file updated with new data

- [ ] **Step 3: Check metrics are reasonable**

```bash
# Quick sanity check on recall/precision ranges
tail -n +2 results/use_case/stvar/stratified_metrics.csv | cut -d, -f5-6 | sort -t, -k1 -n
```

Expected: Recall and precision values between 0.0 and 1.0

- [ ] **Step 4: Update TODO.md**

Remove or mark completed:
```markdown
- [x] Integrate SV use-case evaluation into Snakemake pipeline (deferred to future work)
```

Add future work item:
```markdown
### Pipeline Integration (Future)

- [ ] Integrate SV use-case evaluation into Snakemake pipeline (Approach 3)
  - Add `workflow/rules/use_case_evaluation.smk`
  - Add config entries for callsets, parameters, stratifications
  - Automate callset download and benchmarking
```

- [ ] **Step 5: Commit TODO updates**

```bash
git add TODO.md
git commit -m "docs: update TODO for SV use-case completion"
```

- [ ] **Step 6: Final verification - render all analysis notebooks**

```bash
quarto render analysis/
```

Expected: All notebooks render without errors

- [ ] **Step 7: Review git log for commit quality**

```bash
git log --oneline -7
```

Expected commit messages:
1. `docs: update TODO for SV use-case completion`
2. `feat: update SV use-case notebook for ONT callset comparison`
3. `feat: generate SV use-case metrics with ONT callsets`
4. `feat: use Truvari API for SV metrics extraction`
5. `feat: update SV use-case script for ONT callsets and Truvari stratify`
6. `chore: remove old SV use-case results for ONT callset update`

---

## Completion Criteria

**All tasks complete when:**

1. ✅ Old Baylor/DRAGEN results and callsets removed
2. ✅ `run_sv_use_case.sh` updated with:
   - ONT-Sniffles and ONT-Verkko callset paths
   - Google Drive copy step
   - Truvari stratify loop for HP/TR/SD/MAP
3. ✅ `extract_sv_metrics.py` updated with:
   - Truvari API imports (`truvari.performance_metrics`, `truvari.get_sizebin`, `truvari.vcf_to_df`)
   - Summary.json parsing for overall metrics
   - Stratify text file parsing for per-context metrics
   - Size binning via Truvari API
4. ✅ Benchmarking pipeline runs successfully for both callsets
5. ✅ Three CSV files generated with expected structure
6. ✅ `use_case_evaluation.qmd` updated with new callset labels and colors
7. ✅ Quarto notebook renders without errors
8. ✅ Figures show reasonable performance metrics
9. ✅ All changes committed with conventional commit messages
10. ✅ TODO.md updated with future Snakemake integration task

**Success metrics:**
- Pipeline completes in <2 hours total runtime
- CSV files have non-zero rows and valid metrics (0.0-1.0 range)
- Quarto HTML output displays updated callset labels
- Git history shows 6-7 focused commits

**Known limitations:**
- Size bin labels may differ from previous format (brackets vs ranges)
- Results text needs manual review and update after metrics inspection
- Small variant section unchanged (rtg vcfeval approach retained)
