# Manuscript TODO

## Pre-Draft Sharing (co-author review)

### Figures

- [ ] **Update stats table** (Table 1) — consolidate summary stats
- [ ] **Update variant size distribution figure** (Figure 2) — verify v0.6 bars look correct relative to v5
- [ ] **Simplify variant context figure** (Figure 3) — emphasize points made in text
- [ ] **Combine v5 coverage change with variant count comparison** (Figure 4) — merge into single figure
- [ ] **Better version of pangene graphs**
- [ ] **Better presentation of external evaluation strata** (Figure 5 / S3)
- [ ] **Clean-up chr8 inversion figure** (Figure S1) — use PAV callset inversion coordinates for specifics
- [ ] **Create/locate Figure 1** — DeFrABB workflow diagram

### Text

- [ ] **Fill and verify highlighted values** — all placeholders (##, ??, XYZ, xx) in manuscript
- [ ] **Add chr8 inversion specifics** — coordinates from PAV callset
- [ ] **Benchmarking use-case section** — similar to pFDA manuscript
- [ ] **Replace Discussion section** with content from `manuscript/discussion.qmd`
- [ ] **Remove Sina's embedded email** (paras 352–366 in docx)
- [ ] **Remove informal notes and outline fragments** (18 items per cleanup checklist)
- [ ] **Apply figure/table numbering scheme** (cross-reference map in numbering doc)
- [ ] **Address JZ's Word comments** (14 comments)

### Metadata / Admin

- [ ] **Finalize author list**
- [ ] **Fill iPSC cell line identifiers** (check GIAB/Coriell docs)
- [ ] **Decide: add UCLA stvar callsets to Table 1** or keep at 6?
- [ ] **Verify Snakemake provenance archive URL/DOI** on data.nist.gov

## Completed

### 2026-03-18

- [x] Merged PR #47: correct variant type classification for UNK variants
- [x] Fixed CI (ruff lint/format, air R formatting)
- [x] Regenerated analysis notebooks and figures after variant type fix
- [x] Added `R/bed_helpers.R` (extracted from benchmark_unique_regions.qmd)

### 2026-02-17

- [x] Verified all placeholder numbers from pipeline results
- [x] Drafted co-author email → see `docs/coauthor_email_draft.md`
- [x] Wrote full Discussion section → `manuscript/discussion.qmd`
- [x] Fixed FTP version mismatch (V0.019 → V0.020)
- [x] Reconciled stvar callset count discrepancy
- [x] Updated results.qmd display items to match numbering scheme
- [x] Removed pipeline-internal section from data.qmd
- [x] Deleted stale diff_variant_counts.png

## Post-Draft (lower priority)

### Analysis

- [ ] Add variants to exclusion tables in benchmarkset_characterization.qmd
- [ ] Small variant counts breakdown by <15bp, 15–49bp
- [ ] Look into small intervals (50–100bp) in v5q
- [ ] Validate summary outputs against old_only status CSVs
- [ ] Add compact figure for v5-only vs previous-only base/variant deltas
- [ ] Linear model for confidence intervals on curations (external_evaluation.qmd)
- [ ] Annotate V5 benchmark VCFs with exclusion-based FILTER values
- [ ] Update difficult variant counts to include complex variants

### Infrastructure

- [ ] Transition from regex-based `parse_benchmark_id` to structured metadata sidecars
- [ ] Evaluate replacing custom caching with R `targets` package
- [ ] Address documentation debt in `docs/` files
- [ ] Cleanup `data_loading.R` script
- [ ] More sensible results and log directory structures
- [ ] Look into intervals <100bp (potential exclusions bug)
- [ ] Remove genome context worktree and local branch
