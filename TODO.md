# Codebase Review Findings (Feb 2026)
Remove genome context worktree and local branch

# Manuscript Preparation

## Completed (2026-02-17)

- [x] Verified all placeholder numbers from pipeline results
      → see `docs/manuscript_number_verification.md`
- [x] Drafted co-author email → see `docs/coauthor_email_draft.md`
- [x] Wrote full Discussion section → `manuscript/discussion.qmd`
- [x] Created docx cleanup checklist → `docs/docx_cleanup_checklist.md`
- [x] Created figure/table numbering scheme → `docs/figure_table_numbering.md`
- [x] Fixed FTP version mismatch (V0.019 → V0.020 in data.qmd, methods.qmd)
- [x] Reconciled stvar callset count discrepancy (updated number_verification.md)
- [x] Updated results.qmd display items to match numbering scheme
- [x] Removed pipeline-internal section from data.qmd
- [x] Deleted stale diff_variant_counts.png

## Remaining — Docx Edits (manual)

- [ ] Fill all placeholder values (##, ??, XYZ, xx) per cleanup checklist
- [ ] Remove Sina's embedded email (paras 352–366)
- [ ] Remove informal notes and outline fragments (18 items in checklist)
- [ ] Apply figure/table numbering scheme (cross-reference map in numbering doc)
- [ ] Address Word comments from JZ (14 comments)
- [ ] Replace Discussion section with content from discussion.qmd
- [ ] Create/locate Figure 1 (DeFrABB workflow diagram)
- [ ] Create Figure S1 (chr8 inversion visualization)
- [ ] Finalize author list
- [ ] Fill iPSC cell line identifiers (check GIAB/Coriell docs)
- [ ] Decide: add UCLA stvar callsets to Table 1 or keep at 6?
- [ ] Verify "Snakemake provenance archive: data.nist.gov" has actual URL/DOI

# Analysis

- v0.6 size distribution compared to v5 the bars don't look quite right, v0.6
  has more variants than expected compared to v5 based on previous plots and
  expectations. Need to verify.

### R Analysis & Data Loading

- [ ] **Robust Metadata Handling**: Transition from regex-based filename parsing (`parse_benchmark_id`) to structured metadata sidecars (JSON/YAML) for pipeline outputs.
- [ ] **Simplify Data Loading**: Evaluate replacing the custom caching/validation logic in `R/data_loading.R` with a standard framework like the R `targets` package.
- [ ] **Documentation Debt**: Address internal TODOs and "documentation debt" found within `docs/` files and notebook comments.


## Analysis Revisions

### benchmarkset_characterization.qmd

- Add variants to the exclusion tables
- Small variant counts breakdown by <15bp , 15 - 49bp
- Look into number and distribution of small intervals (50 - 100bp) in v5q

### benchmark_unique_regions.qmd

- Validate summary outputs against old_only status CSVs
- Add compact figure for v5-only vs previous-only base/variant deltas


### manuscript/results.qmd

- **Figure 2**: Small and structural variant counts, coverage, and comparison to previous benchmark sets. %%TODO%% Stratified, annotation based???

### analysis/external_evaluation.qmd

- **TODO**: linear model to get confidence intervals for curations

### Low priority

- more sensible results and log directory structures
- cleanup data_loading script
- look into intervals, less than 100bp (potential bug in exclusions) in small variant benchmark regions bed files for v5 benchmarkset.

## Other

- Annotating V5 benchmark vcfs based on exclusions and add filter field values for variants outside the benchmark regions.
- Updating difficult variant counts to include complex variants
