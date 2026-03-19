# Manuscript TODO

## New list

- [ ] Strech goal: Linear model for confidence intervals on curations (external_evaluation.qmd)
- [ ] Update difficult variant counts to include complex variants



## Pre-Draft Sharing (co-author review)

- Curate unsure SVs

### Figures

- [ ] Only include UCLA ensamble results in manuscript
- [ ] **Better version of pangene graphs** (stretch goal)
- [ ] workflow figure legend
- [ ] write tables to docx for inclusion
- [x] **Update stats table** (Table 1) — consolidated flextable
- [x] **Update variant size distribution figure** (Figure 2) — overlaid histograms, GRCh37 for SV
- [x] **Simplify variant context figure** (Figure 3) — fold change bars with exploded context_ids
- [x] **Combine v5 coverage change with variant count comparison** (Figure 4) — merge into single figure
- [x] **Better presentation of external evaluation strata** (Figure 5 / S3) — readable labels, UCLA ensemble only
- [x] **Clean-up chr8 inversion figure** (Figure S1) — PAV excluded region added, panels stacked
- [x] **Create/locate Figure 1** — DeFrABB workflow diagram

### Text
- [ ] figure legends
- [ ] revise use-case results
- [ ] **Fill and verify highlighted values** — all placeholders (##, ??, XYZ, xx) in manuscript
- [x] **Add chr8 inversion specifics** — coordinates from PAV callset (chr8:8,237,843–12,234,345)
- [x] **Benchmarking use-case section** — draft text in `analysis/use_case_evaluation.qmd`
- [ ] **Replace Discussion section** with content from `manuscript/discussion.qmd`
- [ ] **Remove Sina's embedded email** (paras 352–366 in docx)
- [ ] **Remove informal notes and outline fragments** (18 items per cleanup checklist)
- [ ] **Apply figure/table numbering scheme** (cross-reference map in numbering doc)
- [ ] **Address JZ's Word comments** (14 comments)

### Metadata / Admin

- [x] **Finalize author list**
- [ ] **Fill iPSC cell line identifiers** (check GIAB/Coriell docs)
- [ ] OTHER: update script/make_chr8_figure.py - comments consistent with script. 

## Post-Draft (lower priority)

- [ ] **Verify Snakemake provenance archive URL/DOI** on data.nist.gov

### Analysis

- [ ] Add variants to exclusion tables in benchmarkset_characterization.qmd
- [ ] Small variant counts breakdown by <15bp, 15–49bp
- [ ] Look into small intervals (50–100bp) in v5q
- [ ] Add compact figure for v5-only vs previous-only base/variant deltas
- [ ] Annotate V5 benchmark VCFs with exclusion-based FILTER values

### Pipeline Integration

- [ ] Integrate SV use-case evaluation into Snakemake pipeline (new rules + config entries for callsets, truvari params, metrics extraction)

### Infrastructure

- [ ] Address documentation debt in `docs/` files
- [ ] Look into intervals <100bp (potential exclusions bug)

