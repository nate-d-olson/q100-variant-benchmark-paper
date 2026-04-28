# Manuscript TODO

## Manuscript revisions

- [ ] text revisions
- [ ] confirm highlighted values
- [ ] update and add new references
- [ ] Workflow figure 1 revise based on feedback: Alternatively focus on v5 genchmark generation, move to supplemental? Bigger text
- [ ] exclusion figure move to methods or supplemental. Bigger text
- [ ] Ideogram
  - [ ] tracks: variant density, benchmark coverage, HG002 - ref alignment?, comparison to v0.6 and v4.2.1? - not sure how to present for GRCh37 and GRCh38
- Table v5 stats: total bp variation included, maybe a baseline comparison of the HG002 assembly compared to the ref genome and previous benchmarksets.
- [ ] figure v5 size, update legend, clearer differentiation between greys
- [ ] Chr 8 inversion, larger text, cut top to 100Mb, fix truncated legend.
- [ ] v5 coverage change figure: superseeded by ideogram???
- [ ] figure v5 var context, larger text
- [ ] Figure supplemental var size fold-change, log2 space
- [ ] figure pangene - add legends, better visualizations
- [ ] Figure eval curations: more consistent layout for B and D (either use facets or not)
- [ ] Figure RIDE-CI cleaer presentation of RIDE-CI
- [ ] Figure use_case_smvar - clearer presentation, revise based on feedback
- [ ] Figure use_case_stvar - clearer presentation, revise based on feedback
- [ ] Miqa platform evaluation figure



- [ ] Curate unsure SVs

### Figures

- [ ] **Better version of pangene graphs** (stretch goal)
- [ ] workflow figure legend

### Text

- [ ] figure legends (table legends drafted in `manuscript/tables.docx`)
- [ ] **Fill and verify highlighted values** — all placeholders (##, ??, XYZ, xx) in manuscript
- [ ] **Address JZ's Word comments** (14 comments)

## Pre-Submission

- [ ] **Verify Snakemake provenance archive URL/DOI** on data.nist.gov

### Analysis

- [ ] Add variants to exclusion tables in benchmarkset_characterization.qmd
- [ ] Small variant counts breakdown by <15bp, 15–49bp
- [ ] Look into small intervals (50–100bp) in v5q
- [ ] Add compact figure for v5-only vs previous-only base/variant deltas
- [ ] Annotate V5 benchmark VCFs with exclusion-based FILTER values

### Pipeline Integration

- [ ] Integrate SV use-case evaluation into Snakemake pipeline
  - Add `workflow/rules/use_case_evaluation.smk`
  - Add config entries for callsets, parameters, stratifications
  - Automate callset download and benchmarking

### Infrastructure

- [ ] Address documentation debt in `docs/` files
- [ ] Look into intervals <100bp (potential exclusions bug)