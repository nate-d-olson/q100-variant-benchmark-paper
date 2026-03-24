# Manuscript TODO

- [ ] Curate unsure SVs

### Figures

- [ ] **Better version of pangene graphs** (stretch goal)
- [ ] workflow figure legend

### Text

- [ ] figure legends (table legends drafted in `manuscript/tables.docx`)
- [ ] **Fill and verify highlighted values** — all placeholders (##, ??, XYZ, xx) in manuscript
- [ ] **Address JZ's Word comments** (14 comments)

## Post-Draft (lower priority)

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