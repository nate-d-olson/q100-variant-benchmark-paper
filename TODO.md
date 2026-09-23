# Manuscript TODO

Trimmed 2026-09-23 after submission to *Cell Genomics* and bioRxiv
(2026-09-22). Pre-submission history is in git (tag `v5-submission-raw`).

## Revision notes

- **Code state at submission:** tag `v5-submission` (clean) /
  `v5-submission-raw` (as-submitted working tree). The figure and table map
  is in `README.md` and `docs/figure-map.csv`.
- **Hand-edited items that need manual updates if the underlying data change:**
  - Fig 1 (`figures/manual/fig1_development_cycle.af`)
  - Fig 4 panels (`figures/manual/fig4*.af`, starting from the scripted SVGs)
  - Tables 2 and 3, which were typed in the manuscript. Table 3 mirrors the
    `exclusions` blocks in `config/config.yaml`.
- **Manuscript values:** re-run `analysis/manuscript_value_verification.qmd`
  after any pipeline change. The verification results against the 20260922
  text are below.
- **Fig 3 panel B:** SVbyEye inputs come from `scripts/prep_svbyeye_*.sh`,
  outside Snakemake. See `CLAUDE.md`, "Manuscript Figures Outside the
  Notebooks".

### Value verification against the 20260922 text

Re-checked 2026-09-23. All 24 in-scope values match the submitted text,
including the Fig 2 fold changes: SD 33.4× (1,868 vs 56 SVs) and TR 4.4×
(20,655 vs 4,693). Minor differences to fix at revision:

- [ ] Table 1 region Δ%: the text says +7.7 %, but it computes to +7.8 %
  (rounding).
- [ ] Table 4 labels: the text uses "TSPY2" and "SV regions (small variant
  benchmark only)", while the pipeline uses `TSPY2-segdups` and a different
  name for the SV regions. Align the wording.
- `guppy_version` and `snakemake_archive_doi` are no longer in the text; they
  are marked `out_of_scope`.

## Data release

Packaging requirements and validation checklists are in
`docs/dataset-release-plan.md`.

- [ ] Choose the hosting platform (Zenodo versus a GIAB-coordinated repository).
- [ ] Decide whether to publish the large variant Parquet files or only
  aggregated metrics plus regeneration instructions.
- [ ] Confirm that large, regenerable coverage BEDs should be excluded.
- [ ] Decide whether external-evaluation TSVs can be released as-is or require
  further curation/attribution.
- [ ] Confirm that pinned GIAB input URLs are stable, or mirror those inputs in
  the deposit.
- [ ] Decide whether to include CSV versions of smaller Parquet outputs.
- [ ] Build the staged release package (manifests, READMEs, citations,
  license). Verify `input_manifest/pipeline_inputs.tsv` coverage and the
  `MANIFEST.tsv` checksums from a clean download.
- [ ] Verify the Snakemake provenance archive URL/DOI on data.nist.gov.
- [ ] Add the bioRxiv and data release DOIs to `README.md` and `CITATION.cff`.

## Candidate analyses for revision

- [ ] Compact figure summarizing v5-only versus previous-only base and variant
  deltas.
- [ ] Annotate v5 benchmark VCFs with exclusion-based FILTER values.
- [ ] Add total included variation to the v5 statistics table and compare it
  against a baseline HG002 assembly, the reference, and the previous
  benchmarks.
- [ ] Numerical claims in `use_case_evaluation.qmd` (F1, recall, precision,
  false-negative counts) are hard-coded prose, not inline R. Verify them
  manually if that analysis is used in a response.

## Pipeline and repository maintenance

- [ ] Evaluate Snakevision for a pipeline diagram:
  <https://github.com/OpenOmics/snakevision>.
- [ ] Decide whether to extend the SV use-case Snakemake integration with
  download, native `truvari bench/refine/stratify`, and configuration rules.
- [ ] Optional: align bcftools versions across `workflow/envs/biotools.yaml`
  (1.22) and `truvari.yaml` (1.20).
