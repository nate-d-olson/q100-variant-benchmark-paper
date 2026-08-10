# Manuscript TODO

Last consolidated: 2026-07-16. Sources reviewed: repository notes and code
comments, open GitHub issues, figure review decisions, and
`docs/agent_work/dataset-release-plan.md`.

## Next steps (ordered)

1. [ ] Resolve the scientific-content blockers: curate unsure SVs, investigate
   50–100 bp benchmark intervals, and verify every highlighted manuscript value.
2. [ ] Revise the main exclusion-categories figure and the Chr 8 figure; neither
   is ready for final manuscript placement.
3. [ ] Finalize the selected external-evaluation, ideogram, and combined use-case
   figures at publication dimensions, then update their legends and manuscript
   references.
4. [ ] Address the 14 Word comments and complete the general text/reference pass.
5. [ ] Make the data-release decisions listed below, then build and validate the
   release package.

## Manuscript text and scientific review

- [ ] Complete the manuscript text revision pass.
- [ ] Fill and verify all highlighted values and placeholders in the manuscript.
- [ ] Address JZ's 14 Word comments.
- [ ] Update and add references.
- [ ] Complete all figure legends; table legends are drafted in
  `tables/tables.docx`.
- [ ] Curate unsure structural variants.
- [ ] Verify the Snakemake provenance archive URL/DOI on data.nist.gov.

## Main and supplemental figures

### Revisions required

- [ ] Revise workflow Figure 1 based on feedback: focus it on v5 benchmark
  generation, enlarge the text, and decide whether it belongs in the supplement.
- [ ] Simplify the exclusion categories and finalize the main-text figure.
- [ ] Finalize the BED-operations figure number and legend for Supplemental
  Methods.
- [ ] Revise the Chr 8 inversion figure substantially. Improve text, annotations
  (MAT, PAT, and GRCh38), legend, and layout before deciding whether to retain the
  full-chromosome panel. None of the reviewed versions is publication quality.
- [ ] Update the v5 size figure legend and distinguish the gray series more
  clearly.
- [ ] Enlarge text in the v5 variant-context figure.
- [ ] Show supplemental variant-size fold change in log2 space.
- [ ] Make the evaluation-curation panels B and D consistent (both faceted or
  both unfaceted).
- [ ] Present the RIDE-CI result more clearly.
- [ ] Revise the combined use-case figure based on small- and structural-variant
  feedback; retain `use_case.png` as the alternate portrait layout and identify
  its source.
- [ ] Create the MIQA platform-evaluation figure.
- [ ] Improve Pangene graphs and legends (stretch goal). The implemented vector
  export is documented in `docs/pangene-vector-export.md`; the proposed
  reference-walk layout is in `docs/pangene-reference-walk-layout-plan.md`.

### Selected or superseded

- [x] Select the detailed vertical `ideogram_main` layout.
- [ ] Verify `ideogram_main` labels at final manuscript width and settle the
  final tracks (variant density, benchmark coverage, assembly/reference
  alignment, and comparisons with v0.6/v4.2.1 as scientifically appropriate).
- [x] Select `combined_eval_strata` and `combined_eval_callset` as the canonical
  external-evaluation figures.
- [ ] Check manuscript references before retiring legacy `combined_eval` and
  `combined_callset_eval` exports.
- [x] Select `use_case_combined` as the canonical use-case figure.
- [x] Treat the standalone v5 coverage-change figure as superseded by the
  ideogram/combined coverage-and-variant-change figure.

## Analyses

- [ ] Add variant counts to the exclusion tables in
  `analysis/benchmarkset_characterization.qmd`.
- [ ] Add a small-variant count breakdown for <15 bp and 15–49 bp.
- [ ] Investigate 50–100 bp v5.0q benchmark intervals as a possible exclusion
  bug; document whether the intervals are expected and fix the pipeline if not.
- [ ] Add a compact figure summarizing v5-only versus previous-only base and
  variant deltas.
- [ ] Annotate v5 benchmark VCFs with exclusion-based FILTER values.
- [ ] Add total included variation to the v5 statistics table and evaluate a
  baseline HG002-assembly comparison with the reference and previous benchmarks.

## Data release

- [ ] Detailed packaging requirements and validation checklists live in
`docs/agent_work/dataset-release-plan.md`.

### Author decisions needed

- [ ] Choose the hosting platform (Zenodo versus a GIAB-coordinated repository).
- [ ] Decide whether to publish the large detailed variant Parquet files or only
  aggregated metrics plus regeneration instructions.
- [ ] Confirm that large, regenerable coverage BEDs should be excluded.
- [ ] Decide whether external-evaluation TSVs can be released as-is or require
  further curation/attribution.
- [ ] Confirm that pinned GIAB input URLs are stable or mirror those inputs in
  the deposit.
- [ ] Decide whether to include CSV versions of smaller Parquet outputs.

### Packaging and validation

- [ ] Complete the staged release package, manifests, README files, citations,
  and CC0 licensing described in the release plan.
- [ ] Verify required benchmark, exclusion-impact, exclusion-interaction,
  cross-version, reference-size, and external-evaluation files.
- [ ] Verify `input_manifest/pipeline_inputs.tsv` coverage and all `MANIFEST.tsv`
  checksums using a clean download.

## Pipeline and repository maintenance

- [ ] Repair or remove `tests/unit/test_common_helpers.py`. It imports the
  non-importable Snakemake `common` module and tests two helpers that no longer
  exist; two tests fail, and CI currently ignores the file.
- [ ] Migrate the remaining `gt()` tables in
  `analysis/external_evaluation.qmd` to flextable. These are the last blocker to
  removing `gt` from the project.
- [ ] Address documentation debt in `docs/`; start by reconciling historical
  plans/specifications with current implementation status.
- [ ] Evaluate Snakevision for a Snakemake pipeline diagram:
  <https://github.com/OpenOmics/snakevision>.
- [ ] Decide whether to extend the light SV use-case Snakemake integration with
  download, native `truvari bench/refine/stratify`, and configuration rules.

## GitHub issue reconciliation

- [ ] Close [#49, “Create highlevel ideogram figure”](https://github.com/nate-d-olson/q100-variant-benchmark-paper/issues/49),
  after noting that `scripts/make_ideogram.R` now produces the selected
  `figures/ideogram_main` figure. This is the repository's only open
  issue as of 2026-07-16.

## Recently verified or completed

- [x] `analysis/external_evaluation.qmd` renders successfully (verified
  2026-07-16); the former missing-input/error note was stale.
- [x] Integrate SV use-case metric extraction into Snakemake (light scope).
- [x] Fix R tests and add the locked-renv CI job.
