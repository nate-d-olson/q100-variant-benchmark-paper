# Manuscript TODO

Last consolidated: 2026-08-11. Sources reviewed: full codebase audit (no embedded
TODO/FIXME markers found in tracked files), TODO.md history, git log, and repo
cleanup completed 2026-08-10.

## Next steps (ordered)

1. [ ] Resolve the scientific-content blockers: curate unsure SVs and verify
   every highlighted manuscript value.
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
  Note: numerical claims in `use_case_evaluation.qmd` (F1, recall, precision,
  false-negative counts) are hardcoded prose strings, not live inline R
  expressions — verify manually against current pipeline outputs.
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
- [ ] Simplify the exclusion categories and finalize the main-text figure
  (`figures/exclusion_diagram.png`).
- [ ] Finalize the BED-operations figure number and legend for Supplemental
  Methods (`figures/exclusion_bed_operations.png`).
- [ ] Revise the Chr 8 inversion figure substantially. Improve text, annotations
  (MAT, PAT, and GRCh38), legend, and layout before deciding whether to retain the
  full-chromosome panel. None of the reviewed versions is publication quality.
- [ ] Update the v5 size figure legend and distinguish the gray series more
  clearly.
- [ ] Review the R 4.6 regenerated manuscript figures, especially
  `figures/var_len_dist.png`, and decide whether to retain the current PNGs and
  version their PDF counterparts.
- [ ] Enlarge text in the v5 variant-context figure.
- [ ] Show supplemental variant-size fold change in log2 space.
- [ ] Make the evaluation-curation panels B and D consistent (both faceted or
  both unfaceted).
- [x] Present the RIDE-CI result more clearly. Replaced geom_pointrange with
  dumbbell plot (grey segment + colored endpoints) and moved Overall facet to
  top in both `ci_forest_plot` and `cs_forest_plot` chunks. 2026-08-11.
- [ ] Revise the combined use-case figure based on small- and structural-variant
  feedback. `figures/use_case.png` is an alternate portrait layout with no
  generative code in the repo (manually created, confirmed orphan); decide whether
  to keep, regenerate, or retire it.
- [ ] Create the MIQA platform-evaluation figure.
- [ ] Improve Pangene graphs and legends (stretch goal). The implemented vector
  export is documented in `docs/pangene-vector-export.md`; the proposed
  reference-walk layout is in `docs/pangene-reference-walk-layout-plan.md`.

### Selected or superseded

- [x] Select the detailed vertical `ideogram_main` layout.
- [ ] Verify `figures/ideogram_main` labels at final manuscript width and settle
  the final tracks (variant density, benchmark coverage, assembly/reference
  alignment, and comparisons with v0.6/v4.2.1 as scientifically appropriate).
- [x] Select `combined_eval_strata` and `combined_eval_callset` as the canonical
  external-evaluation figures.
- [ ] Check manuscript references before retiring legacy `combined_eval` and
  `combined_callset_eval` exports.
- [x] Select `use_case_combined` as the canonical use-case figure.
- [x] Treat the standalone v5 coverage-change figure as superseded by the
  ideogram/combined coverage-and-variant-change figure.

## Analyses

- [x] Add variant counts to the exclusion tables in
  `analysis/benchmarkset_characterization.qmd`. Added `## Exclusion Variant Counts`
  section with smvar and stvar flextables (GRCh38 primary). 2026-08-10.
- [x] Add a small-variant count breakdown for <15 bp and 15–49 bp. Added
  `### Supplemental: Small Variant Count by Size Bin` to `benchmark_difficult.qmd`
  (GRCh38, v4.2.1 vs v5.0q, SNV/INDEL). Corrected the type filter after
  rendering on 2026-08-11.
- [x] Investigate 50–100 bp v5.0q benchmark intervals as a possible exclusion
  bug. Investigated 2026-08-10: **expected, not a bug.** The benchmark BED is
  downloaded directly from GIAB; this pipeline does not apply exclusions to
  produce it. v5.0q has only 15 intervals in the 50–99 bp range (0.05% of
  28,797 total), far fewer than v4.2.1 (18,474 of 481,622). They are legitimate
  narrow confident regions from the upstream DeFrABB pipeline. Documented in
  `docs/troubleshooting.md` under "50-100 bp Benchmark Intervals".
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

- [ ] Confirm the R 4.6.1 / Bioconductor 3.23 lockfile as the project baseline
  after reviewing rendered results; the former R 4.5.3 lockfile is retained as
  `renv.lock.r-4.5.3.bak`.
- [x] Remove `tests/unit/test_common_helpers.py`. Removed 2026-08-10.
- [x] Migrate the 2 remaining `gt()` tables in `analysis/external_evaluation.qmd`
  to flextable. Removed `theme_gt_manuscript`, `library(gt)` from `plot_themes.R`,
  and `load_gt` param from `_notebook_setup.R`. gt fully removed 2026-08-10.
- [ ] Evaluate Snakevision for a Snakemake pipeline diagram:
  <https://github.com/OpenOmics/snakevision>.
- [ ] Decide whether to extend the light SV use-case Snakemake integration with
  download, native `truvari bench/refine/stratify`, and configuration rules.

## GitHub issue reconciliation

- [x] Close [#49, "Create highlevel ideogram figure"](https://github.com/nate-d-olson/q100-variant-benchmark-paper/issues/49).
  Closed 2026-08-10.

## Recently verified or completed

- [x] Repo cleanup completed 2026-08-10: removed worktrees (reclaimed ~32 GB),
  deleted build artifacts, consolidated all figures into `figures/` + `figures/vector/`,
  moved tables to `tables/`, added `FIG_DIR` constant to `analysis/_notebook_setup.R`,
  updated all generators and docs, deleted stale branches, pruned remote refs.
- [x] Documentation updated 2026-08-10: CLAUDE.md, docs/architecture.md,
  workflow/README.md, docs/figure-manifest.csv, analysis/README.md all reflect
  the new `figures/` layout.
- [x] `analysis/external_evaluation.qmd` renders successfully (verified
  2026-07-16); the former missing-input/error note was stale.
- [x] Integrate SV use-case metric extraction into Snakemake (light scope).
- [x] Fix R tests and add the locked-renv CI job.
- [x] No embedded TODO/FIXME markers exist in any tracked source file (confirmed
  2026-08-10 audit); all open tasks are tracked here.
