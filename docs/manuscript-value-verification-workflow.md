# Manuscript Value Verification — Workflow and Lessons Learned

Record of the 2026-08-14 verification of all quantitative claims in the v5.0q
manuscript (`HG002v5-variant-benchmark-manuscript_20260814_revised.docx`),
kept as a reusable playbook for future manuscript checks. Task spec:
`manuscript-value-verification-prompt.md` (repo root). Deliverables:
`analysis/manuscript_value_verification.qmd` (renders end-to-end) and
`docs/agent_work/manuscript_value_updates.json` (25-entry substitution
manifest written by the notebook itself at render time).

## Process that worked

1. **Inventory first.** Write a task spec listing every value with its exact
   manuscript text, the comparison it implies, and the likely source notebook.
   Known-suspect inconsistencies get their own section so they are
   investigated, not just recomputed.
2. **Fan out by value group, replicate existing notebook logic.** Six
   parallel agents (one each for: summary table, exclusions table,
   coverage/contexts, fold-changes, external evaluation/GLM, Platinum
   Pedigree comparison). Each agent read the relevant `analysis/*.qmd` +
   `R/data_loading.R` first and replicated that exact logic (filters,
   groupings, rounding) rather than inventing methodology.
3. **Two-language cross-check.** Agents verified values in Python (pandas/
   pyarrow, directly against `results/**`); the final deliverable is an R
   Quarto notebook that recomputes everything through the repo's own loaders.
   Agreement between the two is a genuine independent check — it caught one
   real bug in proposed R code (a `system2()` shell-quoting failure that
   silently returned 1 instead of 4.3M).
4. **Notebook writes the JSON.** The `record_value()` accumulator pattern
   (one call per inventory item) builds both the human-review summary table
   (sorted mismatches-first) and the machine-readable JSON in the same
   render, so the two deliverables cannot drift.
5. **Batch the human decisions.** Items with more than one defensible answer
   were flagged `needs_human_review` with all candidates and evidence; the
   author resolved 9 of them in a single pass at the end. Do not silently
   pick a candidate.

## Why manuscript values mismatch — diagnostic checklist

Every mismatch found traced to one of these patterns; check them in order:

1. **Stale methodology.** If a value will not reproduce, `git log -p` the
   source notebook. Three mismatches traced to earlier code: the
   altered-bases formula changed (SNVs 0 bp → 1 bp, and a `/1e6` scaling was
   dropped); the abstract fold-changes (5.6×/4.6×/2×) predated the
   `BMKREGIONS`+`is_pass` filter tightening; the 69.1 % SV-in-repeats figure
   predated the Feb-2026 `normalize_annotation()` context-ID bugfix.
2. **Filtered vs unfiltered variant populations.** Pipeline per-context
   counts (`variants_by_genomic_context.parquet`) ignore `region_ids`/
   `is_pass`; notebook analyses filter to `BMKREGIONS` + PASS. Any % or
   fold-change must state which population it uses.
3. **Chromosome scope.** Autosome-only vs all-chromosome totals diverge badly
   when a legacy benchmark lacks chrX/Y entirely (v4.2.1: autosome-only
   increase 1.5 % vs whole-benchmark 7.8 %). Interval counts: manuscript used
   all-chromosome (28 797/9 658) while the notebook table is autosome-only.
4. **Region definition.** "Added sequence" candidates: final benchmark-region
   delta (+197 Mb) vs dip.bed-minus-legacy (+298 Mb). Name the regions.
5. **Count semantics.** Callsets *submitted* (11) vs *curated* (9) vs
   *analyzed* (7); data rows (70/40) vs header-inclusive line counts (71/41).
6. **Rounding provenance.** Δ% from full-precision inputs (+7.8 %) vs from
   already-rounded table values (197/2542 = +7.7 %).

## Environment notes for reruns

- `results/` and `resources/` must be populated (pipeline outputs are not in
  git); the verification notebook renders against them via the R loaders.
- renv activation hangs under the Claude Code sandbox — run `Rscript`/
  `quarto render` with the sandbox bypass and `RENV_CONSENT=yes`; never run R
  in the background.
- For sandbox-safe Python checks of parquet outputs, use the truvari env
  under `.snakemake/conda/` (has pandas + pyarrow).
- R gotcha: `system2("bash", c("-c", "<piped command>"))` does not quote the
  command string; use `system(sprintf(...), intern = TRUE)` for shell
  pipelines, and guard counts with a `stopifnot()` sanity check.

## JSON manifest schema

One entry per inventory row (including confirmations and out-of-scope items):
`id`, `section`, `current_text_snippet` (enough words to locate by search),
`current_value`, `calculated_value` (plain number, NIST formatting applied at
substitution time), `unit`, `status` (`placeholder_resolved` |
`confirmed_match` | `inconsistency_found` | `needs_human_review` |
`out_of_scope`), `confidence`, `needs_human_review`, `source`, `notes`
(candidates, decisions, and provenance of stale values).
