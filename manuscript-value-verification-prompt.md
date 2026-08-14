# Verify and Fill In Manuscript Result Values (HG002 v5.0q Benchmark Paper)

## Purpose

The HG002 v5.0q variant benchmark manuscript (`HG002v5-variant-benchmark-manuscript_20260814_revised.docx`, lives outside this repo in the parent `v5-manuscript/` directory) has gone through NIST formatting review and contains a mix of:

1. **Unresolved placeholder values** (`##`) that were never filled in.
2. **Concrete values** that were hand-transcribed from earlier pipeline runs and may now be stale relative to the current `results/` outputs and `analysis/*.qmd` notebooks in this repo.

Your job is to use this repository's existing pipeline outputs and analysis notebooks to calculate/confirm every quantitative claim listed below, and produce two deliverables (see "Deliverables" section). **Do not edit the manuscript `.docx` file** — that substitution step happens afterward, by hand, using your JSON output.

## Ground rules (read first)

1. **Prefer existing notebook logic over new code.** `analysis/*.qmd` already computes most of these numbers for the figures/tables currently in the manuscript. Re-deriving a value with new ad hoc code risks silently using different methodology (different exclusion set, different rounding, different reference) than what's already in the manuscript's figures — which would create a *new* inconsistency instead of resolving one. For each value, first check whether an existing `.qmd` chunk, `analysis/cache/`, or `results/**` file already produced it, and reuse/extend that rather than rewriting from scratch.
2. **Flag ambiguity instead of guessing.** A few items below (callset counts, in particular) have genuine ambiguity in what should be counted. If more than one defensible number exists, report all candidates with your reasoning in the notebook and mark the JSON entry `"needs_human_review": true` rather than silently picking one. **Stop and ask the user before finalizing any item marked `needs_human_review` — do not silently resolve it and move on.**
3. **Don't recompute externally-provided provenance.** Sequencing depth/coverage figures, software versions, and basecaller settings quoted in the Methods section for externally-submitted datasets (e.g., "downsampled to ≈ 30×", "80× ONT coverage") are self-reported by the contributing groups, not something this pipeline calculates. Don't try to verify these by downloading/re-analyzing raw BAMs — out of scope. These are listed under "Out of scope" below.
4. **Match existing NIST number formatting in your output values.** The manuscript has already been passed through NIST SI formatting (thin non-breaking space digit grouping instead of commas, `%`/unit spacing, `≈` for approximations, `×` for fold-change/coverage, en-dash ranges converted to "N unit to M unit"). When you report a calculated value, give the plain numeric value (e.g. `699449`) in JSON — do not hand-format it — but in the notebook's human-readable table, present it NIST-style (e.g. `699 449`) so it's easy for a human to visually diff against manuscript text.
5. **Two outputs only.** A notebook (for human review) and a JSON file (for programmatic substitution). Do not touch the `.docx`.

## Repository orientation

- Environment: `mamba activate q100-smk` (Snakemake 8.x env). Quarto notebooks in `analysis/` already render against existing `results/**` outputs — you should **not** need to re-run the Snakemake pipeline; `results/` is already populated (last generated 2026-08-10).
- Key notebooks and what they already compute:
  | Notebook | Relevant existing sections |
  |---|---|
  | `analysis/benchmarkset_characterization.qmd` | Main-text summary table (source of **Table v5.0 stats**), variant size distributions, benchmark region coverage change by chromosome |
  | `analysis/benchmark_difficult.qmd` | Variant fold-change by genomic context (source of abstract's `5.6×`/`4.6×`/`2×` claims), fold-change by variant size bin |
  | `analysis/benchmark_exclusions.qmd` | Bases removed by exclusion (source of **Table XYZ**), exclusion interactions/upset plots |
  | `analysis/benchmark_interval_size_distributions.qmd` | Interval counts and sizes (source of the `9 658` / `28 797` region-count claims vs. Platinum Pedigree) |
  | `analysis/benchmark_unique_regions.qmd` | v5-only vs. previous-only region/variant breakdowns |
  | `analysis/genomic_context_analysis.qmd` | % of genome/context covered by benchmark, coverage change across versions, variant counts by context |
  | `analysis/external_evaluation.qmd` | External curation counts, and the **binomial GLM** (bias-reduced logistic regression via `brglm2`, marginal estimates via `emmeans`) — this is the exact model behind the "one-sided 95 % lower confidence bound" claims |
  | `analysis/use_case_evaluation.qmd` | NeuSomatic/Sniffles use-case metrics — **note:** the corresponding manuscript section and figures (`use_case_smvar`, `use_case_stvar`) were removed from the 20260814 draft, so this notebook's output is likely no longer needed for the current manuscript text. Confirm before spending time here. |
  | `R/data_loading.R` | Caching helpers (`cache_info()`, `invalidate_cache()`) if cached values seem stale |
- Reference sizes: `results/ref_genome_sizes/{GRCh37,GRCh38,CHM13v2.0}_size.tsv`
- External evaluation raw curation files: `data/external-evaluations/*.csv` and `data/external-evaluations/Q100-ext-evals-2025-04-11-{smvars,stvars}.tsv`
- Exclusion metrics (already computed, source of truth for Table XYZ): `results/exclusions/*/exclusion_impact.csv`

## Value inventory

Each item below includes the manuscript section, the exact current text, and where in this repo to find/confirm the number. Treat the "current value" as what's *in the text right now* — your job is to state what you calculate and whether it matches.

### A. Unresolved placeholders (`##`) — must be resolved

| ID | Section | Current text (exact) | What to calculate | Likely source |
|---|---|---|---|---|
| `autosome_pct_increase` | Results — Benchmark content | "The v5.0q benchmark set includes ## % more of the autosomes than v4.2.1" | % increase in autosome bp covered, v5.0q GRCh38 vs. v4.2.1 GRCh38 | `genomic_context_analysis.qmd` or `benchmarkset_characterization.qmd` (region coverage by chromosome) |
| `sv_difficult_context_pct_increase` | Results — Benchmark content | "The v5.0 SV benchmark set includes ## % more variants in difficult genomic context" | % increase in SV variant count within difficult genomic context, v5.0q GRCh37 vs. v0.6 GRCh37 | `benchmark_difficult.qmd` (Variant Fold Change by Genomic Context, stvar) |
| `sv_reference_coverage_pct_increase` | Results — Benchmark content | "and covers ## % more of the reference genome" | % increase in reference genome bp covered by SV benchmark regions, v5.0q GRCh37 vs. v0.6 GRCh37 | `genomic_context_analysis.qmd` / `benchmarkset_characterization.qmd` |
| `small_var_glm_lcb_pct` | External Evaluation and Validation | "greater than ## % for the small variant benchmark set" | Minimum one-sided 95 % lower confidence bound across all strata, small-variant GLM | `external_evaluation.qmd` "Binomial Model" section |
| `sv_glm_lcb_pct` | External Evaluation and Validation | "and ## % for the structural variant benchmark set" | Same, structural-variant GLM | `external_evaluation.qmd` "Binomial Model" section |
| `ext_eval_smvar_callset_count` | External Evaluation Manual Curation (Methods) | "compared to ## and ## small and structural variant callsets" (first `##`) | Count of small variant callsets in external evaluation | See Section B below — cross-check against Table EvalCallsets and other mentions first |
| `ext_eval_stvar_callset_count` | External Evaluation Manual Curation (Methods) | same sentence, second `##` | Count of structural variant callsets in external evaluation | See Section B below |
| `guppy_version` | Methods — heading | "Input data: 60× Guppy (??version) SUPv5 model" | Not calculable from this repo — basecaller version used by the data-submitting group. Look for it in `data/` provenance/README for that ONT dataset, or flag for the manuscript author to supply. | N/A — likely **out of scope**, see below |
| `snakemake_archive_doi` | Data and Code availability | "Snakemake provenance archive: TODO data.nist.gov" | Not calculable — this is a pending data.nist.gov deposit handle. Flag as a manual action, do not attempt to compute. | N/A — **out of scope** |

### B. Known candidate inconsistencies (already spotted — please resolve, don't just recompute blindly)

These came up while cross-referencing the manuscript text against itself and against this repo. Each is a real discrepancy between two things that should agree.

1. **Small-variant benchmark size increase: "≈ 600 000" vs. table's "+699 449".**
   Results text (para "Benchmark content and genome coverage"): *"The small variant benchmark set includes ≈ 600 000 more variants in GRCh38 compared to v4.2.1."*
   **Table v5.0 stats**, Total row, small variants: GRCh38 Δn = **+699 449** (Δ% = +17.8 %), per the table's own footnote comparing v5.0q(GRCh38) vs. v4.2.1(GRCh38).
   These don't match (≈600K vs. ~700K). Recompute the small-variant total delta directly from `results/variant_tables/v5.0q_GRCh38_smvar/variants.parquet` vs. the v4.2.1 equivalent and report which figure (≈600K, ≈700K, or something else) is correct, and where the discrepancy likely originated.

2. **Abstract's "2-fold" SV increase vs. Table's "+192.5 %".**
   Abstract: *"...and 2-fold for structural variants..."* Results text: *"...more than twice the number of SVs..."* **Table v5.0 stats**, Total SVs row: Δ% = **+192.5 %** (v5.0q GRCh37 vs. v0.6 GRCh37), which implies v5.0q has **~2.9×** as many SVs as v0.6 — closer to "3-fold" than "2-fold". Confirm the +192.5 % figure against `results/variant_tables/`, and flag whether "2-fold" in the abstract should be "~3-fold" (or whatever the confirmed multiple is) for consistency with the table.

3. **External SV callset count: "7" vs. "11".**
   Results/External Evaluation: *"...8 small variant callsets and 7 structural variant callsets..."*
   Methods/External Callset Generation Methods: *"...compared to 8 small variant and 11 structural variant callsets..."*
   **Table EvalCallsets** (Table 1 in the docx) currently lists exactly **7 unique structural-variant Callset IDs** (`comenius-svdss2`, `illumina-dragen`, `iter-dragen`, `baylor-ont`, `pacbio-sawfish`, `baylor-hifi`, `UCLA Ensemble`) and **8 unique small-variant Callset IDs**. The Methods GLM section separately states the SV GLM used **N = 280 observations from 7 callsets**, which agrees with "7", not "11".
   However, `data/external-evaluations/` contains curation files for `Ucla-Delly-Newvcf` and `Ucla-Manta-Newvcf` in addition to `Ucla-Ensamble`, which aren't reflected in Table EvalCallsets — these may be earlier/superseded submissions that were later consolidated into "UCLA Ensemble", which could explain how an older draft arrived at "11". Please: (a) determine whether Delly/Manta were merged into "UCLA Ensemble" or are genuinely separate additional callsets that should be counted, and (b) report both the small- and structural-variant callset counts with your reasoning. This one likely needs a human decision (flag `needs_human_review: true`) rather than a pure recount, since it depends on which historical submissions are still "live."

4. **Curation sample size: "40"/"70" (manuscript) vs. "41"/"71" (code comment).**
   Manuscript: *"manually curated 40 variants for each submitted SV callset and 70 variants for small variant callsets."* Manuscript also separately states *"Five variants from the 8 strata were randomly sampled from each external structural variant callset"* — 8 × 5 = 40, internally consistent with "40".
   However, `analysis/external_evaluation.qmd` (Data Loading section) has a code comment: *"small variants expect 71, structural variants 41... a few files have less as less than 5 variants per strata to subset"* — i.e., the code's expected row counts are 71/41, one more than the manuscript's 70/40. This is very likely a header-row-vs-data-row off-by-one (71 including a header = 70 data rows) rather than a real discrepancy, but please confirm by actually counting data rows (excluding header) in `data/external-evaluations/*.csv` / the `Q100-ext-evals-2025-04-11-*.tsv` files and report the confirmed per-callset curated counts (ideally the actual min/max/mean across callsets, not just the nominal target, since the same qmd comment notes some files have fewer).

5. **Platinum Pedigree comparison: "≈ 4.5 million" vs. v5.0q's actual small variant total.**
   Text: *"each including ≈ 4.5 million small variants spanning ≈ 2.7 Gbp of GRCh38."* Table v5.0 stats shows v5.0q GRCh38 small variant total = 4 625 918 (≈ 4.6 million). Confirm whether "≈ 4.5 million" was meant to describe the Platinum Pedigree benchmark specifically (which may indeed be ≈4.5M, a different, external number) rather than v5.0q, and that the sentence isn't conflating the two. If PP's actual count is available anywhere in `resources/platinum-pedigree-data/`, extract and report it alongside v5.0q's actual count for direct comparison.

### C. Headline results to re-verify (tables + abstract) — high scrutiny, values already exist but should be reconfirmed against current pipeline outputs

| Location | Current value(s) | Source to confirm against |
|---|---|---|
| **Table v5.0 stats** — all cells (SNV, INDEL<50bp, Total, Altered bases, Regions incl., %Reference for small variants; DEL≥50bp, INS≥50bp, Total SVs, Total Variation, Region, %Reference for SVs; across GRCh37/GRCh38/CHM13v2.0, plus Δn/Δ% vs. legacy) | See full current values transcribed below | `benchmarkset_characterization.qmd` "Main Text: Summary (Table 1)" — this is the direct source; recompute and diff cell-by-cell |
| **Table XYZ** (Bases removed by exclusion, GRCh38) — all 16 exclusion categories, bp removed + % of dip.bed | See full current values transcribed below | `results/exclusions/v5.0q_GRCh38_{smvar,stvar}/exclusion_impact.csv` — direct source of truth |
| Abstract: `≈ 300 Mbp` added sequence | vs. v4.2.1 | Region size difference, `benchmarkset_characterization.qmd` |
| Abstract: `18 %` small variant increase | Table shows Total Δ% = +17.8 % (rounds to 18 %, consistent) | Confirm rounding is intentional |
| Abstract: `5.6×` SD, `4.6×` TR, `2×` homopolymer fold-change | | `benchmark_difficult.qmd` fold-change by genomic context, v5.0q vs. v4.2.1 |
| Results: `1.2 %` of genome in pangenome-derived LC regions (excl. satellites) | | `genomic_context_analysis.qmd` coverage summary |
| Results: `69.1 %` of v5.0q SVs in those repetitive regions | | `genomic_context_analysis.qmd` variant counts by context (stvar, % of total) |
| Results: chr8 (and now also chr1, chr16) show lower v5.0q coverage vs. v0.6 | text was recently edited to add chr1 and chr16 to the exception list — confirm all three (and only these three) chromosomes actually show a coverage decrease | `benchmarkset_characterization.qmd` "Benchmark Region Coverage Change by Chromosome" |
| Comparison to PP: `9 658` SV benchmark regions, `28 797` small variant benchmark regions | | `benchmark_interval_size_distributions.qmd` "Interval Counts"/"Interval Summary Table" |
| Methods: GLM `N = 350` (5 callsets, small variant), `N = 280` (7 callsets, SV) | | `external_evaluation.qmd` — confirm observation counts and callset counts match the fitted model's input data |
| Methods: "5 failures across 350 observations" (Firth sensitivity analysis) | | `external_evaluation.qmd` sensitivity analysis |

**Current Table v5.0 stats values (for reference — recompute independently, don't just copy these):**

| Metric | GRCh37 | GRCh38 | CHM13v2.0 | Δn (vs. legacy) | Δ% |
|---|---|---|---|---|---|
| SNV | 3651855 | 3671436 | 3395941 | +303973 | +9.0% |
| INDEL <50bp | 951308 | 954482 | 885454 | +395476 | +70.7% |
| Total (small var) | 4603163 | 4625918 | 4281395 | +699449 | +17.8% |
| Altered bases (Mb) | 3.6 | 3.6 | 3.4 | +2.0 Mb | +118.9% |
| Regions incl. (Mb) | 2726 | 2739 | 2768 | +197 Mb | +7.7% |
| % Reference (small var) | 95.4% | 93.7% | 88.8% | — | +6.7 pp |
| DEL ≥50bp | 10561 | 10852 | 12841 | +6491 | +159.5% |
| INS ≥50bp | 17571 | 17763 | 12204 | +12024 | +216.8% |
| Total SVs | 28132 | 28615 | 25045 | +18515 | +192.5% |
| Total Variation (Mb) | 19.1 | 18.6 | 16.1 | +11.8 Mb | +161.3% |
| Region (Mb, SV) | 2742 | 2756 | 2820 | +229 Mb | +9.1% |
| % Reference (SV) | 95.9% | 94.3% | 90.5% | — | +8.0 pp |

*(Footnote: Δ compares v5.0q(GRCh38) vs. v4.2.1(GRCh38) for small variants; v5.0q(GRCh37) vs. v0.6(GRCh37) for SVs.)*

**Current Table XYZ values (Bases removed by exclusion, GRCh38 — for reference):**

| Exclusion | Bases excluded | % of dip.bed |
|---|---|---|
| Segmental Duplications | 52311253 | 1.841 |
| SV regions | 29069865 | 1.023 |
| PAV-inversions | 23518680 | 0.828 |
| Satellites | 16627547 | 0.585 |
| Flanks | 10220302 | 0.360 |
| Tandem Repeats | 6553662 | 0.231 |
| Gaps | 6172312 | 0.217 |
| HG002Q100v1.1 errors | 5039842 | 0.177 |
| VDJ | 3105977 | 0.109 |
| HG002-mosaic | 1178158 | 0.041 |
| TSPY2-segdups | 589097 | 0.021 |
| dipcall-pav_discrep-smvar | 333608 | 0.012 |
| Consecutive SVs | 251531 | 0.009 |
| dipcall-pav_discrep-stvar | 127132 | 0.004 |
| Dipcall-bugs T2TACE | 45604 | 0.002 |
| Self-Discrepancies | 15177 | 0.001 |

### D. Out of scope — do not attempt to (re)calculate these

- Author names/affiliations, NPS system upload status — administrative, not data.
- Live URL check for all hyperlinks — separate task, no calculation involved.
- NIST Disclaimer placeholder text — boilerplate insertion, not a calculated value.
- `guppy_version` and `snakemake_archive_doi` placeholders (see Section A) — external provenance/deposit status, not derivable from this repo.
- Externally-submitted sequencing depth/coverage claims in Methods (e.g., "≈ 30×" PacBio downsampling, "80×" ONT coverage, "35×" Illumina coverage, "≈ 1000 bp" Element insert size) — self-reported by contributing groups for their own submitted datasets; not computed by this pipeline. Do not attempt to verify by re-analyzing raw BAMs/FASTQs.
- Software/tool version numbers throughout Methods (e.g., `DRAGEN 4.4.x`, `Sniffles2 v2.5.3`, `pbmm2 v1.16.0`) — provenance facts, not calculated results.
- "DRAGEN has a 30 % increase in recall" and "pangenome reduces errors by 5 %" in Discussion — these are citations to *other* published work (refs 40, 24), not results from this pipeline. Feel free to sanity-check that the citation is described accurately, but don't try to recompute them from this repo's data.

## Deliverables

### 1. Notebook (primary human-review output)

A single, fully-rendered Quarto or Jupyter notebook (your choice — match this repo's convention and use `.qmd` with `_notebook_setup.R` if you're doing this in R, consistent with existing analysis notebooks) that:

- Walks through **every item in Sections A, B, and C above**, one at a time, showing:
  - The exact question being answered / claim being checked.
  - The code used to calculate it, referencing which existing notebook/data file it's built on (or explicitly noting where you had to write new code and why no existing computation covered it).
  - The calculated result.
- Ends with a **single summary table** with (at minimum) these columns:
  | Value ID | Manuscript location | Current text value | Calculated value | Match? (Yes/No/Needs review) | Notes |
  Sort this table so **mismatches and `##` placeholders appear first**, followed by confirmed matches.
- Clearly separates "placeholder — newly calculated" rows from "existing value — confirmed or flagged as inconsistent" rows (e.g., a status column or grouping).

### 2. JSON key-value file (for programmatic substitution)

A single JSON file, e.g. `docs/agent_work/manuscript_value_updates.json`, structured as an array of objects — one per value — with this schema:

```json
{
  "id": "sv_difficult_context_pct_increase",
  "section": "Results — Benchmark content and genome coverage",
  "current_text_snippet": "The v5.0 SV benchmark set includes ## % more variants in difficult genomic context",
  "current_value": null,
  "calculated_value": 63.4,
  "unit": "%",
  "status": "placeholder_resolved",
  "confidence": "high",
  "needs_human_review": false,
  "source": "benchmark_difficult.qmd — Variant Fold Change by Genomic Context, stvar, v5.0q(GRCh37) vs v0.6(GRCh37)",
  "notes": ""
}
```

Field notes:
- `current_text_snippet`: enough surrounding text (10-20 words) to unambiguously locate the value with a text search in the manuscript, since paragraph numbers shift between drafts.
- `current_value`: the value currently in the text, if any (`null` for `##` placeholders); for numeric text values pull out just the number (e.g. `699449` not `"≈ 600 000"`).
- `calculated_value`: plain number (not NIST-formatted — that's applied at substitution time).
- `status`: one of `placeholder_resolved`, `confirmed_match`, `inconsistency_found`, `needs_human_review`, `out_of_scope`.
- `needs_human_review`: `true` for anything like the callset-count ambiguity in Section B, item 3.
- Include one entry for **every row** in Sections A, B, and C's tables (including the ones you confirm match — don't only report mismatches in the JSON; the notebook summary table can lead with mismatches, but the JSON should be complete so it can be used as a full substitution manifest).

## Definition of done

- [ ] Every `##` placeholder in Section A has a `calculated_value` or an explicit `needs_human_review`/`out_of_scope` status with reasoning.
- [ ] All five items in Section B have been investigated and reported with reasoning (not just recomputed in isolation).
- [ ] All rows in Section C have a confirmed or flagged status.
- [ ] Notebook renders cleanly end-to-end with no errors.
- [ ] JSON validates as an array of objects matching the schema above, one entry per inventory row.
- [ ] Neither deliverable modifies any `.docx` file.
