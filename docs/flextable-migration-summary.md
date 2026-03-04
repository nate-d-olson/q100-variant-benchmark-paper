---
title: Flextable Migration Summary
created: 2026-03-04
tags:
  - migration
  - r
  - quarto
  - tables
  - flextable
  - gt
---

## Flextable Migration Summary

## Goal

Migrate table workflows from `gt`-only assumptions to a `flextable`-compatible path that works in Word and Google Docs, while preserving reproducible Quarto rendering.

## Decision Outcome

- **Primary package for manuscript-friendly tables:** `flextable` + `officer`
- **`gt` remains useful for HTML-first outputs**
- For mixed-output notebooks, explicitly handle HTML-only table engines in non-HTML targets.

## What Was Implemented

1. Added/updated shared table theming and formatting helpers in `R/plot_themes.R`:
   - `theme_flextable_manuscript()`
   - `fmt_integer_flextable()`
   - `fmt_number_flextable()`
   - `cols_label_flextable()`
2. Added migration validation notebooks/scripts in `tests/`:
   - `tests/table_compatibility_test.qmd`
   - `tests/flextable_theme_implementation.qmd`
   - `tests/run_table_tests.R`
3. Removed migration sprawl documents/scripts and consolidated into this single summary.

## Critical Fixes Applied During Validation

### 1) Border namespace fix

- **Issue:** `flextable::fp_border()` is not exported.
- **Fix:** use `officer::fp_border()` in flextable border calls.

### 2) Word vertical alignment compatibility

- **Issue:** Word render failed with `vertical.align must be one of 'top', 'center', 'bottom'`.
- **Fix:** replace `valign = "middle"` with `valign = "center"`.

### 3) HTML-only chunks in DOCX targets

- **Issue:** DOCX renders failed when chunks produced HTML-only output.
- **Fix:** gate those chunks with `eval = knitr::is_html_output()`.

### 4) `gt` theme helper bug

- **Issue:** `theme_gt_manuscript()` used piped `do.call(...)` incorrectly.
- **Fix:** call `do.call(gt::tab_options, c(list(gt_object), merged_options))` outside the pipe.

### 5) Notebook syntax cleanups

- Removed malformed repeated namespace chains (e.g. `pkg::pkg::pkg::fn`) introduced during refactor.
- Corrected helper pipeline code that used `...` in an invalid expression context.

## Validation Evidence

The following validations were executed successfully:

- `Rscript tests/run_table_tests.R`
  - `tests/table_compatibility_test.qmd` render: ✅
  - `tests/flextable_theme_implementation.qmd` render: ✅
- Direct Quarto renders:
  - `quarto render tests/table_compatibility_test.qmd --to html`: ✅
  - `quarto render tests/table_compatibility_test.qmd --to docx`: ✅
  - `quarto render tests/flextable_theme_implementation.qmd --to html`: ✅
  - `quarto render tests/flextable_theme_implementation.qmd --to docx`: ✅
- Production helper smoke test:
  - Minimal `flextable` export to DOCX using `theme_flextable_manuscript()`: ✅

## Lessons Learned

1. **Word rendering paths are stricter than HTML paths**—validate both early.
2. **Namespace correctness matters** for R table stacks (`flextable` + `officer`).
3. **Quarto multi-format notebooks need explicit chunk gating** when packages output format-specific objects.
4. **Consolidated documentation reduces maintenance burden**; avoid creating many one-off migration docs.
5. **Keep generated artifacts out of commits** (`*.docx`, notebook `_files/` assets, local compiler config files).

## Recommended Ongoing Pattern

- Keep manuscript-style table helpers centralized in `R/plot_themes.R`.
- Use `tests/run_table_tests.R` as a quick migration regression check.
- For any notebook targeting both HTML and DOCX:
  - ensure flextable code is DOCX-safe,
  - gate HTML-only chunks with `knitr::is_html_output()`.

## Follow-up (optional)

- Convert one production analysis notebook as a pilot and compare rendered HTML vs DOCX outputs side-by-side before full migration.
