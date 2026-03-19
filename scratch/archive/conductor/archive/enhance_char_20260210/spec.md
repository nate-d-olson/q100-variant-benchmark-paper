# Track Specification: Enhance Benchmark Characterization Analysis

## Overview

This track focuses on refining the benchmark characterization analysis to support the scientific manuscript. The goal is to provide granular insights into variant distributions, especially for small variants and exclusion regions.

## Goals

1. **Refine Variant Caching Schema:** Update the data schema to support filtering variants by precise length bins (e.g., <15bp, 15-49bp).
2. **Enhance Exclusion Analysis:** Update `benchmarkset_characterization.qmd` to include variant counts within exclusion regions.
3. **Fix Metrics:** Correct the calculation of "difficult region bases excluded" and "difficult regions included".
4. **Visualize Small Variants:** Generate breakdown tables/plots for small variants (<15bp vs 15-49bp).

## Requirements

### Data Schema

- **Variant Cache:** Must include a `variant_length` or `is_small_variant` (boolean) and `variant_type` field that allows efficient filtering.
- **Exclusion Tables:** Must support joining variant counts with exclusion regions.

### Analysis Notebook (`benchmarkset_characterization.qmd`)

- **Input:** Load cached variant tables and exclusion BED metrics.
- **Output:**
  - Table: Variant counts in exclusion regions.
  - Table: Small variant counts breakdown (<15bp, 15-49bp).
  - Fixed Metrics: Correct values for difficult region overlaps.

## Acceptance Criteria

- [ ] `benchmarkset_characterization.qmd` renders without error.
- [ ] Output tables correctly reflect the new small variant bins.
- [ ] Exclusion analysis includes variant counts, not just base pairs.
- [ ] All tests for data loading and schema validation pass.
