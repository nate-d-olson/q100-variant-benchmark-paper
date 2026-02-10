# Implementation Plan - Enhance Benchmark Characterization Analysis

## Phase 1: Schema Updates and Data Loading
This phase ensures the underlying data structures support the required analysis.

- [x] **Task: Update Variant Cache Schema** [51edaf3]
    - [x] Write test for new schema fields (variant length bins).
    - [x] Update `R/schemas.R` to include variant length/type definitions.
    - [x] Update `R/cache.R` (if needed) to ensure these fields are preserved.
- [x] **Task: Refactor Data Loading for Exclusions** [57ecf79]
    - [x] Sub-task: Write test for loading exclusion variant counts.
    - [x] Sub-task: Update `R/data_loading.R` to load `results/exclusions/{benchmark}/exclusion_impact.csv` with variant counts.
- [x] **Task: Conductor - User Manual Verification 'Phase 1' (Protocol in workflow.md)** [checkpoint: dacccf6]

## Phase 2: Notebook Refinement (Characterization)
This phase implements the analytical logic and visualizations in the Quarto notebook.

- [x] **Task: Implement Small Variant Breakdowns** [fc0335b]
    - [x] Sub-task: Write test for small variant binning logic (R unit test).
    - [x] Sub-task: Update `benchmarkset_characterization.qmd` to calculate and display <15bp vs 15-49bp counts.
- [x] **Task: Incorporate Exclusion Variant Counts** [38eedad]
    - [x] Sub-task: Update `benchmarkset_characterization.qmd` to join exclusion regions with variant counts.
    - [x] Sub-task: Generate summary table of variants removed by each exclusion.
- [x] **Task: Fix Difficult Region Metrics** [29c3b1d]
    - [x] Sub-task: Debug and fix the logic for "difficult regions included" vs "excluded".
    - [x] Sub-task: Verify against known ground truth or manual spot check (documented in notebook).
- [ ] **Task: Conductor - User Manual Verification 'Phase 2' (Protocol in workflow.md)**
