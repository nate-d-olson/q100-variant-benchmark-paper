# Codebase Review Findings (Feb 2026)

## Snakemake Workflow Improvements

- [ ] **Remove dead code**: remove rule, functions, and variables not used in current codebase.
- [ ] **Re-organize** snakemake modules to have a more logical structure grouping of rules. Consolidate rule files when appropriate.
- [ ] **Decouple from Human-Specific Logic**: Move hardcoded chromosome lists (`1-22, X, Y`) from `common.smk` to `config.yaml` for better portability.
- [ ] **Streamline Rules**: Evaluate custom Python scripts (e.g., `extract_info_fields.py`) for potential replacement with standard tools like `bcftools query`.
- [ ] **Refactor `common.smk`**: Split the large helper file into separate modules for configuration, path generation, and functional logic.

### R Analysis & Data Loading

- [ ] **Robust Metadata Handling**: Transition from regex-based filename parsing (`parse_benchmark_id`) to structured metadata sidecars (JSON/YAML) for pipeline outputs.
- [ ] **Simplify Data Loading**: Evaluate replacing the custom caching/validation logic in `R/data_loading.R` with a standard framework like the R `targets` package.
- [ ] **Documentation Debt**: Address internal TODOs and "documentation debt" found within `docs/` files and notebook comments.


## Analysis Revisions

### benchmarkset_characterization.qmd

- Add variants to the exclusion tables
- Small variant counts breakdown by <15bp , 15 - 49bp
- Look into number and distribution of small intervals (50 - 100bp) in v5q

### benchmark_unique_regions.qmd

- Validate summary outputs against old_only status CSVs
- Add compact figure for v5-only vs previous-only base/variant deltas


### manuscript/results.qmd

- **Figure 2**: Small and structural variant counts, coverage, and comparison to previous benchmark sets. %%TODO%% Stratified, annotation based???

### analysis/external_evaluation.qmd

- **TODO**: linear model to get confidence intervals for curations

### Low priority

- more sensible results and log directory structures
- cleanup data_loading script
- look into intervals, less than 100bp (potential bug in exclusions) in small variant benchmark regions bed files for v5 benchmarkset.

## Other

- Annotating V5 benchmark vcfs based on exclusions and add filter field values for variants outside the benchmark regions.
- Updating difficult variant counts to include complex variants
