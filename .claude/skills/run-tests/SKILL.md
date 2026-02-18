---
name: run-tests
description: Run the full test suite — R testthat and Python pytest — reporting only failures
disable-model-invocation: true
---

# Run Tests

Run all unit tests for R and Python code.

## Steps

1. **R cache/schema tests** (most reliable — no pipeline outputs required):
   ```
   Rscript -e 'testthat::test_file("tests/test_cache.R")'
   Rscript -e 'testthat::test_file("tests/test_schema_update.R")'
   ```

2. **R data loading tests** (has known pre-existing failures without `results/`):
   ```
   Rscript -e 'testthat::test_file("tests/test_data_loading.R")'
   Rscript -e 'testthat::test_file("tests/test_exclusion_loading.R")'
   ```

3. **Python unit tests**:
   ```
   micromamba run -n q100-smk pytest tests/unit/ -v --tb=short
   ```

## Known Pre-existing Failures (not caused by recent changes)

- `parse_benchmark_id()` regex captures `"5.0q"` not `"v5.0q"` — tests expect `"v"` prefix
- `test_data_loading.R` references old column names (`strat_name`, `strat_filter`) from pre-refactor
- Any test that loads real data will fail without `results/` directory (pipeline outputs not in repo)

## Reporting

Report results grouped by test file. Distinguish between:
- **New failures**: Likely caused by recent changes — investigate these
- **Pre-existing failures**: Listed above — note them but don't treat as regressions
- **Passes**: Confirm count
