# Test Suite

This directory contains R and Python tests for data loading logic, workflow helper logic, and selected workflow scripts.

## Structure

```text
tests/
├── conftest.py
├── fixtures/
│   └── sample.bed
├── test_cache.R
├── test_data_loading.R
├── test_exclusion_loading.R
├── test_schema_update.R
├── unit/
│   ├── test_common_helpers.py
│   ├── test_find_chr8_inversion.py
│   └── test_variant_binning.R
└── README.md
```

## Run Tests

From repository root.

## R tests

```bash
Rscript -e 'testthat::test_dir("tests")'
```

Run a single file:

```bash
Rscript -e 'testthat::test_file("tests/test_cache.R")'
```

## Python tests

```bash
pytest tests/ -v
```

Run a single file:

```bash
pytest tests/unit/test_common_helpers.py -v
```

## Workflow-level validation

Use Makefile targets for lint + formatting checks + workflow dry-run:

```bash
make test
```

## Test Data

- Small local fixture(s): `tests/fixtures/`
- Debug/integration-like workflow fixtures are configured in `config/config.test_grch38_debug.yaml`

Example dry-run with test config:

```bash
snakemake -s workflow/Snakefile \
  --replace-workflow-config \
  --configfile config/config.test_grch38_debug.yaml \
  -n results/variant_tables/v5.0q_GRCh38_smvar/variants.parquet
```

## Notes

- Cache-related R tests use temporary directories to avoid mutating the real cache.
- Pytest fixtures are defined in `tests/conftest.py`.
