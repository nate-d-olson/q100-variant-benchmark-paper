# Q100 Variant Benchmark - Test Suite

This directory contains tests for the Q100 variant benchmark pipeline.

## Structure

```
tests/
├── test_cache.R             # R tests for schemas and caching (45 tests)
├── test_data_loading.R      # R tests for data loading functions
├── unit/                    # Python unit tests for individual functions/modules
│   ├── test_common_helpers.py      # Tests for workflow/rules/common.smk
│   └── test_stratify_comparison.py # Tests for stratify_comparison.py
├── integration/             # Integration tests for workflow rules
├── fixtures/                # Test data files
├── conftest.py             # Pytest configuration and shared fixtures
└── README.md               # This file
```

## Running Tests

### R Tests

R tests use `testthat` and can be run from the repository root:

```bash
# Schema registry and Parquet caching tests (45 tests, self-contained with synthetic data)
Rscript -e 'testthat::test_file("tests/test_cache.R")'

# Data loading function tests (requires pipeline outputs in results/)
Rscript -e 'testthat::test_file("tests/test_data_loading.R")'
```

The cache tests use `withr::local_options(q100.cache_dir = tempdir)` for isolation -- they create temporary directories and never touch the real cache.

### Python Tests

#### Prerequisites

Install test dependencies:

```bash
pip install pytest pytest-cov
```

#### Run All Tests

```bash
# From repository root
pytest tests/ -v

# With coverage report
pytest tests/ -v --cov=workflow/scripts --cov=workflow/rules
```

#### Run Specific Test Files

```bash
# Unit tests only
pytest tests/unit/ -v

# Specific test file
pytest tests/unit/test_common_helpers.py -v

# Specific test class or function
pytest tests/unit/test_common_helpers.py::TestExclusionHelpers -v
pytest tests/unit/test_common_helpers.py::TestExclusionHelpers::test_get_exclusion_file_path -v
```

## Writing Tests

### Unit Tests

Unit tests should:

- Test individual functions in isolation
- Use mocks for external dependencies
- Be fast (<1s per test)
- Have clear assertions

Example:

```python
def test_get_exclusion_file_path():
    """Test exclusion file path generation."""
    result = get_exclusion_file_path("benchmark", "satellites", 0)
    expected = "resources/exclusions/benchmark/satellites_0.bed"
    assert result == expected
```

### Integration Tests

Integration tests should:

- Test complete workflows end-to-end
- Use real (but minimal) test data
- Verify expected outputs are created
- Clean up temporary files

## Test Data

Test fixtures are stored in `tests/fixtures/`. Keep test data:

- Small (prefer <1KB files when possible)
- Representative (cover edge cases)
- Well-documented (add comments explaining what each file tests)

### GRCh38 debug subset fixture

For integration-style debugging against realistic benchmark content, use:

- `tests/fixtures/grch38_debug_subset/`

This fixture includes clipped benchmark VCF/BED files, exclusion BEDs, genomic
context BEDs, and a selected region set (5 autosomes + X/Y) with both
overlap-rich regions and a no-overlap control region.

Regenerate it with:

```bash
python scripts/create_grch38_debug_subset.py --force
```

Run Snakemake against the subset fixture config with workflow replacement:

```bash
snakemake -s workflow/Snakefile \
  --replace-workflow-config \
  --configfile config/config.test_grch38_debug.yaml \
  -n results/variant_tables/v5q_GRCh38_smvar/variants.parquet
```

### Local File Support

The pipeline now supports local file paths in configuration files for test fixtures:

**Supported formats:**
- `file://` URLs: `file:///absolute/path/to/file.vcf.gz` or `file://$PWD/relative/path/file.bed`
- Absolute paths: `/absolute/path/to/file.vcf.gz`
- Relative paths: `tests/fixtures/file.bed` (converted to absolute internally)

**How it works:**
- Download rules automatically detect local file paths vs. remote URLs
- Local files are copied using `cp` instead of `wget`
- Checksum validation still performed (SHA256/MD5 as configured)
- Environment variable expansion supported (e.g., `$PWD`, `$HOME`)

**Example test config:**
```yaml
benchmarksets:
  v5q_GRCh38_smvar:
    vcf:
      url: "file://$PWD/tests/fixtures/grch38_debug_subset/benchmarksets/v5.0q_GRCh38_smvar_benchmark.vcf.gz"
      sha256: "b6ad3c3956ec1270a159e1f3bb3ca0d589c3c5b6d0338ab207947cf5e70628bf"
```

See `config/config.test_grch38_debug.yaml` for a complete example.

## Continuous Integration

Tests run automatically on:

- Pull requests to main branch
- Pushes to main branch

See `.github/workflows/tests.yml` for CI configuration.

## Coverage Goals

Target >80% code coverage for:

- `workflow/scripts/*.py` - Python analysis scripts
- `workflow/rules/common.smk` - Helper functions

Current coverage: Run `pytest --cov` to see current coverage report.
