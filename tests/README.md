# Q100 Variant Benchmark - Test Suite

This directory contains tests for the Q100 variant benchmark pipeline.

## Structure

```
tests/
├── unit/                    # Unit tests for individual functions/modules
│   ├── test_common_helpers.py      # Tests for workflow/rules/common.smk
│   └── test_stratify_comparison.py # Tests for stratify_comparison.py
├── integration/             # Integration tests for workflow rules
├── fixtures/                # Test data files
├── conftest.py             # Pytest configuration and shared fixtures
└── README.md               # This file
```

## Running Tests

### Prerequisites

Install test dependencies:

```bash
pip install pytest pytest-cov
```

### Run All Tests

```bash
# From repository root
pytest tests/ -v

# With coverage report
pytest tests/ -v --cov=workflow/scripts --cov=workflow/rules
```

### Run Specific Test Files

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
