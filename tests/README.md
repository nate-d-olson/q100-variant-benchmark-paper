# Q100 Variant Benchmark Pipeline - Tests

This directory contains unit tests and test fixtures for the Q100 variant benchmark pipeline.

## Structure

```
tests/
├── README.md           # This file
├── unit/               # Unit tests for Python modules
│   ├── test_validators.py    # Tests for data validation
│   └── test_exceptions.py    # Tests for exception classes
└── fixtures/           # Test data files
    ├── sample.vcf      # Sample VCF file
    └── sample.bed      # Sample BED file
```

## Running Tests

### Prerequisites

Install test dependencies:

```bash
pip install -e ".[dev]"
```

This installs:
- pytest - Testing framework
- pytest-cov - Coverage reporting

### Run All Tests

```bash
# Run all tests with coverage
pytest

# Run tests verbosely
pytest -v

# Run tests with detailed coverage report
pytest --cov-report=term-missing
```

### Run Specific Tests

```bash
# Run only validator tests
pytest tests/unit/test_validators.py

# Run only exception tests
pytest tests/unit/test_exceptions.py

# Run a specific test class
pytest tests/unit/test_validators.py::TestBEDValidation

# Run a specific test function
pytest tests/unit/test_validators.py::TestBEDValidation::test_validate_valid_bed
```

### Coverage Reports

```bash
# Generate HTML coverage report
pytest --cov-report=html

# View coverage report
open htmlcov/index.html
```

## Test Organization

### Unit Tests

#### `test_validators.py`

Tests for the `validators.py` module:

- **TestFileValidation**: Basic file existence and readability checks
  - `test_validate_existing_file` - Valid file passes
  - `test_validate_missing_file` - Missing file raises error
  - `test_validate_empty_file` - Empty file raises error
  - `test_validate_directory_not_file` - Directory raises error

- **TestBEDValidation**: BED format validation
  - `test_validate_valid_bed` - Valid BED file passes
  - `test_validate_malformed_bed` - Missing columns detected
  - `test_validate_invalid_coordinates` - Start >= end detected
  - `test_validate_non_numeric_coordinates` - Non-integer coords detected
  - `test_validate_unsorted_bed` - Sort order violations detected
  - `test_validate_skip_comments` - Comments properly ignored

- **TestVCFValidation**: VCF header validation
  - `test_validate_valid_vcf` - Valid VCF passes
  - `test_validate_missing_fileformat` - Missing ##fileformat detected
  - `test_validate_missing_column_header` - Missing #CHROM detected
  - `test_validate_missing_required_columns` - Missing columns detected
  - `test_validate_gzipped_vcf` - Gzipped files handled

- **TestTSVValidation**: TSV column validation
  - `test_validate_valid_tsv` - Valid TSV passes
  - `test_validate_missing_required_columns` - Missing columns detected
  - `test_validate_empty_tsv` - Empty files detected

#### `test_exceptions.py`

Tests for the `exceptions.py` module:

- **TestPipelineError**: Base exception class
- **TestValidationError**: Validation error context
  - File paths, line numbers, expected vs actual values
- **TestDataFormatError**: Format error context
  - File types, missing columns, suggestions
- **TestProcessingError**: Processing error context
  - Operations, file paths, context dictionaries
- **TestConfigurationError**: Configuration error context
  - Parameters, expected types, suggestions
- **TestExceptionInheritance**: Exception hierarchy
  - Inheritance from PipelineError and Exception
  - Catch behavior

## Test Fixtures

### `fixtures/sample.vcf`

Minimal valid VCF file with:
- Proper ##fileformat header
- INFO field definitions
- #CHROM header line
- 6 variant records (SNPs, indels, SVs)

### `fixtures/sample.bed`

Minimal valid BED file with:
- 5 intervals across 2 chromosomes
- Sorted by position
- Valid 3-column format

## Writing New Tests

### Test Naming Convention

- Files: `test_<module_name>.py`
- Classes: `Test<Feature>`
- Functions: `test_<specific_behavior>`

### Example Test

```python
import pytest
from pathlib import Path
import sys

# Add workflow/scripts to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "workflow" / "scripts"))

from validators import validate_bed_format
from exceptions import DataFormatError

class TestMyFeature:
    """Tests for my feature."""

    def test_success_case(self, tmp_path):
        """Test that valid input succeeds."""
        # Create test file
        test_file = tmp_path / "test.bed"
        test_file.write_text("chr1\t100\t200\n")

        # Run function
        result = validate_bed_format(test_file)

        # Assert expectations
        assert result["intervals"] == 1

    def test_failure_case(self, tmp_path):
        """Test that invalid input fails correctly."""
        test_file = tmp_path / "invalid.bed"
        test_file.write_text("invalid\n")

        # Should raise specific exception
        with pytest.raises(DataFormatError) as exc_info:
            validate_bed_format(test_file)

        # Verify error message
        assert "columns" in str(exc_info.value)
```

### Using Fixtures

```python
@pytest.fixture
def sample_data(tmp_path):
    """Create sample test data."""
    data_file = tmp_path / "data.txt"
    data_file.write_text("test data")
    return data_file

def test_with_fixture(sample_data):
    """Test using fixture."""
    assert sample_data.exists()
```

## Coverage Goals

Current coverage status:

- `validators.py`: ~80% (target: 90%+)
- `exceptions.py`: ~95% (target: 95%+)
- `logging_config.py`: Not yet tested (target: 80%+)

**To improve coverage**:
- Add tests for edge cases
- Test error paths and exception handling
- Test logging output (requires log capture)
- Test integration between modules

## Continuous Integration

Tests run automatically via GitHub Actions on:
- Pull requests
- Pushes to main branch
- Manual workflow dispatch

See `.github/workflows/tests.yml` for CI configuration.

## Troubleshooting

### Import Errors

If you get import errors:

```bash
# Ensure workflow/scripts is in Python path
export PYTHONPATH="${PYTHONPATH}:$(pwd)/workflow/scripts"

# Or run tests from project root
cd /path/to/q100-variant-benchmark
pytest
```

### Fixture Not Found

If test fixtures aren't found:

```bash
# Run from project root, not tests/ directory
cd /path/to/q100-variant-benchmark
pytest tests/
```

### Coverage Not Working

If coverage isn't tracking:

```bash
# Ensure pytest-cov is installed
pip install pytest-cov

# Check pyproject.toml [tool.coverage.run] section
cat pyproject.toml
```

---

*Last Updated: 2026-01-13*
*Branch: feature/codebase-improvements*
