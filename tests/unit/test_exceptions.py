#!/usr/bin/env python3
"""
Unit tests for custom exception classes.
"""

import pytest
from pathlib import Path
import sys

# Add workflow/scripts to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "workflow" / "scripts"))

from exceptions import (
    PipelineError,
    ValidationError,
    DataFormatError,
    ProcessingError,
    ConfigurationError,
)


class TestPipelineError:
    """Tests for base PipelineError class."""

    def test_basic_exception(self):
        """Test basic exception creation."""
        err = PipelineError("Test error")
        assert str(err) == "Test error"
        assert isinstance(err, Exception)


class TestValidationError:
    """Tests for ValidationError class."""

    def test_simple_validation_error(self):
        """Test simple validation error."""
        err = ValidationError("File validation failed")
        assert "File validation failed" in str(err)

    def test_validation_error_with_file(self):
        """Test validation error with file path."""
        err = ValidationError(
            "Invalid format", file_path=Path("/path/to/file.txt"), line_number=42
        )
        assert "Invalid format" in str(err)
        assert "/path/to/file.txt" in str(err)
        assert "Line: 42" in str(err)

    def test_validation_error_with_expected_actual(self):
        """Test validation error with expected vs actual."""
        err = ValidationError(
            "Value mismatch", expected="3 columns", actual="2 columns"
        )
        assert "Value mismatch" in str(err)
        assert "Expected: 3 columns" in str(err)
        assert "Actual: 2 columns" in str(err)

    def test_validation_error_full_context(self):
        """Test validation error with full context."""
        err = ValidationError(
            "Coordinate error",
            file_path=Path("test.bed"),
            line_number=10,
            expected="start < end",
            actual="start >= end",
        )
        error_str = str(err)
        assert "Coordinate error" in error_str
        assert "test.bed" in error_str
        assert "Line: 10" in error_str
        assert "Expected:" in error_str
        assert "Actual:" in error_str


class TestDataFormatError:
    """Tests for DataFormatError class."""

    def test_simple_format_error(self):
        """Test simple format error."""
        err = DataFormatError("Invalid BED format")
        assert "Invalid BED format" in str(err)

    def test_format_error_with_file_type(self):
        """Test format error with file and format type."""
        err = DataFormatError(
            "Missing columns", file_path=Path("data.tsv"), format_type="TSV"
        )
        assert "Missing columns" in str(err)
        assert "data.tsv" in str(err)
        assert "Expected format: TSV" in str(err)

    def test_format_error_with_missing_columns(self):
        """Test format error with column details."""
        err = DataFormatError(
            "Column validation failed",
            file_path=Path("test.vcf"),
            format_type="VCF",
            required_columns=["CHROM", "POS", "REF"],
            missing_columns=["REF"],
        )
        error_str = str(err)
        assert "Column validation failed" in error_str
        assert "Required columns: CHROM, POS, REF" in error_str
        assert "Missing columns: REF" in error_str
        assert "Suggestion:" in error_str


class TestProcessingError:
    """Tests for ProcessingError class."""

    def test_simple_processing_error(self):
        """Test simple processing error."""
        err = ProcessingError("Processing failed")
        assert "Processing failed" in str(err)

    def test_processing_error_with_operation(self):
        """Test processing error with operation."""
        err = ProcessingError("Failed to parse data", operation="variant counting")
        assert "Failed to parse data" in str(err)
        assert "Operation: variant counting" in str(err)

    def test_processing_error_with_context(self):
        """Test processing error with context dict."""
        err = ProcessingError(
            "Processing error",
            operation="annotation expansion",
            file_path=Path("input.tsv"),
            context={"row": 150, "column": "STRAT_IDS"},
        )
        error_str = str(err)
        assert "Processing error" in error_str
        assert "Operation: annotation expansion" in error_str
        assert "File: input.tsv" in error_str
        assert "Context:" in error_str
        assert "row=150" in error_str
        assert "column=STRAT_IDS" in error_str


class TestConfigurationError:
    """Tests for ConfigurationError class."""

    def test_simple_config_error(self):
        """Test simple configuration error."""
        err = ConfigurationError("Invalid configuration")
        assert "Invalid configuration" in str(err)

    def test_config_error_with_parameter(self):
        """Test configuration error with parameter."""
        err = ConfigurationError(
            "Missing required parameter", parameter="benchmark.ref"
        )
        assert "Missing required parameter" in str(err)
        assert "Parameter: benchmark.ref" in str(err)

    def test_config_error_with_expected_type(self):
        """Test configuration error with expected type."""
        err = ConfigurationError(
            "Type mismatch", parameter="cores", expected_type="integer"
        )
        assert "Type mismatch" in str(err)
        assert "Parameter: cores" in str(err)
        assert "Expected: integer" in str(err)

    def test_config_error_with_suggestion(self):
        """Test configuration error with suggestion."""
        err = ConfigurationError(
            "Invalid value",
            parameter="reference",
            expected_type="GRCh37, GRCh38, or CHM13",
            suggestion="Check config/config.yaml for valid reference names",
        )
        error_str = str(err)
        assert "Invalid value" in error_str
        assert "Suggestion: Check config/config.yaml" in error_str


class TestExceptionInheritance:
    """Tests for exception inheritance."""

    def test_all_exceptions_inherit_from_pipeline_error(self):
        """Test that all custom exceptions inherit from PipelineError."""
        assert issubclass(ValidationError, PipelineError)
        assert issubclass(DataFormatError, PipelineError)
        assert issubclass(ProcessingError, PipelineError)
        assert issubclass(ConfigurationError, PipelineError)

    def test_all_exceptions_inherit_from_exception(self):
        """Test that all custom exceptions inherit from Exception."""
        assert issubclass(PipelineError, Exception)
        assert issubclass(ValidationError, Exception)
        assert issubclass(DataFormatError, Exception)
        assert issubclass(ProcessingError, Exception)
        assert issubclass(ConfigurationError, Exception)

    def test_exception_catching(self):
        """Test that exceptions can be caught as PipelineError."""
        # Should be catchable as specific type
        with pytest.raises(ValidationError):
            raise ValidationError("test")

        # Should also be catchable as PipelineError
        with pytest.raises(PipelineError):
            raise ValidationError("test")

        # Should also be catchable as Exception
        with pytest.raises(Exception):
            raise ValidationError("test")
