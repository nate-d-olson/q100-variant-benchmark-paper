#!/usr/bin/env python3
"""
Custom exception classes for Q100 variant benchmark pipeline.

Provides context-rich error handling for different failure modes.
"""

from pathlib import Path
from typing import Optional, List


class PipelineError(Exception):
    """Base exception for all pipeline errors."""

    pass


class ValidationError(PipelineError):
    """Raised when input data validation fails."""

    def __init__(
        self,
        message: str,
        file_path: Optional[Path] = None,
        line_number: Optional[int] = None,
        expected: Optional[str] = None,
        actual: Optional[str] = None,
    ):
        """
        Initialize validation error with context.

        Args:
            message: Error description
            file_path: Path to file that failed validation
            line_number: Line number where error occurred
            expected: Expected format or value
            actual: Actual value encountered
        """
        self.file_path = file_path
        self.line_number = line_number
        self.expected = expected
        self.actual = actual

        # Build detailed error message
        error_parts = [message]

        if file_path:
            error_parts.append(f"File: {file_path}")

        if line_number is not None:
            error_parts.append(f"Line: {line_number}")

        if expected and actual:
            error_parts.append(f"Expected: {expected}")
            error_parts.append(f"Actual: {actual}")

        super().__init__("\n  ".join(error_parts))


class DataFormatError(PipelineError):
    """Raised when data format is incorrect or malformed."""

    def __init__(
        self,
        message: str,
        file_path: Optional[Path] = None,
        format_type: Optional[str] = None,
        required_columns: Optional[List[str]] = None,
        missing_columns: Optional[List[str]] = None,
    ):
        """
        Initialize data format error with context.

        Args:
            message: Error description
            file_path: Path to file with format error
            format_type: Expected format (e.g., "VCF", "BED", "TSV")
            required_columns: List of required columns
            missing_columns: List of missing columns
        """
        self.file_path = file_path
        self.format_type = format_type
        self.required_columns = required_columns
        self.missing_columns = missing_columns

        error_parts = [message]

        if file_path:
            error_parts.append(f"File: {file_path}")

        if format_type:
            error_parts.append(f"Expected format: {format_type}")

        if required_columns:
            error_parts.append(f"Required columns: {', '.join(required_columns)}")

        if missing_columns:
            error_parts.append(f"Missing columns: {', '.join(missing_columns)}")
            error_parts.append(
                f"Suggestion: Check that input file has correct format and headers"
            )

        super().__init__("\n  ".join(error_parts))


class ProcessingError(PipelineError):
    """Raised when data processing fails."""

    def __init__(
        self,
        message: str,
        operation: Optional[str] = None,
        file_path: Optional[Path] = None,
        context: Optional[dict] = None,
    ):
        """
        Initialize processing error with context.

        Args:
            message: Error description
            operation: Operation that failed (e.g., "counting variants")
            file_path: File being processed
            context: Additional context dictionary
        """
        self.operation = operation
        self.file_path = file_path
        self.context = context

        error_parts = [message]

        if operation:
            error_parts.append(f"Operation: {operation}")

        if file_path:
            error_parts.append(f"File: {file_path}")

        if context:
            context_str = ", ".join(f"{k}={v}" for k, v in context.items())
            error_parts.append(f"Context: {context_str}")

        super().__init__("\n  ".join(error_parts))


class ConfigurationError(PipelineError):
    """Raised when configuration is invalid or missing."""

    def __init__(
        self,
        message: str,
        parameter: Optional[str] = None,
        expected_type: Optional[str] = None,
        suggestion: Optional[str] = None,
    ):
        """
        Initialize configuration error with context.

        Args:
            message: Error description
            parameter: Configuration parameter that caused error
            expected_type: Expected type or format
            suggestion: Suggestion for fixing the error
        """
        self.parameter = parameter
        self.expected_type = expected_type
        self.suggestion = suggestion

        error_parts = [message]

        if parameter:
            error_parts.append(f"Parameter: {parameter}")

        if expected_type:
            error_parts.append(f"Expected: {expected_type}")

        if suggestion:
            error_parts.append(f"Suggestion: {suggestion}")

        super().__init__("\n  ".join(error_parts))
