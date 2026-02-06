#!/usr/bin/env python3
"""
Centralized logging configuration for Q100 variant benchmark pipeline scripts.

Provides consistent logging format and behavior across all Python processing scripts.
"""

import logging
import sys
from pathlib import Path
from typing import Optional


def setup_logger(
    name: str, log_file: Optional[Path] = None, level: int = logging.INFO
) -> logging.Logger:
    """
    Configure and return a logger with consistent formatting.

    Args:
        name: Logger name (typically __name__ from calling module)
        log_file: Optional path to log file. If None, logs only to stderr
        level: Logging level (default: INFO)

    Returns:
        Configured logger instance

    Example:
        >>> logger = setup_logger(__name__)
        >>> logger.info("Processing file", extra={"file": "input.vcf"})
    """
    logger = logging.getLogger(name)
    logger.setLevel(level)

    # Avoid duplicate handlers if logger already configured
    if logger.handlers:
        return logger

    # Consistent format: [TIMESTAMP] [LEVEL] [MODULE] Message
    formatter = logging.Formatter(
        fmt="[%(asctime)s] [%(levelname)s] [%(name)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # Always log to stderr (Snakemake captures this)
    stderr_handler = logging.StreamHandler(sys.stderr)
    stderr_handler.setLevel(level)
    stderr_handler.setFormatter(formatter)
    logger.addHandler(stderr_handler)

    # Optional file logging
    if log_file:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(level)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger


def log_context(**kwargs) -> str:
    """
    Format context dictionary for structured logging.

    Args:
        **kwargs: Key-value pairs to format

    Returns:
        Formatted string with context information

    Example:
        >>> context = log_context(file="input.vcf", variants=1000, time_s=12.5)
        >>> logger.info(f"Processing complete: {context}")
    """
    if not kwargs:
        return ""

    items = [f"{k}={v}" for k, v in kwargs.items()]
    return ", ".join(items)
