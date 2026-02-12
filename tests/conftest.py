"""
Pytest configuration and shared fixtures for Q100 variant benchmark tests.
"""

import pytest
from pathlib import Path


@pytest.fixture
def test_data_dir():
    """Return path to test data directory."""
    return Path(__file__).parent / "fixtures"


@pytest.fixture
def sample_vcf(test_data_dir):
    """Path to sample VCF file for testing."""
    return test_data_dir / "sample.vcf.gz"


@pytest.fixture
def sample_bed(test_data_dir):
    """Path to sample BED file for testing."""
    return test_data_dir / "sample.bed"
