#!/usr/bin/env python3
"""
Unit tests for data validators.
"""

import pytest
from pathlib import Path
import sys
import tempfile
import gzip

# Add workflow/scripts to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "workflow" / "scripts"))

from validators import (
    validate_file_exists,
    validate_bed_format,
    validate_vcf_header,
    validate_tsv_columns,
)
from exceptions import ValidationError, DataFormatError


class TestFileValidation:
    """Tests for basic file validation."""

    def test_validate_existing_file(self, tmp_path):
        """Test validation of existing file."""
        test_file = tmp_path / "test.txt"
        test_file.write_text("content")

        # Should not raise
        validate_file_exists(test_file)

    def test_validate_missing_file(self, tmp_path):
        """Test validation fails for missing file."""
        test_file = tmp_path / "missing.txt"

        with pytest.raises(ValidationError) as exc_info:
            validate_file_exists(test_file)

        assert "not found" in str(exc_info.value)

    def test_validate_empty_file(self, tmp_path):
        """Test validation fails for empty file."""
        test_file = tmp_path / "empty.txt"
        test_file.write_text("")

        with pytest.raises(ValidationError) as exc_info:
            validate_file_exists(test_file)

        assert "empty" in str(exc_info.value)

    def test_validate_directory_not_file(self, tmp_path):
        """Test validation fails for directory."""
        test_dir = tmp_path / "testdir"
        test_dir.mkdir()

        with pytest.raises(ValidationError) as exc_info:
            validate_file_exists(test_dir)

        assert "not a regular file" in str(exc_info.value)


class TestBEDValidation:
    """Tests for BED format validation."""

    def test_validate_valid_bed(self):
        """Test validation of valid BED file."""
        fixtures_dir = Path(__file__).parent.parent / "fixtures"
        bed_file = fixtures_dir / "sample.bed"

        stats = validate_bed_format(bed_file, check_sorted=True)

        assert stats["intervals"] == 5
        assert stats["chromosomes"] == 2
        assert "chr1" in stats["chrom_names"]
        assert "chr2" in stats["chrom_names"]

    def test_validate_malformed_bed(self, tmp_path):
        """Test validation fails for malformed BED."""
        bed_file = tmp_path / "malformed.bed"
        bed_file.write_text("chr1\t100\n")  # Missing end column

        with pytest.raises(DataFormatError) as exc_info:
            validate_bed_format(bed_file)

        assert "Expected at least 3 columns" in str(exc_info.value)

    def test_validate_invalid_coordinates(self, tmp_path):
        """Test validation fails for invalid coordinates."""
        bed_file = tmp_path / "invalid.bed"
        bed_file.write_text("chr1\t200\t100\n")  # Start > end

        with pytest.raises(DataFormatError) as exc_info:
            validate_bed_format(bed_file)

        assert "Start" in str(exc_info.value) and ">=" in str(exc_info.value)

    def test_validate_non_numeric_coordinates(self, tmp_path):
        """Test validation fails for non-numeric coordinates."""
        bed_file = tmp_path / "non_numeric.bed"
        bed_file.write_text("chr1\tabc\t200\n")

        with pytest.raises(DataFormatError) as exc_info:
            validate_bed_format(bed_file)

        assert "not an integer" in str(exc_info.value)

    def test_validate_unsorted_bed(self, tmp_path):
        """Test validation detects unsorted BED."""
        bed_file = tmp_path / "unsorted.bed"
        bed_file.write_text("chr1\t100\t200\nchr1\t50\t100\n")  # Unsorted

        with pytest.raises(DataFormatError) as exc_info:
            validate_bed_format(bed_file, check_sorted=True)

        assert "not sorted" in str(exc_info.value)

    def test_validate_skip_comments(self, tmp_path):
        """Test validation skips comment lines."""
        bed_file = tmp_path / "comments.bed"
        bed_file.write_text("#comment\ntrack name=test\nchr1\t100\t200\n")

        stats = validate_bed_format(bed_file)
        assert stats["intervals"] == 1


class TestVCFValidation:
    """Tests for VCF header validation."""

    def test_validate_valid_vcf(self):
        """Test validation of valid VCF file."""
        fixtures_dir = Path(__file__).parent.parent / "fixtures"
        vcf_file = fixtures_dir / "sample.vcf"

        stats = validate_vcf_header(vcf_file)

        assert stats["has_fileformat"] is True
        assert stats["has_column_header"] is True
        assert stats["header_lines"] >= 1

    def test_validate_missing_fileformat(self, tmp_path):
        """Test validation fails without fileformat line."""
        vcf_file = tmp_path / "no_format.vcf"
        vcf_file.write_text("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        with pytest.raises(DataFormatError) as exc_info:
            validate_vcf_header(vcf_file)

        assert "missing ##fileformat" in str(exc_info.value)

    def test_validate_missing_column_header(self, tmp_path):
        """Test validation fails without column header."""
        vcf_file = tmp_path / "no_header.vcf"
        vcf_file.write_text("##fileformat=VCFv4.2\n")

        with pytest.raises(DataFormatError) as exc_info:
            validate_vcf_header(vcf_file)

        assert "missing #CHROM" in str(exc_info.value)

    def test_validate_missing_required_columns(self, tmp_path):
        """Test validation fails with missing required columns."""
        vcf_file = tmp_path / "missing_cols.vcf"
        vcf_file.write_text(
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\n"  # Missing REF, ALT, etc.
        )

        with pytest.raises(DataFormatError) as exc_info:
            validate_vcf_header(vcf_file)

        assert "missing required columns" in str(exc_info.value)

    def test_validate_gzipped_vcf(self, tmp_path):
        """Test validation works with gzipped VCF."""
        vcf_gz = tmp_path / "sample.vcf.gz"

        content = (
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        )

        with gzip.open(vcf_gz, 'wt') as f:
            f.write(content)

        stats = validate_vcf_header(vcf_gz)
        assert stats["has_fileformat"] is True


class TestTSVValidation:
    """Tests for TSV column validation."""

    def test_validate_valid_tsv(self, tmp_path):
        """Test validation of valid TSV file."""
        tsv_file = tmp_path / "test.tsv"
        tsv_file.write_text("col1\tcol2\tcol3\nval1\tval2\tval3\n")

        stats = validate_tsv_columns(
            tsv_file,
            required_columns=["col1", "col2"],
            optional_columns=["col3", "col4"]
        )

        assert stats["total_columns"] == 3
        assert stats["required_present"] == 2
        assert stats["optional_present"] == 1

    def test_validate_missing_required_columns(self, tmp_path):
        """Test validation fails with missing required columns."""
        tsv_file = tmp_path / "missing.tsv"
        tsv_file.write_text("col1\tcol2\n")

        with pytest.raises(DataFormatError) as exc_info:
            validate_tsv_columns(
                tsv_file,
                required_columns=["col1", "col2", "col3"]
            )

        assert "missing required columns" in str(exc_info.value)
        assert "col3" in str(exc_info.value)

    def test_validate_empty_tsv(self, tmp_path):
        """Test validation fails for TSV with no header."""
        tsv_file = tmp_path / "empty.tsv"
        tsv_file.write_text("")

        with pytest.raises(DataFormatError) as exc_info:
            validate_tsv_columns(tsv_file, required_columns=["col1"])

        assert "no header line" in str(exc_info.value)


# Test fixtures
@pytest.fixture
def sample_bed_file(tmp_path):
    """Create a sample BED file for testing."""
    bed_file = tmp_path / "test.bed"
    bed_file.write_text(
        "chr1\t100\t200\n"
        "chr1\t300\t400\n"
        "chr2\t100\t250\n"
    )
    return bed_file


@pytest.fixture
def sample_vcf_file(tmp_path):
    """Create a sample VCF file for testing."""
    vcf_file = tmp_path / "test.vcf"
    vcf_file.write_text(
        "##fileformat=VCFv4.2\n"
        "##INFO=<ID=TYPE,Number=A,Type=String,Description=\"Type\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t100\t.\tA\tG\t60\tPASS\tTYPE=snp\n"
    )
    return vcf_file
