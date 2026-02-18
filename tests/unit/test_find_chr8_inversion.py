"""
Unit tests for workflow/scripts/find_chr8_inversion.py.

Covers correct SyRI column parsing (Bug 3 fix: query coords are cols 6-7, not 5-6),
largest-inversion selection, and error paths.
"""

import json
import sys
import textwrap
from pathlib import Path

import pytest

# Add workflow/scripts to path so the module can be imported directly
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "workflow" / "scripts"))
from find_chr8_inversion import find_largest_inversion


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _syri_line(
    ref_chrom,
    ref_start,
    ref_end,
    qry_chrom,
    qry_start,
    qry_end,
    ann_type,
):
    """Build a 11-column tab-delimited SyRI .out line."""
    return "\t".join(
        [
            ref_chrom,
            str(ref_start),
            str(ref_end),
            "-",  # col 3 (unused by script)
            "-",  # col 4 (unused by script)
            qry_chrom,  # col 5 — qry chromosome (string, NOT a coordinate)
            str(qry_start),  # col 6
            str(qry_end),   # col 7
            "-",  # col 8
            "-",  # col 9
            ann_type,  # col 10
        ]
    )


@pytest.fixture
def tmp_syri(tmp_path):
    """Factory: writes lines to a temp .out file and returns its path."""

    def _write(lines):
        p = tmp_path / "test.syri.out"
        p.write_text("\n".join(lines) + "\n")
        return str(p)

    return _write


# ---------------------------------------------------------------------------
# Tests: column parsing (Bug 3)
# ---------------------------------------------------------------------------


class TestColumnParsing:
    """Verify that query coordinates come from columns 6-7, not 5-6."""

    def test_qry_coords_not_chrom_string(self, tmp_syri):
        """Column 5 is the qry chromosome string; parsing it as int must NOT happen."""
        # If the old bug were present, int("chr8") would raise ValueError.
        lines = [_syri_line("chr8", 1000, 5000, "chr8", 2000, 6000, "INV")]
        result = find_largest_inversion(tmp_syri(lines), "chr8")
        assert result["qry_start"] == 2000
        assert result["qry_end"] == 6000

    def test_ref_coords_correct(self, tmp_syri):
        """Reference coordinates come from columns 1-2 (0-indexed)."""
        lines = [_syri_line("chr8", 100, 900, "chr8", 200, 800, "INV")]
        result = find_largest_inversion(tmp_syri(lines), "chr8")
        assert result["ref_start"] == 100
        assert result["ref_end"] == 900

    def test_size_uses_qry_interval(self, tmp_syri):
        """Reported size is derived from the query interval (cols 6-7)."""
        lines = [_syri_line("chr8", 0, 1000, "chr8", 0, 4000, "INV")]
        result = find_largest_inversion(tmp_syri(lines), "chr8")
        assert result["size"] == 4000


# ---------------------------------------------------------------------------
# Tests: inversion selection
# ---------------------------------------------------------------------------


class TestInversionSelection:
    def test_picks_largest_by_query_size(self, tmp_syri):
        lines = [
            _syri_line("chr8", 1000, 2000, "chr8", 1000, 2000, "INV"),  # size 1000
            _syri_line("chr8", 5000, 15000, "chr8", 5000, 20000, "INV"),  # size 15000
            _syri_line("chr8", 3000, 6000, "chr8", 3000, 6000, "INV"),  # size 3000
        ]
        result = find_largest_inversion(tmp_syri(lines), "chr8")
        assert result["qry_start"] == 5000
        assert result["qry_end"] == 20000
        assert result["size"] == 15000

    def test_skips_non_inv_annotations(self, tmp_syri):
        lines = [
            _syri_line("chr8", 0, 100000, "chr8", 0, 100000, "SYN"),
            _syri_line("chr8", 0, 5000, "chr8", 0, 5000, "INV"),
        ]
        result = find_largest_inversion(tmp_syri(lines), "chr8")
        assert result["size"] == 5000

    def test_skips_wrong_chromosome(self, tmp_syri):
        lines = [
            _syri_line("chr9", 0, 50000, "chr9", 0, 50000, "INV"),
            _syri_line("chr8", 0, 100, "chr8", 0, 100, "INV"),
        ]
        result = find_largest_inversion(tmp_syri(lines), "chr8")
        assert result["size"] == 100

    def test_min_size_filter(self, tmp_syri):
        lines = [
            _syri_line("chr8", 0, 500, "chr8", 0, 500, "INV"),
            _syri_line("chr8", 1000, 5000, "chr8", 1000, 5000, "INV"),
        ]
        result = find_largest_inversion(tmp_syri(lines), "chr8", min_size=1000)
        assert result["size"] == 4000

    def test_skips_comments_and_blank_lines(self, tmp_syri):
        lines = [
            "# header comment",
            "",
            _syri_line("chr8", 0, 3000, "chr8", 0, 3000, "INV"),
        ]
        result = find_largest_inversion(tmp_syri(lines), "chr8")
        assert result["size"] == 3000

    def test_sorted_coords_handles_reversed_intervals(self, tmp_syri):
        """SyRI can output start > end for inverted regions; sorted() normalises."""
        lines = [_syri_line("chr8", 5000, 1000, "chr8", 6000, 2000, "INV")]
        result = find_largest_inversion(tmp_syri(lines), "chr8")
        assert result["ref_start"] == 1000
        assert result["ref_end"] == 5000
        assert result["qry_start"] == 2000
        assert result["qry_end"] == 6000


# ---------------------------------------------------------------------------
# Tests: error paths
# ---------------------------------------------------------------------------


class TestErrorPaths:
    def test_raises_when_no_inversion_found(self, tmp_syri):
        lines = [_syri_line("chr8", 0, 1000, "chr8", 0, 1000, "SYN")]
        with pytest.raises(RuntimeError, match="No inversion found"):
            find_largest_inversion(tmp_syri(lines), "chr8")

    def test_raises_when_all_below_min_size(self, tmp_syri):
        lines = [_syri_line("chr8", 0, 100, "chr8", 0, 100, "INV")]
        with pytest.raises(RuntimeError, match="No inversion found"):
            find_largest_inversion(tmp_syri(lines), "chr8", min_size=1000)

    def test_raises_on_non_numeric_coord(self, tmp_syri):
        # Corrupt ref_start to be non-numeric
        bad_line = "chr8\tNOT_A_NUM\t5000\t-\t-\tchr8\t1000\t4000\t-\t-\tINV"
        with pytest.raises(ValueError, match="Invalid SyRI coordinates"):
            find_largest_inversion(tmp_syri([bad_line]), "chr8")
