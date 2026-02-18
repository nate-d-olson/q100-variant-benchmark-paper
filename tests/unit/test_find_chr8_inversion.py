"""
Unit tests for workflow/scripts/find_chr8_inversion.py.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "workflow" / "scripts"))
from find_chr8_inversion import find_largest_inversion


def _syri_line(ref_start, ref_end, qry_start, qry_end, ann_type, chrom="chr8"):
    """Build a minimal 11-column SyRI .out line."""
    # SyRI columns (0-indexed): 0=refChr, 1=refStart, 2=refEnd, 5=qryChr,
    # 6=qryStart, 7=qryEnd, 10=type
    return "\t".join([chrom, str(ref_start), str(ref_end), "-", "-",
                      chrom, str(qry_start), str(qry_end), "-", "-", ann_type])


@pytest.fixture
def tmp_syri(tmp_path):
    def _write(lines):
        p = tmp_path / "test.syri.out"
        p.write_text("\n".join(lines) + "\n")
        return str(p)
    return _write


def test_column_index_fix(tmp_syri):
    """Regression: query coords must come from cols 6-7, not 5-6 (chr name string)."""
    lines = [_syri_line(1000, 5000, 2000, 6000, "INV")]
    result = find_largest_inversion(tmp_syri(lines), "chr8")
    assert result["qry_start"] == 2000
    assert result["qry_end"] == 6000


def test_picks_largest_inversion(tmp_syri):
    """Returns the inversion with the largest query interval."""
    lines = [
        _syri_line(0, 1000, 0, 1000, "INV"),
        _syri_line(5000, 15000, 5000, 20000, "INV"),
    ]
    result = find_largest_inversion(tmp_syri(lines), "chr8")
    assert result["size"] == 15000


def test_raises_when_no_inversion_found(tmp_syri):
    lines = [_syri_line(0, 100000, 0, 100000, "SYN")]
    with pytest.raises(RuntimeError, match="No inversion found"):
        find_largest_inversion(tmp_syri(lines), "chr8")
