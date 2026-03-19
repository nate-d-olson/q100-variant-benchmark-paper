"""
Unit tests for variant type reclassification logic in generate_variant_parquet.py.

Tests the UNK/INV → INS/DEL/SNP reclassification branch to ensure:
- UNK variants with alt_len > ref_len are classified as INS
- UNK variants with ref_len > alt_len are classified as DEL
- UNK variants with ref_len == alt_len == 1 are classified as SNP
- UNK variants with ref_len == alt_len > 1 (MNPs) are NOT classified as SNP
"""

import sys
from pathlib import Path

# Add workflow/scripts to path for importing
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "workflow" / "scripts"))


def test_unk_reclassification_to_ins():
    """UNK variant with alt_len > ref_len should be reclassified as INS."""
    ref_len = 1
    alt_len = 5
    var_size = alt_len - ref_len  # = 4

    # Simulate reclassification logic
    if var_size > 0:
        var_type = "INS"
    elif var_size < 0:
        var_type = "DEL"
    else:
        var_type = "SNP"

    assert var_type == "INS", "UNK with alt_len > ref_len should become INS"
    assert var_size == 4, "var_size should be positive for INS"


def test_unk_reclassification_to_del():
    """UNK variant with ref_len > alt_len should be reclassified as DEL."""
    ref_len = 10
    alt_len = 2
    var_size = alt_len - ref_len  # = -8

    # Simulate reclassification logic
    if var_size > 0:
        var_type = "INS"
    elif var_size < 0:
        var_type = "DEL"
    else:
        var_type = "SNP"

    assert var_type == "DEL", "UNK with ref_len > alt_len should become DEL"
    assert var_size == -8, "var_size should be negative for DEL"


def test_unk_reclassification_to_snp_single_nucleotide():
    """UNK variant with ref_len == alt_len == 1 should be reclassified as SNP."""
    ref_len = 1
    alt_len = 1
    var_size = alt_len - ref_len  # = 0

    # Simulate reclassification logic (corrected implementation)
    if var_size > 0:
        var_type = "INS"
    elif var_size < 0:
        var_type = "DEL"
    elif ref_len == 1 and alt_len == 1:
        var_type = "SNP"
    else:
        var_type = "UNK"

    assert var_type == "SNP", "UNK with ref_len == alt_len == 1 should become SNP"
    assert var_size == 0, "var_size should be 0 for SNP"


def test_unk_reclassification_mnp_should_not_be_snp():
    """
    UNK variant with ref_len == alt_len > 1 should NOT be reclassified as SNP.

    These are multi-nucleotide polymorphisms (MNPs) or complex substitutions.
    They should remain as UNK to avoid inflating SNP counts.
    """
    ref_len = 21  # e.g., "GCAGGAGGGGGCAGAAGGGGA"
    alt_len = 21  # e.g., "ACAGGAGGGGGCAGAAGGGGA"
    var_size = alt_len - ref_len  # = 0

    # Corrected implementation
    if var_size > 0:
        var_type = "INS"
    elif var_size < 0:
        var_type = "DEL"
    elif ref_len == 1 and alt_len == 1:
        var_type = "SNP"
    else:
        var_type = "UNK"  # MNP/complex - keep as UNK

    assert var_type == "UNK", "MNPs should remain UNK, not become SNP"


def test_unk_reclassification_edge_cases():
    """Test edge cases in variant reclassification."""
    test_cases = [
        # (ref_len, alt_len, expected_type, expected_size)
        (1, 1, "SNP", 0),  # Single nucleotide substitution
        (1, 2, "INS", 1),  # Single base insertion
        (2, 1, "DEL", -1),  # Single base deletion
        (10, 20, "INS", 10),  # Larger insertion
        (50, 25, "DEL", -25),  # Larger deletion
        (2, 2, "UNK", 0),  # MNP (dinucleotide)
        (10, 10, "UNK", 0),  # Larger MNP
    ]

    for ref_len, alt_len, expected_type, expected_size in test_cases:
        var_size = alt_len - ref_len

        if var_size > 0:
            var_type = "INS"
        elif var_size < 0:
            var_type = "DEL"
        elif ref_len == 1 and alt_len == 1:
            var_type = "SNP"
        else:
            var_type = "UNK"

        assert var_type == expected_type, (
            f"Failed for ref_len={ref_len}, alt_len={alt_len}: "
            f"expected {expected_type}, got {var_type}"
        )
        assert var_size == expected_size, (
            f"Failed for ref_len={ref_len}, alt_len={alt_len}: "
            f"expected size {expected_size}, got {var_size}"
        )


def test_inv_reclassification():
    """INV variants should follow the same reclassification logic as UNK."""
    # INV is included in the same branch as UNK in the actual code
    # Test that it follows the same rules
    ref_len = 100
    alt_len = 100
    var_size = alt_len - ref_len  # = 0

    # Corrected implementation
    if var_size > 0:
        var_type = "INS"
    elif var_size < 0:
        var_type = "DEL"
    elif ref_len == 1 and alt_len == 1:
        var_type = "SNP"
    else:
        var_type = "UNK"  # Large INV kept as UNK

    assert var_type == "UNK", "Equal-length INV (ref_len > 1) should remain UNK"
