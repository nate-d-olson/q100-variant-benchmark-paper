"""
Unit tests for workflow/rules/common.smk helper functions.
"""
import pytest
import sys
from pathlib import Path

# Add workflow directory to path
workflow_dir = Path(__file__).parent.parent.parent / "workflow" / "rules"
sys.path.insert(0, str(workflow_dir))


class TestExclusionHelpers:
    """Test exclusion-related helper functions."""

    def test_get_exclusion_file_path(self):
        """Test exclusion file path generation."""
        from common import get_exclusion_file_path

        result = get_exclusion_file_path("v5.0q_CHM13v2.0_smvar", "satellites", 0)
        expected = "resources/exclusions/v5.0q_CHM13v2.0_smvar/satellites_0.bed"

        assert result == expected

    def test_format_exclusion_name(self):
        """Test exclusion name formatting."""
        from common import _format_exclusion_name

        assert _format_exclusion_name("tandem-repeats") == "TANDEM_REPEATS"
        assert _format_exclusion_name("HG002Q100-errors") == "HG002Q100_ERRORS"


class TestConfigHelpers:
    """Test configuration helper functions."""

    def test_get_chromosomes_grch37(self):
        """Test chromosome list generation for GRCh37."""
        # This would require mocking the config and wildcards
        # Placeholder for demonstration
        pass

    def test_get_chromosomes_grch38(self):
        """Test chromosome list generation for GRCh38."""
        # This would require mocking the config and wildcards
        # Placeholder for demonstration
        pass


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
