"""
Unit and regression test for the bgmol package.
"""

# Import package, test suite, and other packages as needed
import bgmol
import sys
import pytest
from click.testing import CliRunner
from bgmol.cli import main


def test_bgmol_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "bgmol" in sys.modules


@pytest.mark.parametrize("subcommand", ("systems", ))
def test_list_subcommands(subcommand):
    """Test subcommands with no options."""
    runner = CliRunner()
    result = runner.invoke(main, [subcommand])
    assert result.exit_code == 0
    assert len(result.stdout.split('\n')) > 2
