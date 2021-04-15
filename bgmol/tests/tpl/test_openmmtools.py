"""Tests for openmmtools testsystems"""


from bgmol.api import list_openmmtools_systems
from bgmol.systems.base import OpenMMToolsTestSystem
import pytest


@pytest.mark.parametrize("name", list_openmmtools_systems())
def test_openmmtools_system(name):
    """Test if all openmmtools systems work with default parameters."""
    testsystem = OpenMMToolsTestSystem(name)
    assert testsystem.system.getNumParticles() > 0
    assert testsystem.name == name
    assert testsystem.parameter_names == []
