"""Tests for openmmtools testsystems"""


from openmmsystems.api import get_openmmtools_system_names
from openmmsystems.base import OpenMMToolsTestSystem
import pytest


@pytest.mark.parametrize("name", get_openmmtools_system_names())
def test_openmmtools_system(name):
    """Test if all openmmtools systems work with default parameters."""
    testsystem = OpenMMToolsTestSystem(name)
    assert testsystem.system.getNumParticles() > 0
    assert testsystem.name == name
    assert testsystem.parameter_names == []
