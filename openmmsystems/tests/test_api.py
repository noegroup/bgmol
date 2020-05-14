"""Test API functions."""

from openmmsystems.api import get_openmmtools_system_names


def test_get_openmmtools_system_names():
    """Check the number of testsystems in the openmmtools package."""
    assert len(get_openmmtools_system_names()) == 68
    # just to be made aware when testsystems are added or removed
