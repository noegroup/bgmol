"""PyTest Fixtures"""

import pytest
from openmmsystems import systems, api


@pytest.fixture(scope="session", params=api.get_openmmsystems_names())
def system_class(request):
    """All system classes defined in openmmsystems.systems"""
    return getattr(systems, request.param)


@pytest.fixture(scope="session")
def system_instance(system_class):
    """All instances of system classes defined in openmmsystems.systems with default parameters."""
    return system_class()
