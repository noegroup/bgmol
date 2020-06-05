"""PyTest Fixtures"""

import pytest
import os
from openmmsystems import systems, api


@pytest.fixture(scope="session", params=api.get_openmmsystems_names())
def system_class(request):
    """All system classes defined in openmmsystems.systems"""
    return getattr(systems, request.param)


@pytest.fixture(scope="session")
def system_instance(system_class):
    """All instances of system classes defined in openmmsystems.systems with default parameters."""
    return system_class()


@pytest.fixture(scope='session')
def get_fn():
    test_dir = os.path.dirname(os.path.abspath(__file__))
    def _get_fn(fn):
        return '{}/data/{}'.format(test_dir, fn)
    return _get_fn