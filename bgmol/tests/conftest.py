"""PyTest Fixtures"""

import pytest
import os
import torch
from bgmol.datasets import Ala2Implicit1000Test
from itertools import product
from bgmol.systems import ChignolinC22Implicit
from bgmol.systems.fastfolders import FastFolder, FAST_FOLDER_NAMES
from bgmol import systems, api


@pytest.fixture(scope="session", params=api.list_bgmol())
def system_class(request):
    """All system classes defined in bgmol.systems"""
    return getattr(systems, request.param)


@pytest.fixture(scope="session")
def system_instance(system_class):
    """All instances of system classes defined in bgmol.systems with default parameters."""
    return system_class()


@pytest.fixture(scope='session')
def get_fn():
    test_dir = os.path.dirname(os.path.abspath(__file__))
    def _get_fn(fn):
        return '{}/data/{}'.format(test_dir, fn)
    return _get_fn


@pytest.fixture(scope="session")
def ala2dataset(tmpdir_factory):
    return Ala2Implicit1000Test(root=tmpdir_factory.mktemp("ala2dataset"), download=True, read=True)


@pytest.fixture(scope="session")
def chignolin(tmpdir_factory):
    return ChignolinC22Implicit(root=tmpdir_factory.mktemp("chignolin_system"))


@pytest.fixture(scope="session", params=product(
    FAST_FOLDER_NAMES,
    (False, True),
    ("charmm36m", ["amber99sbildn.xml", "tip3p.xml"])
))
def fastfolder_system(request, tmpdir_factory):
    protein, solvated, forcefield = request.param
    tmpdir = tmpdir_factory.mktemp(protein)
    if forcefield != "charmm36m" and protein in ["ntl9", "villin"]:
        pytest.skip("Nonstandard residues.")
    yield FastFolder(protein, download=True, solvated=solvated, forcefield=forcefield, root=str(tmpdir))


# skipping slow tests by default
# ==============================
def pytest_addoption(parser):
    parser.addoption(
        "--runslow", action="store_true", default=False, help="run slow tests"
    )


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: mark test as slow to run")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--runslow"):
        # --runslow given in cli: do not skip slow tests
        return
    skip_slow = pytest.mark.skip(reason="need --runslow option to run")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip_slow)


@pytest.fixture(
    params=[
        "cpu",
        pytest.param(
            "cuda:0",
            marks=pytest.mark.skipif(
                not torch.cuda.is_available(),
                reason="CUDA not available."
            )
        )
    ]
)


def device(request):
    """Run a test case for all available devices."""
    return request.param


@pytest.fixture(params=[torch.float32, torch.float64])
def dtype(request, device):
    """Run a test case in single and double precision."""
    return request.param


@pytest.fixture()
def ctx(dtype, device):
    return {"dtype": dtype, "device": device}

