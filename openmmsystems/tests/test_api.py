"""Test API functions."""

import io

from openmmtools.utils import is_quantity_close

from openmmsystems.api import (
    get_openmmtools_system_names, get_openmmsystems_names, get_system_by_name, get_system_by_yaml
)

from openmmsystems.base import OpenMMToolsTestSystem
from openmmsystems.systems import ImplicitBPTI
from openmmsystems._openmmtools_testsystems import HarmonicOscillator
from simtk import unit


def test_get_openmmtools_system_names():
    """Check the number of testsystems in the openmmtools package."""
    assert len(get_openmmtools_system_names()) == 68
    # just to be made aware when testsystems are added or removed


def test_get_openmmsystems_names():
    """Check the number of testsystems in the openmmtools package."""
    assert len(get_openmmsystems_names()) > 0
    # just to be made aware when testsystems are added or removed


def test_get_system_by_name_or_yaml():
    """Check getting a system by name."""
    # openmmsystems
    ff = ["amber10.xml", "amber10_obc.xml"]
    bpti = get_system_by_name("ImplicitBPTI", forcefield=ff)
    assert isinstance(bpti, ImplicitBPTI)
    assert bpti.forcefield == ff

    # by yaml
    stream = io.StringIO(str(bpti))
    bpti2 = get_system_by_yaml(stream)
    assert isinstance(bpti2, ImplicitBPTI)
    assert bpti2.forcefield == ff

    # openmmtools testsystem
    harmonic = get_system_by_name("HarmonicOscillator", K=4*unit.kilojoule_per_mole/unit.nanometer**2)
    assert isinstance(harmonic, OpenMMToolsTestSystem)
    assert isinstance(harmonic._testsystem, HarmonicOscillator)
    assert harmonic.name == "HarmonicOscillator"
    assert is_quantity_close(harmonic.K, 4*unit.kilojoule_per_mole/unit.nanometer**2)

    # by yaml
    stream = io.StringIO(str(harmonic))
    print(stream.getvalue())
    harmonic2 = get_system_by_yaml(stream)
    assert isinstance(harmonic2, OpenMMToolsTestSystem)
    assert isinstance(harmonic._testsystem, HarmonicOscillator)
    assert harmonic2.name == "HarmonicOscillator"
    assert is_quantity_close(harmonic2.K, 4*unit.kilojoule_per_mole/unit.nanometer**2)