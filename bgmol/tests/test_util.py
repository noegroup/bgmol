
from bgmol.util import parse_quantity, yaml_load, yaml_dump

import io

import numpy as np

from simtk.openmm import unit
from simtk.openmm.app import *

from openmmtools.utils import is_quantity_close


def test_parse_quantity():
    for str, res in [
        ("1 * nanometer", 1 * unit.nanometer)
    ]:
        assert is_quantity_close(parse_quantity(str), res)


def test_quantity_to_yaml():
    """Tests for saving and parsing quantity objects to and from yaml files."""
    stream = io.StringIO()
    q = [
        1*unit.nanometer,
        2.0*unit.dalton/unit.kilojoule_per_mole,
        3.0000000001/unit.elementary_charge,
        2.74*1e-22*unit.radians/unit.elementary_charge/unit.dalton/unit.kilojoule_per_mole,
        np.pi/unit.nanometer**2
    ]
    yaml_dump(q, stream)
    q2 = yaml_load(stream.getvalue())
    for a,b in zip(q, q2):
        print(a,b)
        assert is_quantity_close(a,b)


def test_singleton_to_yaml():
    """Tests for saving and parsing openmm Singleton objects to and from yaml files."""
    stream = io.StringIO()
    q = [HBonds, PME, LJPME]
    yaml_dump(q, stream)
    q2 = yaml_load(stream.getvalue())
    for a,b in zip(q, q2):
        assert a == b

