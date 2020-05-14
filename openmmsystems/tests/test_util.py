
from openmmsystems.util import yaml_load_with_quantity, yaml_dump_with_quantity

import io

import numpy as np

from simtk.openmm import unit


def test_quantity_to_yaml():
    """Tests for saving and parsing quantity objects to and from yaml files."""
    stream = io.StringIO()
    q = [
        1.0*unit.nanometer,
        2.0*unit.dalton/unit.kilojoule_per_mole,
        3.0000000001/unit.elementary_charge,
        2.74*1e-22*unit.radians/unit.elementary_charge/unit.dalton/unit.kilojoule_per_mole,
        np.pi*unit.dimensionless
    ]
    yaml_dump_with_quantity(q, stream)
    q2 = yaml_load_with_quantity(stream.getvalue())
    for a,b in zip(q, q2):
        assert np.isclose(
            a.value_in_unit_system(unit.md_unit_system),
            b.value_in_unit_system(unit.md_unit_system),
            atol=1e-10,
            rtol=1e-10
        )
        # make sure units are consistent
        a.value_in_unit(b.unit)