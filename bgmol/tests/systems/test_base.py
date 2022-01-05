"""Test functionality of the system classes in the base module"""

import os
import packaging

from bgmol.util.importing import import_openmm
mm, unit, app = import_openmm()
import pytest

from bgmol.systems.base import OpenMMSystem, OpenMMToolsTestSystem
from bgmol.systems.ala2 import AlanineDipeptideImplicit
from bgmol.util import BGMolException
from bgmol.tpl.hdf5 import HDF5TrajectoryFile


class SomeSystem(OpenMMSystem):
    def __init__(self, keyword_argument=True):
        super(SomeSystem, self).__init__()
        self.keyword_argument = self.system_parameter("keyword_argument", keyword_argument, default=True)
        # this way of defining the system parameters is a bit bulky; but it has its advantages
        # - autocompletion of all arguments
        # - no messing with inspect, decorators etc. -> makes it more readable


def test_parameters():
    """Test specification of system parameters."""
    s = SomeSystem()
    s.a = s.system_parameter("a", True, False)
    s.system_parameter("b", 1.0, 0.0)
    s.system_parameter("c", "a", "a")
    s.system_parameter("d", None, True)
    s.system_parameter("e", 1.0*unit.nanometer, None)
    s.system_parameter("f", [2, 3, 4], None)
    s.system_parameter("g", {"a": 1, "b": 2}, {"a": 2})
    with pytest.raises(BGMolException):
        s.system_parameter("h", set(), None)
    with pytest.raises(BGMolException):
        s.system_parameter("a", 1.0, 0.0)


def test_create_simulation():
    system = AlanineDipeptideImplicit()
    with pytest.raises(ValueError):
        system.create_openmm_simulation()
    with pytest.raises(ValueError):
         system.create_openmm_simulation(temperature=200., integrator=...)
    simulation = system.create_openmm_simulation(temperature=305.)
    integrator = simulation.context.getIntegrator()
    assert integrator.getTemperature().value_in_unit(unit.kelvin) == pytest.approx(305.)
    simulation.step(10)


def test_create_reporters(tmpdir):
    system = AlanineDipeptideImplicit()
    simulation = system.create_openmm_simulation(temperature=305.)
    simulation.reporters = system.create_openmm_reporters(tmpdir/"sim", interval=10)
    simulation.step(40)
    del simulation  # to close files
    hdf5 = HDF5TrajectoryFile(tmpdir/"sim.h5")
    traj = hdf5.read_as_traj()
    assert traj.n_frames == 4

    if packaging.version.parse(mm.__version__) >= packaging.version.parse("7.7"):
        assert os.path.isfile(tmpdir/"sim.xml")


def test_openmmtools_testsystem():
    """Test conversion of an openmmtools.TestSystem into an OpenMMSystem instance."""
    s = OpenMMToolsTestSystem("AlanineDipeptideImplicit", constraints=app.HBonds)
    assert len(s.positions) == 22
    assert s.system.getNumParticles() == 22
    assert s.mdtraj_topology.n_atoms == 22
    assert s.name == "AlanineDipeptideImplicit"
