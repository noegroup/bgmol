
from simtk import unit
from simtk.openmm import app
import pytest

from openmmsystems.base import OpenMMSystem, OpenMMToolsTestSystem
from openmmsystems.util import OpenMMSystemsException


class SomeSystem(OpenMMSystem):
    def __init__(self, keyword_argument=True):
        super(SomeSystem, self).__init__()
        self.keyword_argument = self.system_parameter("keyword_argument", keyword_argument, default=True)
        # this way of defining the system parameters is a bit bulky; but it has its advantages
        # - autocompletion of all arguments
        # - no messing with inspect, decorators etc. -> makes it more readable


def test_parameters():
    s = SomeSystem()
    s.a = s.system_parameter("a", True, False)
    s.system_parameter("b", 1.0, 0.0)
    s.system_parameter("c", "a", "a")
    s.system_parameter("d", None, True)
    s.system_parameter("e", 1.0*unit.nanometer, None)
    s.system_parameter("f", [2, 3, 4], None)
    s.system_parameter("g", {"a": 1, "b": 2}, {"a": 2})
    with pytest.raises(OpenMMSystemsException):
        s.system_parameter("h", set(), None)
    with pytest.raises(OpenMMSystemsException):
        s.system_parameter("a", 1.0, 0.0)


def test_openmmtools_testsystem():
    """Test conversion of an openmmtools.TestSystem into an OpenMMSystem instance."""
    s = OpenMMToolsTestSystem("AlanineDipeptideImplicit", constraints=app.HBonds)
    assert len(s.positions) == 22
    assert s.system.getNumParticles() == 22
    assert s.mdtraj_topology.n_atoms == 22
    assert s.name == "AlanineDipeptideImplicit"
