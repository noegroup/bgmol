"""Base classes for systems."""


import io

from simtk import openmm, unit
from simtk.openmm import app
import numpy as np

from openmmsystems.util import yaml_dump, OpenMMSystemsException
from openmmsystems.tpl import _openmmtools_testsystems

__all__ = ["BaseSystem", "OpenMMSystem", "OpenMMToolsTestSystem"]


class BaseSystem:
    """Abstract base class for testsystems.

    Attributes
    ----------
    name : str
        The identifier of the system (usually the class name).
    parameter_names : list of str
        Names of parameters that have been registered for this system.
    """
    def __init__(self):
        self._parameter_defaults = {}

    @property
    def name(self):
        """The name of the test system."""
        return self.__class__.__name__

    def system_parameter(self, name, value, default):
        """
        Register a system parameter.
        """
        if name in self._parameter_defaults:
            raise OpenMMSystemsException(f"Parameter {name} already in use.")
        self._validate_parameter_type(value)
        self._validate_parameter_type(default)
        self._parameter_defaults[name] = default
        setattr(self, name, value)
        return value

    @property
    def parameter_names(self):
        return list(self._parameter_defaults.keys())

    def __str__(self):
        stream = io.StringIO()
        parameters = {name: getattr(self, name) for name in self._parameter_defaults}
        yaml_dump(
            {"system": {"identifier": self.name, "parameters": parameters}},
            stream
        )
        return stream.getvalue()

    @staticmethod
    def _validate_parameter_type(value):
        """Allow only some types for parameters."""
        if isinstance(value, app.internal.singleton.Singleton):
            # allow openmm.app.HBonds, ...
            return
        if not type(value) in [type(None), bool, str, float, int, list, dict, unit.Quantity, tuple]:
            raise OpenMMSystemsException(
                f"Parameter type {type(value)} is not allowed for parameter: was {type(value)} ({value})"
            )
        if type(value) in [list, tuple]:
            for i, item in enumerate(value):
                OpenMMSystem._validate_parameter_type(item)
        if type(value) is dict:
            for k,v in value.items():
                OpenMMSystem._validate_parameter_type(v)
        if type(value) is unit.Quantity and type(value._value) not in [float, int]:
            raise OpenMMSystemsException(
                f"Quantity value has to be of type int or float: was {type(value._value)} ({value._value})"
            )

    def reinitialize_energy_model(self, temperature=300, **kwargs):
        pass


class OpenMMSystem(BaseSystem):
    """
    The implementation is based on the openmmtools.TestSystem class.
    It adds storing parameters in order to construct the testsystem from a compact yaml file.

    Attributes
    ----------
    system : simtk.openmm.System
        System object for the test system
    positions : list
        positions of test system
    topology : list
        topology of the test system

    """

    def __init__(self):
        """Abstract base class for test system.
        """
        super(OpenMMSystem, self).__init__()
        # Create an empty system object.
        self._system = openmm.System()

        # Store positions.
        self._positions = unit.Quantity(np.zeros([0, 3], np.float), unit.nanometers)

        # Empty topology.
        self._topology = app.Topology()
        # MDTraj Topology is built on demand.
        self._mdtraj_topology = None
        # BGTorch Energy bridge is built on demand.
        self._energy_model = None

    @property
    def system(self):
        """The simtk.openmm.System object corresponding to the test system."""
        return self._system

    @system.setter
    def system(self, value):
        self._system = value

    @system.deleter
    def system(self):
        del self._system

    @property
    def positions(self):
        """particle positions
        The simtk.unit.Quantity object containing the particle positions,
        with units compatible with simtk.unit.nanometers."""
        return self._positions

    @property
    def dim(self):
        return len(self._positions)*3

    @positions.setter
    def positions(self, value):
        self._positions = value

    @positions.deleter
    def positions(self):
        del self._positions

    @property
    def topology(self):
        """The simtk.openmm.app.Topology object corresponding to the test system."""
        return self._topology

    @topology.setter
    def topology(self, value):
        self._topology = value
        self._mdtraj_topology = None

    @topology.deleter
    def topology(self):
        del self._topology

    @property
    def mdtraj_topology(self):
        """The mdtraj.Topology object corresponding to the test system (read-only)."""
        import mdtraj as md
        if self._mdtraj_topology is None:
            self._mdtraj_topology = md.Topology.from_openmm(self._topology)
        return self._mdtraj_topology

    def serialize(self):
        """Return the System and positions in serialized XML form.

        Returns
        -------

        system_xml : str
            Serialized XML form of System object.

        state_xml : str
            Serialized XML form of State object containing particle positions.

        """

        from simtk.openmm import XmlSerializer

        # Serialize System.
        system_xml = XmlSerializer.serialize(self._system)

        # Serialize positions via State.
        if self._system.getNumParticles() == 0:
            # Cannot serialize the State of a system with no particles.
            state_xml = None
        else:
            platform = openmm.Platform.getPlatformByName('Reference')
            integrator = openmm.VerletIntegrator(1.0 * unit.femtoseconds)
            context = openmm.Context(self._system, integrator, platform)
            context.setPositions(self._positions)
            state = context.getState(getPositions=True)
            del context, integrator
            state_xml = XmlSerializer.serialize(state)
        return system_xml, state_xml

    def reinitialize_energy_model(self, temperature=300, **kwargs):
        """reinitialize the energy bridge"""
        if "integrator" in kwargs:
            integrator = kwargs["integrator"]
        else:
            integrator = openmm.LangevinIntegrator(temperature, 1., 0.002)
        from bgtorch.distribution.energy.openmm import OpenMMBridge, OpenMMEnergy
        energy_bridge = OpenMMBridge(self.system, integrator, **kwargs)
        self._energy_model = OpenMMEnergy(self.dim, energy_bridge)

    @property
    def energy_model(self):
        if self._energy_model is None:
            self.reinitialize_energy_model()
        return self._energy_model

    def energy(self, xyz):
        return self.energy_model.energy(xyz)


class OpenMMToolsTestSystem(OpenMMSystem):
    """An openmmtools.TestSystem in disguise."""
    def __init__(self, name, **kwargs):
        """
        Parameters
        ----------
        name (str) :
            The class name of an openmmtools.TestSystem subclass.
        **kwargs : (optional)
            Keyword arguments that are passed to the constructor of the testsystem

        Notes
        -----
        This package has a local copy of openmmtools.testsystems (in openmmsystems._openmmtools_testsystems)
        in order to avoid inconsistencies between systems and data.
        This copy is pinned to openmmtools version 0.19.0.
        """
        super(OpenMMToolsTestSystem, self).__init__()
        TestsystemClass = getattr(_openmmtools_testsystems, name)
        assert issubclass(TestsystemClass, _openmmtools_testsystems.TestSystem)
        testsystem = TestsystemClass(**kwargs)
        self._testsystem = testsystem
        self._topology = testsystem.topology
        self._system = testsystem.system
        self._positions = testsystem.positions
        self._name = name
        # We register only the non-default arguments.
        # It's hard to track down the actual defaults even using the inspect module.
        for key, value in kwargs.items():
            self.system_parameter(key, value, default=None)

    @property
    def name(self):
        """The name of the test system."""
        return self._name

    def __getattr__(self, item):
        return getattr(self._testsystem, item)


class TorchSystem(OpenMMSystem):
    @property
    def energy_model(self):
        return self