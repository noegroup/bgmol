"""Base classes for systems."""
import datetime
import io
from typing import Union, Optional, List, Sequence

import packaging.version

from ..util.importing import import_openmm

mm, unit, app = import_openmm()

import numpy as np

from bgmol.util import yaml_dump, BGMolException
from bgmol.tpl import _openmmtools_testsystems
from bgmol.tpl.hdf5 import HDF5Reporter

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
            raise BGMolException(f"Parameter {name} already in use.")
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
            # allow mm.app.HBonds, ...
            return
        if not type(value) in [type(None), bool, str, float, int, list, dict, unit.Quantity, tuple]:
            raise BGMolException(
                f"Parameter type {type(value)} is not allowed for parameter: was {type(value)} ({value})"
            )
        if type(value) in [list, tuple]:
            for i, item in enumerate(value):
                OpenMMSystem._validate_parameter_type(item)
        if type(value) is dict:
            for k, v in value.items():
                OpenMMSystem._validate_parameter_type(v)
        if type(value) is unit.Quantity and type(value._value) not in [float, int]:
            raise BGMolException(
                f"Quantity value has to be of type int or float: was {type(value._value)} ({value._value})"
            )

    def reinitialize_energy_model(self, temperature=300, **kwargs):
        pass


class OpenMMSystem(BaseSystem):
    """
    The implementation is based on the openmmtoolsTestSystem class.
    It adds storing parameters in order to construct the testsystem from a compact yaml file.

    Attributes
    ----------
    system : mm.System
        System object for the test system
    positions : np.ndarray
        positions of test system
    topology : app.Topology
        topology of the test system

    """

    def __init__(self):
        """Abstract base class for test system.
        """
        super(OpenMMSystem, self).__init__()
        # Create an empty system object.
        self._system = mm.System()

        # Store positions.
        self._positions = unit.Quantity(np.zeros([0, 3], float), unit.nanometers)

        # Empty topology.
        self._topology = app.Topology()
        # MDTraj Topology is built on demand.
        self._mdtraj_topology = None
        # BGTorch Energy bridge is built on demand.
        self._energy_model = None

    @property
    def system(self):
        """The openmm.System object corresponding to the test system."""
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
        The openmm.unit.Quantity object containing the particle positions,
        with units compatible with openmm.unit.nanometers."""
        return self._positions

    @property
    def dim(self):
        return len(self._positions) * 3

    @positions.setter
    def positions(self, value):
        self._positions = value

    @positions.deleter
    def positions(self):
        del self._positions

    @property
    def topology(self):
        """The openmm.app.Topology object corresponding to the test system."""
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

    def select(self, selection):
        """Return atom indices based on selection string."""
        return self.mdtraj_topology.select(selection)

    def serialize(self):
        """Return the System and positions in serialized XML form.

        Returns
        -------

        system_xml : str
            Serialized XML form of System object.

        state_xml : str
            Serialized XML form of State object containing particle positions.

        """
        # Serialize System.
        system_xml = mm.XmlSerializer.serialize(self._system)

        # Serialize positions via State.
        if self._system.getNumParticles() == 0:
            # Cannot serialize the State of a system with no particles.
            state_xml = None
        else:
            platform = mm.Platform.getPlatformByName('Reference')
            integrator = mm.VerletIntegrator(1.0 * unit.femtoseconds)
            context = mm.Context(self._system, integrator, platform)
            context.setPositions(self._positions)
            state = context.getState(getPositions=True)
            del context, integrator
            state_xml = mm.XmlSerializer.serialize(state)
        return system_xml, state_xml

    def reinitialize_energy_model(self, temperature=300, **kwargs):
        """reinitialize the energy bridge"""
        if "integrator" in kwargs:
            integrator = kwargs["integrator"]
        else:
            integrator = mm.LangevinIntegrator(temperature, 1., 0.002)
        from bgflow.distribution.energy.openmm import OpenMMBridge, OpenMMEnergy
        energy_bridge = OpenMMBridge(self.system, integrator, **kwargs)
        self._energy_model = OpenMMEnergy(self.dim, energy_bridge)

    @property
    def energy_model(self):
        if self._energy_model is None:
            self.reinitialize_energy_model()
        return self._energy_model

    def energy(self, xyz):
        return self.energy_model.energy(xyz)

    def create_openmm_state(
            self,
            temperature: Optional[Union[unit.Quantity, float]] = None,
            velocities: Optional[np.ndarray] = None,
            positions: Optional[np.ndarray] = None,
            unitcell_vectors: Optional[Sequence[np.ndarray]] = None,
            time: Optional[float] = None,
            **parameters
    ) -> mm.State:
        """Create an openmm State object that contains positions, velocities and system parameters.

        Notes
        -----
        - By default, this method uses self.positions as the positions
        - `temperature` and `velocities` are mutually exclusive
        """
        if (velocities is not None) and (temperature is not None):
            raise ValueError("Keyword arguments 'velocities' and 'temperature' are mutually exclusive.")
        context = mm.Context(self.system, mm.VerletIntegrator(1.0 * unit.femtoseconds))
        positions = self.positions if positions is None else positions
        context.setPositions(positions)
        if velocities is not None:
            context.setVelocities(velocities)
        elif temperature is not None:
            context.setVelocitiesToTemperature(temperature)
        if unitcell_vectors is not None:
            context.setPeriodicBoxVectors(*unitcell_vectors)
        if time is not None:
            context.setTime(time)
        for name, value in parameters.items():
            context.setParameter(name, value)
        return context.getState(getPositions=True, getVelocities=True, getParameters=True, getIntegratorParameters=True)

    def create_openmm_reporters(
            self,
            out_stub: Optional[str] = None,
            interval: Optional[int] = 1000,
            **kwargs
    ) -> List[object]:
        """Creates an HDF5Reporter and (for openmm >= 7.6) an XML-CheckpointReporter

        Parameters
        ----------
        out_stub : str
            The output filename without suffix.
        interval : int
            Report interval
        kwargs : dict
            Any keyword arguments to the `HDF5Reporter`

        Returns
        -------
        reporters : List of reporters
        """
        if out_stub is None:
            now = datetime.datetime.now().strftime("%y-%m-%d-%H:%M")
            out_stub = f"{self.name}-{now}"
        reporters = []
        h5file = f"{out_stub}.h5"
        reporters.append(HDF5Reporter(file=h5file, reportInterval=interval, **kwargs))

        if packaging.version.parse(mm.__version__) >= packaging.version.parse("7.6"):
            xmlfile = f"{out_stub}.xml"
            checkpoint_reporter = app.CheckpointReporter(xmlfile, reportInterval=interval, writeState=True)
            reporters.append(checkpoint_reporter)

        return reporters

    def create_openmm_simulation(
            self,
            system: Optional[mm.System] = None,
            platform: Optional[mm.Platform] = None,
            integrator: Optional[mm.Integrator] = None,
            temperature: Optional[Union[unit.Quantity, float]] = None,
            state: Optional[mm.State] = None,
    ) -> app.Simulation:
        """Create an OpenMM simulation based on this system.
        All arguments have reasonable default choices. Either the temperature or an integrator have to be specified.

        Parameters
        ----------
        system
            By default, use self.system.
        platform
            By default, use CUDA in mixed precision on first GPU. If CUDA is not available, fall back to CPU
        integrator
            By default use a `LangevinMiddleIntegrator` with 1 fs time step and 1/ps friction coefficient.
        temperature
            The integration temperature for the default integrator.
        state
            The initial state of the simulation. By default, use self.positions and random velocities.

        Examples
        --------
        >>> system = AlanineDipeptideImplicit()
        >>> simulation = system.create_openmm_simulation(temperature=305.)
        """
        # defaults
        if (integrator is None) == (temperature is None):
            raise ValueError("Requires exactly one of the keyword arguments 'integrator' and 'temperature'.")
        system = self.system if system is None else system
        if platform is None:
            try:
                platform = mm.Platform.getPlatformByName("CUDA")
                platform.setPropertyDefaultValue("DeviceIndex", '0')
                platform.setPropertyDefaultValue("Precision", 'mixed')
            except mm.OpenMMException:
                platform = mm.Platform.getPlatformByName("CPU")
        if integrator is None:
            integrator = (
                mm.LangevinMiddleIntegrator(temperature, 1.0 / unit.picoseconds, 1.0 * unit.femtoseconds)
                if integrator is None else integrator
            )
            if _is_single_precision(platform):
                integrator.setConstraintTolerance(1e-7)
        else:
            temperature = integrator.getTemperature() if hasattr(integrator, "getTemperature") else None
        if state is None:
            state = self.create_openmm_state(temperature=temperature)

        simulation = app.Simulation(
            topology=self.topology,
            system=system,
            integrator=integrator,
            platform=platform,
        )
        simulation.context.setState(state)
        return simulation


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
        This package has a local copy of openmmtools.testsystems (in bgmol._openmmtools_testsystems)
        in order to avoid inconsistencies between systems and data.
        This copy is pinned to mmtools version 0.19.0.
        """
        super(OpenMMToolsTestSystem, self).__init__()
        TestsystemClass = getattr(_openmmtools_testsystems, name)
        assert issubclass(TestsystemClass, _openmmtools_testsystems.TestSystem)
        testsystem = TestsystemClass(**kwargs)
        self._testsystem = testsystem
        self._topology = testsystem.topology
        self._system = testsystem.system
        self._positions = testsystem.positions.value_in_unit(unit.nanometer)
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


def _is_single_precision(platform: mm.Platform) -> bool:
    return (
            "Precision" in platform.getPropertyNames()
            and platform.getPropertyDefaultValue("Precision") == 'single'
    )
