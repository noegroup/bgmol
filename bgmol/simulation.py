

from bgmol.util.importing import import_openmm
mm, unit, app = import_openmm()
from .tpl.hdf5 import HDF5Reporter


class Simulator:
    def __init__(self, system):
        self.integrator = mm.LangevinIntegrator(300*unit.kelvin, 1./unit.picosecond, 2.*unit.femtosecond)
        self.platform = mm.Platform.getPlatformByName("CPU")
        self.equilibration=500*unit.picosecond
        self.simulation=1*unit.microsecond
        self.reporters = [
            HDF5Reporter("traj.h5", 5000, forces=True, velocities=True)
        ]

    def equilibrate(self):
        pass

