
import os

from simtk import unit

from .base import DataSet
from ..systems.minipeptides import MiniPeptide
from ..tpl.hdf5 import load_hdf5, HDF5TrajectoryFile


class Ala2Implicit1000(DataSet):
    """AlanineDipeptideImplicit at 300 K.
    1 ms Langevin dynamics with 1/ps friction coefficient and 2fs time step,
    output spaced in 10 ps intervals.
    """
    url = "ftp://ftp.mi.fu-berlin.de/pub/cmb-data/openmmsystems/Ala2Implicit1000.tgz"
    md5 = "9a90965254dac7440976d0d62687caed"
    num_frames = 100000
    size = 51071715 # in bytes
    selection = "all"
    openmm_version = "7.4.1"
    date = "2020/09/29"

    def __init__(self, root=os.getcwd(), download: bool = False, read: bool = False):
        super(Ala2Implicit1000, self).__init__(root=root, download=download, read=read)
        self._system = MiniPeptide(
            aminoacids="A",
            solvated=True,
            forcefield=["amoeba2013.xml"],
            nonbonded_cutoff=1.2*unit.nanometer,
            switch_distance=1.0*unit.nanometer
        )
        self._temperature = 300.

    @property
    def trajectory_file(self):
        return os.path.join(self.root, "Ala2Implicit1000/traj0.h5")

    def read(self, n_frames=None, stride=None, atom_indices=None):
        self.trajectory = load_hdf5(self.trajectory_file)
        f = HDF5TrajectoryFile(self.trajectory_file)
        frames = f.read(n_frames=n_frames, stride=stride, atom_indices=atom_indices)
        self._energies = frames.potentialEnergy
        self._forces = frames.forces
        f.close()


class ASolvatedAmoeba(DataSet):
    """Capped alanine in explicit water with Amoeba-2013.
    20 ns Langevin dynamics with 1/ps friction coefficient and 0.5 ps time step.
    Samples are spaced in 1 ps intervals.
    The dataset contains positions, forces, and energies.
    """
    url = "ftp://ftp.mi.fu-berlin.de/pub/cmb-data/openmmsystems/minipeptides/ASolvatedAmoeba.tgz"
    md5 = "0f5aab1b946e313ee4139cf1d7645199"
    num_frames = 20000
    size = 1253106
    selection = "all"
    openmm_version = "7.5.0"
    date = "2021/01/15"

    @property
    def trajectory_file(self):
        return os.path.join(self.root, "ASolvatedAmoeba/traj0.h5")

    def __init__(self, root=os.getcwd(), download: bool = False, read: bool = False):
        super(ASolvatedAmoeba, self).__init__(root=root, download=download, read=read)
        self._system = MiniPeptide(
            aminoacids="A",
            solvated=True,
            forcefield=["amoeba2013.xml"],
            nonbonded_cutoff=1.2*unit.nanometer,
            switch_distance=1.0*unit.nanometer
        )
        self._temperature = 300.

    def read(self, n_frames=None, stride=None, atom_indices=None):
        self.trajectory = load_hdf5(self.trajectory_file)
        f = HDF5TrajectoryFile(self.trajectory_file)
        frames = f.read(n_frames=n_frames, stride=stride, atom_indices=atom_indices)
        self._energies = frames.potentialEnergy
        self._forces = frames.forces
        f.close()


class ASolvatedAmber14(DataSet):
    """Capped alanine in explicit water with Amber14 and TIP3P.
    20 ns Langevin dynamics with 1/ps friction coefficient and 1 ps time step.
    Samples are spaced in 1 ps intervals.
    The dataset contains positions, forces, and energies.
    """
    url = "ftp://ftp.mi.fu-berlin.de/pub/cmb-data/openmmsystems/minipeptides/ASolvatedAmber14.tgz"
    md5 = "0bccb8e06dabf8a7fc11a407d7a360b3"
    num_frames = 20000
    size = 1252535
    selection = "all"
    openmm_version = "7.5.0"
    date = "2021/01/15"

    @property
    def trajectory_file(self):
        return os.path.join(self.root, "ASolvatedAmber14/traj0.h5")

    def __init__(self, root=os.getcwd(), download: bool = False, read: bool = False):
        super(ASolvatedAmber14, self).__init__(root=root, download=download, read=read)
        self._system = MiniPeptide(
            aminoacids="A",
            solvated=True,
            forcefield=["amber14/protein.ff14SB.xml", "amber14/tip3p.xml"],
            nonbonded_cutoff=1.2 * unit.nanometer,
            switch_distance=1.0 * unit.nanometer
        )
        self._temperature = 300.

    def read(self, n_frames=None, stride=None, atom_indices=None):
        self.trajectory = load_hdf5(self.trajectory_file)
        f = HDF5TrajectoryFile(self.trajectory_file)
        frames = f.read(n_frames=n_frames, stride=stride, atom_indices=atom_indices)
        self._energies = frames.potentialEnergy
        self._forces = frames.forces
        f.close()
