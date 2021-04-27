
import os
import numpy as np

from simtk import unit
from simtk.openmm.app import HBonds
from simtk.openmm import LangevinIntegrator, Platform

from .base import DataSet
from ..systems.minipeptides import MiniPeptide
from ..tpl.hdf5 import load_hdf5, HDF5TrajectoryFile


__all__ = ["ASolvatedAmoeba", "ASolvatedAmber14", "ASolvatedAmber99"]


class ASolvatedAmoeba(DataSet):
    """Capped alanine in explicit water with Amoeba-2013.
    20 ns Langevin dynamics with 1/ps friction coefficient and 0.5 fs time step.
    Samples are spaced in 1 ps intervals.
    The dataset contains positions, forces, and energies.
    """
    url = "ftp://ftp.mi.fu-berlin.de/pub/cmb-data/bgmol/datasets/minipeptides/ASolvatedAmoeba.tgz"
    md5 = "0f5aab1b946e313ee4139cf1d7645199"
    num_frames = 20000
    size = 1253106
    selection = "all"
    openmm_version = "7.5.0"
    date = "2021/01/15"
    author = "Andreas Krämer"

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
    url = "ftp://ftp.mi.fu-berlin.de/pub/cmb-data/bgmol/datasets/minipeptides/ASolvatedAmber14.tgz"
    md5 = "0bccb8e06dabf8a7fc11a407d7a360b3"
    num_frames = 20000
    size = 1252535
    selection = "all"
    openmm_version = "7.5.0"
    date = "2021/01/15"
    author = "Andreas Krämer"

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


class ASolvatedAmber99(DataSet):
    """Capped alanine in explicit water with Amber99 and TIP3P.
    2 microsecond samples spaced in 2 ps intervals.
    The dataset contains protein all-atom positions, protein all-atom forces, and energies as well as
    specific forces and energies from the solvent environment.
    """
    url = "ftp://ftp.mi.fu-berlin.de/pub/cmb-data/bgmol/datasets/minipeptides/ASolvatedAmber99.tgz"
    md5 = "94b4a5221014ddfb03fb62cc5c7d67df"
    num_frames = 1000000
    size = 749056  # in bytes
    selection = "protein"
    openmm_version = "7.4.1"
    date = "2020/05/30"
    author = "Yaoyi Chen"

    def __init__(self, root=os.getcwd(), download: bool = False, read: bool = False):
        super(ASolvatedAmber99, self).__init__(root=root, download=download, read=read)
        self._system = MiniPeptide(
            aminoacids="A",
            solvated=True,
            forcefield=["amber99sbildn.xml", "tip3p.xml"],
            nonbonded_cutoff=0.9 * unit.nanometer,
            constraints=HBonds,
            hydrogen_mass=4 * unit.amu
        )
        self._temperature = 300.

    def read(self, n_frames=None, stride=None, atom_indices=None):
        npz = np.load(os.path.join(self.root, "ASolvatedAmber99/spep_0000_dataset.npz"))
        self._xyz = npz["coords"]
        self._forces = npz["Fs"]
        self._energies = npz["Es"]
        self.solvation_forces = npz["solvFs"]
        self.solvation_energies = npz["solvEs"]

    @property
    def integrator(self):
        integrator = LangevinIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 0.004 * unit.picoseconds)
        integrator.setConstraintTolerance(1e-5)
        return integrator

    @property
    def platform(self):
        platform = Platform.getPlatformByName("CUDA")
        platform.setPropertyDefaultValue("Precision", "mixed")
        return platform
