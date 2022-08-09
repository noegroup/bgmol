"""Datasets of water systems"""

__all__ = ["WaterDimerFlexibleTIP3P", "WaterDimerRigidTIP3P"]


import os

from .base import DataSet
from ..systems.base import OpenMMToolsTestSystem
from ..util.importing import import_openmm
_, unit, _ = import_openmm()


class WaterDimerFlexibleTIP3P(DataSet):
    """Two flexible water molecules in a harmonic external field.
    Created by Langevin dynamics using a 0.1 fs time step and a 1/ps friction coefficient.
    Frames spaced 1 ps apart.
    The dataset contains positions, forces, energies.
    """
    url = "http://ftp.mi.fu-berlin.de/pub/cmb-data/bgmol/datasets/water/WaterDimerFlexibleTIP3P.tgz"
    md5 = "865fa3abc96bbc9765c56a0d8c94790e"
    num_frames = 50000
    size = 7263250  # in bytes
    selection = "all"
    openmm_version = "7.7"
    date = "2022/08/09"
    author = "Andreas Krämer"

    def __init__(self, root=os.getcwd(), download=True,  read=True):
        self.trajectory_file = os.path.join(root, "WaterDimerFlexibleTIP3P/water_unconstrained_0.h5")
        super().__init__(root=root, download=download, read=read)
        force_constant = 3 * unit.kilojoule_per_mole / unit.nanometer ** 2
        self._system = OpenMMToolsTestSystem(
            'WaterCluster',
            n_waters=2,
            K=force_constant,
            constrained=False
        )

    def read(self, n_frames=None, stride=None, atom_indices=None):
        self.load_hdf5()


class WaterDimerRigidTIP3P(DataSet):
    """Two rigid water molecules in a harmonic external field.
    Created by Langevin dynamics using a 1 fs time step and a 1/ps friction coefficient.
    Frames spaced 1 ps apart.
    The dataset contains positions, forces, energies.
    """
    url = "http://ftp.mi.fu-berlin.de/pub/cmb-data/bgmol/datasets/water/WaterDimerRigidTIP3P.tgz"
    md5 = "611ab61c3a9e850927cf4b520af26f45"
    num_frames = 50000
    size = 7076924  # in bytes
    selection = "all"
    openmm_version = "7.7"
    date = "2022/08/09"
    author = "Andreas Krämer"

    def __init__(self, root=os.getcwd(), download=True,  read=True):
        self.trajectory_file = os.path.join(root, "WaterDimerRigidTIP3P/water_constrained_0.h5")
        super().__init__(root=root, download=download, read=read)
        force_constant = 3 * unit.kilojoule_per_mole / unit.nanometer ** 2
        self._system = OpenMMToolsTestSystem(
            'WaterCluster',
            n_waters=2,
            K=force_constant,
            constrained=True
        )

    def read(self, n_frames=None, stride=None, atom_indices=None):
        self.load_hdf5()

