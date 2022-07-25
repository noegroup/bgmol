import os

from .base import DataSet
from ..systems.base import OpenMMToolsTestSystem
from ..util.importing import import_openmm
_, unit, _ = import_openmm()


class DimerFlexibleTIP3P(DataSet):
    """Two flexible water molecule in a harmonic external field.
    50000 steps of Langevin dynamics using a 0.5 fs time step.
    The dataset contains positions, forces, energies.
    """
    url = "http://ftp.mi.fu-berlin.de/pub/cmb-data/bgmol/datasets/water/DimerFlexibleTIP3P.tgz"
    md5 = "89d94c626d73ba857eb9c04504b2c537"
    num_frames = 50000
    size = 7048084  # in bytes
    selection = "all"
    openmm_version = "7.7"
    date = "2022/07/25"
    author = "Andreas Krämer"

    def __init__(self, root=os.getcwd(), download=True,  read=True):
        self.trajectory_file = os.path.join(root, "DimerFlexibleTIP3P/waters_flexible.h5")
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


