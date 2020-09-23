import os
from openmmsystems.datasets.base import DataSet
from openmmsystems.api import system_by_name
from openmmsystems.tpl.hdf5 import HDF5TrajectoryFile, load_hdf5

__all__ = ["Ala2Implicit300"]


class Ala2Implicit300(DataSet):
    """AlanineDipeptideImplicit at 300 K.
    1 ms Langevin dynamics with 1/ps friction coefficient and 2fs time step,
    output spaced in 10 ps intervals.
    """
    url = "ftp://ftp.mi.fu-berlin.de/pub/cmb-data/openmmsystems/Ala2Implicit300.tgz"
    md5 = "b85d1a2fdd76e2c9b47404855cf27c79"
    num_frames = 99999
    size = "46MB"
    selection = "all"
    openmm_version = "7.4.1"
    date = "2020/09/18"

    def __init__(self, root=os.getcwd(), download: bool = False, read: bool = False):
        super(Ala2Implicit300, self).__init__(root=root, download=download, read=read)
        self._system = system_by_name("AlanineDipeptideImplicit")
        self._temperature = 300

    @property
    def trajectory_file(self):
        return os.path.join(self.root, "300K/traj0.h5")

    def read(self, n_frames=None, stride=None, atom_indices=None):
        self.trajectory = load_hdf5(self.trajectory_file)
        f = HDF5TrajectoryFile(self.trajectory_file)
        frames = f.read(n_frames=n_frames, stride=stride, atom_indices=atom_indices)
        self._energies = frames.potentialEnergy
        self._forces = frames.forces
        f.close()
