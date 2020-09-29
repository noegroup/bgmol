import os
from openmmsystems.datasets.base import DataSet
from openmmsystems.api import system_by_name
from openmmsystems.tpl.hdf5 import HDF5TrajectoryFile, load_hdf5

__all__ = ["Ala2Implicit300", "Ala2Implicit1000"]


class Ala2Implicit300(DataSet):
    """AlanineDipeptideImplicit at 300 K.
    1 ms Langevin dynamics with 1/ps friction coefficient and 2fs time step,
    output spaced in 10 ps intervals.
    """
    url = "ftp://ftp.mi.fu-berlin.de/pub/cmb-data/openmmsystems/Ala2Implicit300.tgz"
    md5 = "ce9d6f6aa214f3eb773d52255aeaeacb"
    num_frames = 99999
    size = 49692000 # in bytes
    selection = "all"
    openmm_version = "7.4.1"
    date = "2020/09/18"

    def __init__(self, root=os.getcwd(), download: bool = False, read: bool = False):
        super(Ala2Implicit300, self).__init__(root=root, download=download, read=read)
        self._system = system_by_name("AlanineDipeptideImplicit")
        self._temperature = 300

    @property
    def trajectory_file(self):
        return os.path.join(self.root, "Ala2Implicit300/traj0.h5")

    def read(self, n_frames=None, stride=None, atom_indices=None):
        self.trajectory = load_hdf5(self.trajectory_file)
        f = HDF5TrajectoryFile(self.trajectory_file)
        frames = f.read(n_frames=n_frames, stride=stride, atom_indices=atom_indices)
        self._energies = frames.potentialEnergy
        self._forces = frames.forces
        f.close()


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
        self._system = system_by_name("AlanineDipeptideImplicit")
        self._temperature = 1000

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
