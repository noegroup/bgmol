
import os
import numpy as np
import mdtraj as md
from openmmsystems.tpl.download import download_and_extract_archive
from simtk.openmm import LangevinIntegrator

__all__ = ["DataSet"]


class DataSet:
    """Base class for datasets

    Parameters
    ----------
    name : str
    root : str
    download : bool
    read : bool
    """
    def __init__(self, root=os.getcwd(), download: bool = False, read: bool = False):
        # download info
        assert os.path.isdir(root), NotADirectoryError(f"root ({root}) is not a directory")
        self.root = root
        if download:
            self.download()

        # read info
        self._system = None
        self._xyz = None
        self._energies = None
        self._forces = None
        self._temperature = None
        self._trajectory = None

        if read:
            self.read()

    def download(self):
        download_and_extract_archive(
            url=self.url,
            download_root=self.root,
            extract_root=self.root,
            md5=self.md5,
            remove_finished=True
        )

    def read(self, indices=None):
        files = [os.path.join(self.root, f) for f in self.datafiles]
        for category in self._cfg["datafiles"]:
            if isinstance(self._cfg["datafiles"][category]) is str:
                ...

    @property
    def dim(self):
        return np.prod(self.xyz[0].shape)

    @property
    def xyz(self):
        return self._xyz

    @property
    def coordinates(self):
        return self._xyz

    @property
    def energies(self):
        if self._energies is None:
            raise AttributeError("This dataset contains no energies")
        return self._energies

    @property
    def forces(self):
        if self._forces is None:
            raise AttributeError("This dataset contains no forces")
        return self._forces

    @property
    def trajectory(self):
        return self._trajectory

    @trajectory.setter
    def trajectory(self, traj):
        self._trajectory = traj
        self._xyz = self._trajectory.xyz

    @property
    def system(self):
        return self._system

    @property
    def temperature(self):
        return self._temperature

    def __len__(self):
        return self.num_frames

    @property
    def energy_model(self, **kwargs):
        self.system.reinitialize_energy_model(temperature=self.temperature, **kwargs)
        return self.system.energy_model


