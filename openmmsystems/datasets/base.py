
import os
import mdtraj as md
from openmmsystems.tpl.download import download_and_extract_archive

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
        self._coordinates = None
        self._energies = None
        self._forces = None
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
    def coordinates(self):
        return self._coordinates

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
    def system(self):
        return self._system

    def as_mdtraj(self):
        md.Trajectory(xyz=self.coordinates, topology=self.system.mdtraj_trajectory)

    def __len__(self):
        return self.num_frames



