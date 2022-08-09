
import os
from typing import Sequence, Callable
import numpy as np
import mdtraj as md
import torch
import torch.utils.data

from bgmol.tpl.hdf5 import load_hdf5, HDF5TrajectoryFile

from ..util.importing import import_openmm
mm, unit, app = import_openmm()

from ..tpl.download import download_and_extract_archive

__all__ = ["DataSet"]


class DataSet:
    """Base class for datasets

    Parameters
    ----------
    root : str
        The root directory where the dataset is stored locally.
    download : bool
        Whether to download the dataset.
    read : bool
        Whether to read the dataset into memory.

    Attributes
    ----------
    xyz : np.ndarray
        The coordinates with shape (num_frames, num_atoms_in_selection, 3)
    coordinates : np.ndarray
        An alias for DataSet.xyz
    forces: np.ndarray or None
         The forces with shape (num_frames, num_atoms_in_selection, 3)
    energies: np.ndarray or None
         The potential energies with shape (num_frames)
    trajectory: mdtraj.Trajectory
        The trajectory.
    temperature: float
        Temperature in Kelvin.
    unitcell_vectors: np.ndarray or None
        The box vectors of the simulation cell with shape (num_frames, 3, 3).
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
        self._unitcell_vectors = None

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

    def read(self, n_frames=None, stride=None, atom_indices=None):
        files = [os.path.join(self.root, f) for f in self.datafiles]
        for category in self._cfg["datafiles"]:
            if isinstance(self._cfg["datafiles"][category]) is str:
                ...  # TODO

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
        return self._energies
        #if self._energies is None:
        #    raise AttributeError("This dataset contains no energies")
        #return self._energies

    @property
    def forces(self):
        return self._forces
        #if self._forces is None:
        #    raise AttributeError("This dataset contains no forces")
        #return self._forces

    @property
    def trajectory(self):
        if self._trajectory is None:
            trajectory = md.Trajectory(
                xyz=self.xyz,
                topology=self.system.mdtraj_topology
            )
            if self.unitcell_vectors is not None:
                trajectory.unitcell_vectors = self.unitcell_vectors
            # TODO: set time when possible
            self._trajectory = trajectory
        return self._trajectory

    @trajectory.setter
    def trajectory(self, traj):
        self._trajectory = traj
        self._xyz = self._trajectory.xyz
        self._unitcell_vectors = self.trajectory.unitcell_vectors

    @property
    def unitcell_vectors(self):
        return self._unitcell_vectors

    @property
    def system(self):
        return self._system

    @property
    def temperature(self):
        return self._temperature

    def __len__(self):
        return self.num_frames

    def get_energy_model(self, **kwargs):
        self.system.reinitialize_energy_model(temperature=self.temperature, **kwargs)
        return self.system.energy_model

    def load_hdf5(self, n_frames=None, stride=None, atom_indices=None):
        self.trajectory = load_hdf5(self.trajectory_file, stride=stride, atom_indices=atom_indices)
        f = HDF5TrajectoryFile(self.trajectory_file)
        frames = f.read(n_frames=n_frames, stride=stride, atom_indices=atom_indices)
        self._energies = frames.potentialEnergy
        self._forces = frames.forces
        f.close()

    def torch_datasets(
            self,
            fields: Sequence[str] = ("coordinates",),
            shuffle: bool = True,
            random_seed: int = 1,
            val_fraction: float = 0.2,
            test_fraction: float = 0.0,
            process_array_fn: Callable = lambda x: torch.tensor(x)
    ) -> Sequence[torch.utils.data.TensorDataset]:
        """Get torch datasets for these data.

        Parameters
        ----------
            fields :
                A list of strings that correspond to attributes of this dataset.
            shuffle :
                Whether to shuffle each field
            random_seed :
                The random seed to be used for shuffling
            val_fraction :
                The fraction of data to be used as validation data
            test_fraction :
                The fraction of data to be used as test data
            process_array_fn :
                A function to turn the numpy array into a tensor

        Returns
        -------
            train_dataset :  torch.utils.data.TensorDataset
            validation_dataset : torch.utils.data.TensorDataset
            test_dataset : torch.utils.data.TensorDataset


        """
        if shuffle:
            rg_state = np.random.get_state()
            np.random.seed(random_seed)
            indices = np.random.permutation(len(self))
            np.random.set_state(rg_state)
        else:
            indices = np.arange(len(self))
        n_test = int(test_fraction * len(self))
        n_val = int(val_fraction * len(self))
        n_train = len(self) - n_test - n_val
        train_arrays = [getattr(self, key)[indices[:n_train]] for key in fields]
        val_arrays = [getattr(self, key)[indices[n_train:n_train + n_val]] for key in fields]
        test_arrays = [getattr(self, key)[indices[n_train+n_val:]] for key in fields]
        trainset = torch.utils.data.TensorDataset(*[process_array_fn(arr) for arr in train_arrays])
        valset = torch.utils.data.TensorDataset(*[process_array_fn(arr) for arr in val_arrays])
        testset = torch.utils.data.TensorDataset(*[process_array_fn(arr) for arr in test_arrays])
        return trainset, valset, testset
