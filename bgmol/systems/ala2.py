import os
import tempfile
import numpy as np
from bgmol.util.importing import import_openmm
_, _, app = import_openmm()
import mdtraj as md
from ..systems.base import OpenMMToolsTestSystem, OpenMMSystem
from torchvision.datasets.utils import download_url


__all__ = ["AlanineDipeptideImplicit", "AlanineDipeptideTSF"]


def compute_phi_psi(traj):
    """Compute backbone dihedrals.

    Parameters
    ----------
    traj : mdtraj.Trajectory
    """
    phi_atoms = [4, 6, 8, 14]
    phi = md.compute_dihedrals(traj, indices=[phi_atoms])[:, 0]
    psi_atoms = [6, 8, 14, 16]
    psi = md.compute_dihedrals(traj, indices=[psi_atoms])[:, 0]
    return phi, psi


DEFAULT_RIGID_BLOCK = np.array([6, 8, 9, 10, 14])


DEFAULT_Z_MATRIX = np.array([
    [0, 1, 4, 6],
    [1, 4, 6, 8],
    [2, 1, 4, 0],
    [3, 1, 4, 0],
    [4, 6, 8, 14],
    [5, 4, 6, 8],
    [7, 6, 8, 4],
    [11, 10, 8, 6],
    [12, 10, 8, 11],
    [13, 10, 8, 11],
    [15, 14, 8, 16],
    [16, 14, 8, 6],
    [17, 16, 14, 15],
    [18, 16, 14, 8],
    [19, 18, 16, 14],
    [20, 18, 16, 19],
    [21, 18, 16, 19]
])


DEFAULT_GLOBAL_Z_MATRIX = np.row_stack([
    DEFAULT_Z_MATRIX,
    np.array([
        [9, 8, 6, 14],
        [10, 8, 14, 6],
        [6, 8, 14, -1],
        [8, 14, -1, -1],
        [14, -1, -1, -1]
    ])
])


class AlanineDipeptideTSF(OpenMMSystem):
    """Alanine Dipeptide from the Temperature-Steering Flows paper,
    Dibak, Klein, Noé (2020): https://arxiv.org/abs/2012.00429

    Notes
    -----
    Requires an internet connection to download the initial structure.
    """

    url = "http://ftp.mi.fu-berlin.de/pub/cmb-data/bgmol/systems/ala2/"

    def __init__(self,  root=tempfile.gettempdir(), download=True):
        super(AlanineDipeptideTSF, self).__init__()

        # download pdb file
        filename = "alanine-dipeptide-nowater.pdb"
        full_filename = os.path.join(root, filename)
        if download:
            download_url(self.url + filename, root, filename, md5="728635667ed4937cf4a0e5b7c801d9ea")
        assert os.path.isfile(full_filename)

        pdb = app.PDBFile(full_filename)
        ff = app.ForceField("amber99sbildn.xml", "amber96_obc.xml")
        self._system = ff.createSystem(
            pdb.getTopology(),
            removeCMMotion=True,
            nonbondedMethod=app.NoCutoff,
            constraints=app.HBonds,
            rigidWater=True
        )
        self._positions = pdb.getPositions(asNumpy=True)
        self._topology = pdb.getTopology()
        self.z_matrix = DEFAULT_Z_MATRIX.copy()
        self.rigid_block = DEFAULT_RIGID_BLOCK.copy()

    @staticmethod
    def compute_phi_psi(traj):
        return compute_phi_psi(traj)


class AlanineDipeptideImplicit(OpenMMToolsTestSystem):
    def __init__(self, constraints=app.HBonds, hydrogenMass=None):
        super().__init__("AlanineDipeptideImplicit")
        self.constraints = self.system_parameter("constraints", constraints, default=app.HBonds)
        self.hydrogenMass = self.system_parameter("hydrogenMass", hydrogenMass, default=None)
        self.z_matrix = DEFAULT_Z_MATRIX.copy()
        self.rigid_block = DEFAULT_RIGID_BLOCK.copy()

    @staticmethod
    def compute_phi_psi(traj):
        return compute_phi_psi(traj)

