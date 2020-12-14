import os
import tempfile
import numpy as np
from simtk.openmm import app
import mdtraj as md
from openmmsystems.systems.base import OpenMMToolsTestSystem, OpenMMSystem
from openmmsystems.tpl.download import download_url
from openmmsystems.util import get_data_file


__all__ = ["AlanineDipeptideImplicit", "AlanineDipeptideTSF"]


class AlanineDipeptideTSF(OpenMMSystem):
    """Alanine Dipeptide from the Temperature-Steering Flows paper,
    Dibak, Klein, Noé (2020): https://arxiv.org/abs/2012.00429

    Notes
    -----
    Requires an internet connection to download the initial structure.
    """

    url = "ftp://ftp.mi.fu-berlin.de/pub/cmb-data/openmmsystems/ala2/"

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

    @property
    def rigid_block(self):
        return np.array([6, 8, 9, 10, 14])

    @property
    def z_matrix(self):
        return np.array([
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
            [21, 18, 16, 19],
        ])


class AlanineDipeptideImplicit(OpenMMToolsTestSystem):
    def __init__(self, constraints=app.HBonds, hydrogenMass=None):
        super(AlanineDipeptideImplicit, self).__init__("AlanineDipeptideImplicit")
        self.constraints = self.system_parameter("constraints", constraints, default=app.HBonds)
        self.hydrogenMass = self.system_parameter("hydrogenMass", hydrogenMass, default=None)

    @staticmethod
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

    @property
    def z_matrix(self):
        return np.array([
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
