import numpy as np
from simtk.openmm import app
import mdtraj as md
from openmmsystems.systems.base import OpenMMToolsTestSystem


__all__ = ["AlanineDipeptideImplicit"]


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
