import mdtraj as md
from .minipeptides import MiniPeptide


__all__ = ["AlanineTetrapeptideImplicit", "AlanineTetrapeptideVacuum"]


def compute_phi_psi(traj):
    """Compute backbone dihedrals.

    Parameters
    ----------
    traj : mdtraj.Trajectory
    """
    phi = md.compute_phi(traj)[1]
    psi = md.compute_psi(traj)[1]
    return phi, psi


class AlanineTetrapeptideImplicit(MiniPeptide):
    def __init__(self, **kwargs):
        super().__init__(aminoacids="AAA", **kwargs)

    @staticmethod
    def compute_phi_psi(traj):
        return compute_phi_psi(traj)


class AlanineTetrapeptideVacuum(MiniPeptide):
    def __init__(self, forcefield=['amber99sbildn.xml'], **kwargs):
        super().__init__(aminoacids="AAA", forcefield=forcefield, **kwargs)

    @staticmethod
    def compute_phi_psi(traj):
        return compute_phi_psi(traj)
