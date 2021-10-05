
import numpy as np
from bgmol.util.importing import import_openmm
_, unit, _ = import_openmm()


__all__ = ["project_forces_onto_constraints"]


def project_forces_onto_constraints(forces, positions, openmm_system, atom_indices=None):
    """Project forces onto the tangent space of the constraint manifold.

    The background is that simulation engines usually apply the constraints to velocities and positions
    but not forces.

    Parameters
    ----------
    forces : np.ndarray
        Force array of shape (n_samples, n_selected_atoms, 3)
    positions : np.ndarray
        Position array of shape (n_samples, n_selected_atoms, 3)
    openmm_system : openmm.System
        The system that stores the constraint and mass information
    atom_indices : iterable of int
        The selected atom_ids. If None, assume that the force array contains forces for all atoms in the system.

    Returns
    -------
    projected_forces : np.ndarray
        Force array of shape (n_samples, n_selected_atoms, 3).
    """
    atom_indices = np.arange(openmm_system.getNumParticles()) if atom_indices is None else atom_indices
    n_particles = len(atom_indices)
    assert n_particles == forces.shape[1]
    assert n_particles == positions.shape[1]
    constraints = _get_constraints(openmm_system, atom_indices)
    masses = _get_masses(openmm_system, atom_indices)
    jacobian = _build_constraint_jacobian(constraints, positions)
    matrix = np.einsum("bcrx,r,bdrx -> bcd", jacobian, 1/masses, jacobian)
    rhs = - np.einsum("bcrx,r,brx->bc", jacobian, 1/masses, forces)
    correction_scale = np.linalg.solve(matrix, rhs)
    projected_forces = forces.copy()
    correction = np.einsum("bjrx,bj->brx", jacobian, correction_scale)
    return projected_forces + correction


def _build_constraint_jacobian(constraints, positions):
    batchsize = positions.shape[0]
    jacobian = np.zeros((batchsize, len(constraints), positions.shape[1], 3))
    distance_vector = positions[:, constraints[:,0], :] - positions[:, constraints[:,1], :]
    jacobian[:, np.arange(len(constraints)), constraints[:, 0], :] = 2 * distance_vector
    jacobian[:, np.arange(len(constraints)), constraints[:, 1], :] = -2 * distance_vector
    return jacobian


def _get_constraints(system, atom_indices):
    constraints = []
    indices = [i for i in atom_indices]
    for constraint in range(system.getNumConstraints()):
        p1, p2, distance = system.getConstraintParameters(constraint)
        if p1 in atom_indices and p2 in atom_indices:
            constraints.append([indices.index(p1), indices.index(p2)])
    return np.array(constraints)


def _get_masses(system, atom_indices):
    masses = []
    for particle in atom_indices:
        masses.append(system.getParticleMass(int(particle)).value_in_unit(unit.amu))
    return np.array(masses)
