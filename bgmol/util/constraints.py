
import numpy as np

__all__ = ["project_forces_onto_constraints"]


def project_forces_onto_constraints(forces, openmm_system, atom_indices):
    n_particles = len(atom_indices)
    assert n_particles == forces.shape[1]
    constraints = _get_constraints(openmm_system, atom_indices)
    masses = _get_masses(openmm_system, atom_indices)
    jacobian =  _build_constraint_jacobian(constraints, n_particles)
    matrix = jacobian @ np.diag(1/masses) @ jacobian.transpose()
    rhs = - jacobian @ np.diag(1/masses) @ forces
    correction_scale = np.linalg.solve(matrix, rhs)
    projected_forces = forces.copy()
    projected_forces += jacobian.transpose() @ correction_scale
    return projected_forces


def _build_constraint_jacobian(constraints, num_particles):
    jacobian = np.zeros((len(constraints), num_particles))
    jacobian[np.arange(len(constraints)), constraints[:,0]] = 2
    jacobian[np.arange(len(constraints)), constraints[:,1]] = -2
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

