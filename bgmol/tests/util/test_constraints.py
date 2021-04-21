
import numpy as np
from simtk.openmm import LangevinIntegrator, Platform
from simtk.openmm.app import Simulation
from bgmol.util.constraints import project_forces_onto_constraints, _get_constraints, _get_masses


def force_correction(forces, positions, system, atom_indices):
    """manual force correction"""
    constraints = _get_constraints(system, atom_indices)
    masses = _get_masses(system, atom_indices)
    distance_vector = positions[:,constraints[:,1],:] - positions[:,constraints[:,0],:]
    # make unit vector
    distance_vector /= np.linalg.norm(distance_vector, axis=-1)[...,None]
    # project forces onto distance vector
    radial_force0 = np.sum(forces[:,constraints[:,0],:] * distance_vector, axis=-1)
    radial_force1 = np.sum(forces[:,constraints[:,1],:] * distance_vector, axis=-1)
    center_of_mass_radial_force = radial_force0 + radial_force1
    # both experience the same displacement; so distribute the force according to mass
    mass_sum = masses[constraints[:,0]] + masses[constraints[:,1]]
    mass_ratio0 = masses[constraints[:,0]]/mass_sum
    mass_ratio1 = masses[constraints[:,1]]/mass_sum
    force_correction0 = (mass_ratio0 * center_of_mass_radial_force - radial_force0)[...,None] * distance_vector
    force_correction1 = (mass_ratio1 * center_of_mass_radial_force - radial_force1)[...,None] * distance_vector
    return force_correction0, force_correction1


def test_ala2_1000(ala2dataset):
    projected_forces = project_forces_onto_constraints(
        ala2dataset.forces,
        ala2dataset.xyz,
        ala2dataset.system.system,
        ala2dataset.system.mdtraj_topology.select(ala2dataset.selection)
    )
    assert not np.allclose(ala2dataset.forces, projected_forces)
    assert (  # average change smaller than 10%
            np.abs(ala2dataset.forces - projected_forces).mean() /
            np.abs(ala2dataset.forces).mean()
    ) < 0.1

    # check that the force correction for the original forces is big
    correction1, correction2 = force_correction(
        ala2dataset.forces,
        ala2dataset.xyz,
        ala2dataset.system.system,
        ala2dataset.system.mdtraj_topology.select(ala2dataset.selection)
    )
    assert np.abs(correction1).max() > 10.
    assert np.abs(correction2).max() > 10.

    # check that the force correction for the projected forces is small
    correction3, correction4 = force_correction(
        projected_forces,
        ala2dataset.xyz,
        ala2dataset.system.system,
        ala2dataset.system.mdtraj_topology.select(ala2dataset.selection)
    )
    assert np.abs(correction3).max() < 1e-1
    assert np.abs(correction4).max() < 1e-1


