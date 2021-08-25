import torch
from simtk import openmm as mm
from simtk import unit
import numpy as np

__all__ = [
    "bond_constraints", "bond_marginal_estimate", "angle_marginal_estimate",
    "bond_forces", "angle_forces", "torsion_forces", "bond_parameters",
    "angle_parameters", "torsion_parameters", "torsions", "constraint_parameters",
    "lookup_bonds", "lookup_angles"
]


def bond_constraints(system, coordinate_transform):
    """Parse constrained bonds.

    Parameters
    ----------
    coordinate_transform : Union[
        bgflow.GlobalInternalCoordinateTransformation,
        bgflow.MixedCoordinateTransformation,
        bgflow.RelativeInternalCoordinateTransformation
    ]

    Returns
    -------
    constrained_bond_indices : np.ndarray
        Array of shape (n_constrained_bonds, ).
        Indices of constrained bonds in the output of the coordinate transform.
    constrained_bond_lengths : np.ndarray
        Array of shape (n_constrained_bonds, ).
        Lengths of the constrained bonds in nm.
    """
    bonds = coordinate_transform.bond_indices
    lengths, force_constants = lookup_bonds(system, bonds, temperature=300.)
    constrained_bond_indices = np.where(np.isinf(force_constants))[0]
    constrained_bond_lenghts = lengths[constrained_bond_indices]
    return constrained_bond_indices, constrained_bond_lenghts


def bond_marginal_estimate(
        system,
        coordinate_transform,
        temperature,
        min_bond_length=0.01,
        max_bond_length=np.infty,
        device=torch.device("cpu"),
        dtype=torch.get_default_dtype()
):
    import bgflow as bg
    bonds = coordinate_transform.bond_indices
    lengths, force_constants = lookup_bonds(system, bonds, temperature=temperature)
    unconstrained_bond_indices = np.where(np.isfinite(force_constants))[0]
    lengths = lengths[unconstrained_bond_indices]
    force_constants = force_constants[unconstrained_bond_indices]
    sigma = 1.0/np.sqrt(force_constants)
    distribution = bg.TruncatedNormalDistribution(
        mu=torch.tensor(lengths, device=device, dtype=dtype),
        sigma=torch.tensor(sigma, device=device, dtype=dtype),
        lower_bound=torch.tensor(min_bond_length, device=device, dtype=dtype),
        upper_bound=torch.tensor(max_bond_length, device=device, dtype=dtype)
    )
    return distribution


def angle_marginal_estimate(
        system,
        coordinate_transform,
        temperature,
        min_angle=0.05,
        max_angle=None,
        device=torch.device("cpu"),
        dtype=torch.get_default_dtype()
):
    import bgflow as bg
    if max_angle is None:
        max_angle = 1.0 if coordinate_transform.normalize_angles else np.pi
    angles = coordinate_transform.angle_indices
    equilibria, force_constants = lookup_angles(system, angles, temperature=temperature)
    if coordinate_transform.normalize_angles:
        equilibria = equilibria / np.pi
        force_constants = force_constants * np.pi**2
    sigma = 1.0 / np.sqrt(force_constants)
    distribution = bg.TruncatedNormalDistribution(
        mu=torch.tensor(equilibria, device=device, dtype=dtype),
        sigma=torch.tensor(sigma, device=device, dtype=dtype),
        lower_bound=torch.tensor(min_angle, device=device, dtype=dtype),
        upper_bound=torch.tensor(max_angle, device=device, dtype=dtype)
    )
    return distribution


def bond_forces(system):
    for f in system.getForces():
        if isinstance(f, mm.HarmonicBondForce):
            yield f


def angle_forces(system):
    for f in system.getForces():
        if isinstance(f, mm.HarmonicAngleForce):
            yield f


def torsion_forces(system):
    for f in system.getForces():
        if isinstance(f, mm.PeriodicTorsionForce):
            yield f


def bond_parameters(system):
    for f in bond_forces(system):
        for i in range(f.getNumBonds()):
            yield f.getBondParameters(i)


def angle_parameters(system):
    for f in angle_forces(system):
        for i in range(f.getNumAngles()):
            yield f.getAngleParameters(i)


def torsion_parameters(system):
    for f in torsion_forces(system):
        for i in range(f.getNumTorsions()):
            yield f.getTorsionParameters(i)


def torsions(system):
    for t in torsion_parameters(system):
        yield (t[:4])
        yield (t[3::-1])


def constraint_parameters(system):
    for i in range(system.getNumConstraints()):
        yield system.getConstraintParameters(i)


def lookup_bonds(system, pairs, temperature):
    """Parse the equilibrium lengths and force constants of specified bonds
    from a simtk.openmm.System.

    Parameters
    ----------
    system : simtk.openmm.System
        The system object that contains all potential and constraint definitions.
    pairs : np.ndarray
        Atom ids of the bonds of shape (n_bonds_to_look_up, 2).
    temperature : float
        Temperature in Kelvin.

    Returns
    -------
    lengths : np.ndarray
        bond lengths in nanometers
    force_constants : np.ndarray
        dimensionless bond force constants in kBT
    """
    bondlength_lookup = {}
    for bond_params in bond_parameters(system):
        atom1, atom2, length, force_constant = bond_params
        thermodynamic_beta = 1.0 / (unit.constants.MOLAR_GAS_CONSTANT_R * temperature * unit.kelvin)
        bondlength_lookup[frozenset([atom1, atom2])] = (
            length.value_in_unit(unit.nanometer),
            (thermodynamic_beta * force_constant).value_in_unit(unit.nanometer**-2)
        )
    for constraint_params in constraint_parameters(system):
        atom1, atom2, length = constraint_params
        bondlength_lookup[frozenset([atom1, atom2])] = (
            length.value_in_unit(unit.nanometer),
            np.inf
        )
    lengths = np.array(
        [bondlength_lookup[frozenset([pair[0], pair[1]])][0] for pair in pairs]
    )
    force_constants = np.array(
        [bondlength_lookup[frozenset([pair[0], pair[1]])][1] for pair in pairs]
    )
    return lengths, force_constants


def lookup_angles(system, angles, temperature):
    angle_lookup = {}
    for angle_params in angle_parameters(system):
        atom1, atom2, atom3, angle, force_constant = angle_params
        a = angle.value_in_unit(unit.radian)
        thermodynamic_beta = 1.0 / (unit.constants.MOLAR_GAS_CONSTANT_R * temperature * unit.kelvin)
        k = (thermodynamic_beta *force_constant).value_in_unit(unit.radian**-2)
        angle_lookup[(atom1, atom2, atom3)] = (a ,k)
        angle_lookup[(atom3, atom2, atom1)] = (a ,k)
    equilibria = np.array(
        [angle_lookup[(angle[0], angle[1], angle[2])][0] for angle in angles]
    )
    force_constants = np.array(
        [angle_lookup[(angle[0], angle[1], angle[2])][1] for angle in angles]
    )
    return equilibria, force_constants

