import warnings
from typing import Sequence

import torch
from ..util.importing import import_openmm
mm, unit, _ = import_openmm()
import numpy as np

__all__ = [
    "bond_constraints", "bond_marginal_estimate", "angle_marginal_estimate",
    "bond_forces", "angle_forces", "torsion_forces", "bond_parameters",
    "angle_parameters", "torsion_parameters", "torsions", "constraint_parameters",
    "lookup_bonds", "lookup_angles", "torsion_energies_from_ff_parameters",
    "torsion_marginal_cdf_estimate"
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
    """Estimate the marginal distribution of (unconstrained) bonds
    from force field parameters.

    Parameters
    ----------
    system : openmm.System
        The openmm system that contains all information about the energy function and constraints.
    coordinate_transform : bgflow CoordinateTransformation
    temperature : float
        Temperature in Kelvin.
    min_bond_length : float, optional
        The minimum bond length of the returned distribution.
    max_bond_length : float, optional
        The maximum bond length of the returned distribution.
    device : torch.device, optional
        The device on which the returned distribution operates.
    dtype : torch.dtype, optional
        The data type on which the returned distribution operates.

    Returns
    -------
    distribution : bgflow.TruncatedNormalDistribution
        The estimated bond marginal distribution.
        The dimension of this distribution is the number of unconstrained bonds in the coordinate transform.
    """
    import bgflow as bg
    bonds = coordinate_transform.bond_indices
    lengths, force_constants = lookup_bonds(system, bonds, temperature=temperature)
    unconstrained_bond_indices = np.where(np.isfinite(force_constants))[0]
    lengths = lengths[unconstrained_bond_indices]
    force_constants = force_constants[unconstrained_bond_indices]
    sigma = 1.0 / np.sqrt(force_constants)
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
    """Estimate the marginal distribution of angles
    from force field parameters.

    Parameters
    ----------
    system : openmm.System
        The openmm system that contains all information about the energy function and constraints.
    coordinate_transform : bgflow CoordinateTransformation
    temperature : float
        Temperature in Kelvin.
    min_angle : float, optional
        The minimum angle of the returned distribution.
    max_angle : float, optional
        The maximum angle of the returned distribution.
    device : torch.device, optional
        The device on which the returned distribution operates.
    dtype : torch.dtype, optional
        The data type on which the returned distribution operates.

    Returns
    -------
    distribution : bgflow.TruncatedNormalDistribution
        The estimated marginal distribution of angles.
        If the coordinate transform produces normalized angles in [0,1],
        the marginal distribution will produce normalized angles, too.
    """
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


def torsion_marginal_cdf_estimate(
        system,
        coordinate_transform,
        temperature,
        discrete_torsions=64,
        device=torch.device("cpu"),
        dtype=torch.get_default_dtype(),
        method="ff", # "scan",
        max_energy=1e3 # in kT,
):
    """discrete_torsions can be an int or an array of floats"""
    if isinstance(discrete_torsions, int):
        discrete_torsions = np.linspace(-np.pi, np.pi, discrete_torsions + 1)
    if isinstance(discrete_torsions, np.ndarray):
        if len(discrete_torsions.shape) == 1:
            discrete_torsions = np.tile(discrete_torsions, [coordinate_transform.dim_torsions, 1])
    assert len(discrete_torsions.shape) == 2
    assert np.allclose(discrete_torsions[:, 0], - np.pi * np.ones_like(discrete_torsions[:, 0]))
    assert np.allclose(discrete_torsions[:, -1], np.pi * np.ones_like(discrete_torsions[:, 0]))
    assert discrete_torsions.shape[0] == coordinate_transform.dim_torsions

    if method == "ff":
        energies = torsion_energies_from_ff_parameters(system, coordinate_transform, temperature, discrete_torsions)
    else:
        raise NotImplementedError(f"Method {method} not implemented for torsion_marginal_icdf_estimate")

    # clip and normalize energies, so that \int e^(-u) = 1
    energies = torch.tensor(energies)
    energies = energies - energies.min(dim=-1, keepdim=True).values
    energies = energies.clip(0, 20)
    marginal_free_energy = - torch.logsumexp(-energies[...,:-1], dim=-1, keepdim=True)
    energies -= marginal_free_energy
    # check that energies are normalized
    assert torch.allclose(torch.logsumexp(-energies[...,:-1], dim=-1), torch.zeros_like(energies[..., 0]), atol=1e-3)

    support_points = torch.tensor(discrete_torsions)
    if coordinate_transform.normalize_angles:
        support_points = (support_points + np.pi) / (2*np.pi)
    probabilities = torch.exp(- 0.5 * (energies[..., 1:] + energies[..., :-1]))  # midpoint rule
    probabilities = torch.cat([torch.zeros_like(probabilities[:,[0]]), probabilities], dim=-1)
    probabilities /= probabilities.sum(dim=-1, keepdim=True)
    support_values = torch.cumsum(probabilities, dim=-1)
    slopes = torch.exp(-energies)

    assert support_points.shape == (coordinate_transform.dim_torsions,  discrete_torsions.shape[-1])
    assert support_values.shape == (coordinate_transform.dim_torsions, discrete_torsions.shape[-1])
    assert slopes.shape == (coordinate_transform.dim_torsions,  discrete_torsions.shape[-1])

    from bgflow.nn.flow.spline import PeriodicTabulatedTransform
    cdf = PeriodicTabulatedTransform(support_points, support_values, slopes)
    return cdf.to(device=device, dtype=dtype)


def torsion_energies_from_ff_parameters(
        system,
        coordinate_transform,
        temperature,
        discrete_torsions,
):
    # discrete_torsions has shape (n_torsions, n_discrete_torsions)

    torsions = coordinate_transform.torsion_indices
    periodicities, phases, force_constants, chargeprods, sigmas, epsilons = lookup_torsions(system, torsions, temperature)

    # compute 1-4 distances
    d12, _ = lookup_bonds(system, torsions[:, :2], temperature=temperature)
    d23, _ = lookup_bonds(system, torsions[:, 1:3], temperature=temperature)
    d34, _ = lookup_bonds(system, torsions[:, 2:], temperature=temperature)
    a123, _ = lookup_angles(system, torsions[:, :3], temperature=temperature)
    a234, _ = lookup_angles(system, torsions[:, 1:], temperature=temperature)


    distances14 = _torsions_to_distances14(
        discrete_torsions,
        bonds=np.stack([d12, d23, d34], axis=-1),
        angles=np.stack([a123, a234], axis=-1)
    )

    energies = evaluate_torsion_potential(
        torsions=discrete_torsions,
        distances14=distances14,
        periodicities=periodicities,
        phases=phases,
        force_constants=force_constants,
        chargeprods=chargeprods,
        sigmas=sigmas,
        epsilons=epsilons
    )
    return energies


def _torsions_to_distances14(torsions, bonds, angles):
    from bgflow import GlobalInternalCoordinateTransformation
    with torch.no_grad():
        ictrafo = GlobalInternalCoordinateTransformation(
            z_matrix=np.array([[0, -1, -1, -1], [1, 0, -1, -1], [2, 1, 0, -1], [3, 2, 1, 0], ]),
            normalize_angles=False
        )
        bonds = torch.tensor(bonds)
        angles = torch.tensor(angles)
        torsions = torch.tensor(torsions)
        n_torsions = len(bonds)
        assert bonds.shape == (n_torsions, 3)
        assert angles.shape == (n_torsions, 2)
        assert len(torsions.shape) in {1, 2}
        n_discrete_torsions = torsions.shape[-1]
        if len(torsions.shape) == 1:
            # all torsions have the same discretization
            torsions = torsions[None, :, None].repeat(n_torsions, 1, 1)
        else:
            torsions = torsions[..., None]
        bonds = bonds[:, None, :].repeat(1, n_discrete_torsions, 1)
        angles = angles[:, None, :].repeat(1, n_discrete_torsions, 1)
        assert bonds.shape == (n_torsions, n_discrete_torsions, 3)
        assert angles.shape == (n_torsions, n_discrete_torsions, 2)
        assert torsions.shape == (n_torsions, n_discrete_torsions, 1)
        r, dlogp = ictrafo.forward(
            bonds.reshape(-1, 3),
            angles.reshape(-1, 2),
            torsions.reshape(-1, 1),
            x0=torch.zeros((n_discrete_torsions * n_torsions, 1, 3)),
            R=0.5*torch.ones((n_discrete_torsions * n_torsions, 3)),
            inverse=True
        )
        r = r.reshape(-1, 4, 3)
        distances14 = torch.linalg.norm(r[:, 0, :] - r[:, 3, :], dim=-1)
        return distances14.reshape(n_torsions, n_discrete_torsions).numpy()


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


def nonbonded_forces(system):
    for f in system.getForces():
        if isinstance(f, mm.NonbondedForce):
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
            p1, p2, p3, p4, period, phase, k = f.getTorsionParameters(i)
            yield p1, p2, p3, p4, period, phase, k
            yield p4, p3, p2, p1, period, phase, k


def torsions(system):
    for t in torsion_parameters(system):
        yield (t[:4])


def constraint_parameters(system):
    for i in range(system.getNumConstraints()):
        yield system.getConstraintParameters(i)


def exception_parameters(nonbonded_force : mm.NonbondedForce):
    assert nonbonded_force.getNumExceptionParameterOffsets() == 0
    for i in range(nonbonded_force.getNumExceptions()):
        p1, p2, *params = nonbonded_force.getExceptionParameters(i)
        yield (p1, p2, *params)
        yield (p2, p1, *params)


def nonbonded_parameters(
        system : mm.System,
        pairs : Sequence[Sequence[int]]
):
    # only implement it for one nonbonded force right now
    nbforces = list(nonbonded_forces(system))
    assert len(nbforces) <= 1
    if len(nbforces) == 0:
        return

    nbforce = nbforces[0]
    # non-mixed 1-4 interactions
    exceptions = {
        frozenset([p1, p2]): (qq, sigma, epsilon)
        for p1, p2, qq, sigma, epsilon in exception_parameters(nbforce)
    }
    for pair in pairs:
        pair = [int(pair[0]), int(pair[1])]
        if frozenset(pair) in exceptions:
            qq, sigma, epsilon = exceptions[frozenset(pair)]
        else:
            # mixing rules
            q1, s1, e1 = nbforce.getParticleParameters(pair[0])
            q2, s2, e2 = nbforce.getParticleParameters(pair[1])
            qq = q1*q2
            sigma = 0.5 * (s1 + s2)
            epsilon = np.sqrt(e1 * e2)
        yield qq, sigma, epsilon


def evaluate_torsion_potential(
        torsions, distances14, periodicities, phases,
        force_constants, chargeprods, sigmas, epsilons
):
    cosines = evaluate_raw_torsion_potential(torsions, periodicities, phases, force_constants)
    nonbonded_14 = evaluate_14_potential(distances14, chargeprods, sigmas, epsilons)
    assert cosines.shape == nonbonded_14.shape
    return cosines + nonbonded_14


def evaluate_raw_torsion_potential(torsions, periodicities, phases, force_constants):
    # torsions has shape (n_torsions_in_system, n_discrete_torsions)
    # periodicities, phases, force constants have shape (n_torsions_in_system, N_MAX_COSINE_TERMS)
    terms = (
            force_constants[:, None, :] *
            (1.0 + np.cos(periodicities[:, None, :] * torsions[:, :, None] - phases[:, None, :]))
    )
    # shape (n_torsions_in_system, n_discrete_torsions, N_MAX_COSINE_TERMS)
    # shape (... batch idxs ..., n_torsions, N_MAX_COSINE_TERMS, n_bins)
    return terms.sum(-1)


def evaluate_14_potential(distances, chargeprods, sigmas, epsilons):
    sig_by_r6 = (sigmas[...,None]/distances)**6
    lj14 = 4 * epsilons[...,None] * (sig_by_r6**2 - sig_by_r6)
    coulomb14 = chargeprods[...,None] / distances
    return lj14 + coulomb14


def lookup_bonds(system, pairs, temperature):
    """Parse the equilibrium lengths and force constants of specified bonds
    from a openmm.System.

    Parameters
    ----------
    system : openmm.System
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
        dimensionless bond force constants in kBT;
        any constrained bond is assigned an infinte force constant.
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
    lengths = []
    force_constants = []
    for pair in pairs:
        bond = frozenset([pair[0], pair[1]])
        if bond in bondlength_lookup:
            eq, k = bondlength_lookup[bond]
            lengths.append(eq)
            force_constants.append(k)
        else:
            warnings.warn(f"Bond {bond} not found in force field.", UserWarning)
            lengths.append(0.2)
            force_constants.append(1e-5)
    return np.array(lengths), np.array(force_constants)


def lookup_angles(system, angles, temperature):
    """Parse the equilibrium angles and force constants of specified angles
    from a openmm.System.

    Parameters
    ----------
    system : openmm.System
        The system object that contains all potential and constraint definitions.
    angles : np.ndarray
        Atom ids of the angles of shape (n_bonds_to_look_up, 3).
    temperature : float
        Temperature in Kelvin.

    Returns
    -------
    equilibria : np.ndarray
        equilibrium angles in radians
    force_constants : np.ndarray
        dimensionless force constants in kBT
    """
    angle_lookup = {}
    for angle_params in angle_parameters(system):
        atom1, atom2, atom3, angle, force_constant = angle_params
        a = angle.value_in_unit(unit.radian)
        thermodynamic_beta = 1.0 / (unit.constants.MOLAR_GAS_CONSTANT_R * temperature * unit.kelvin)
        k = (thermodynamic_beta * force_constant).value_in_unit(unit.radian**-2)
        angle_lookup[(atom1, atom2, atom3)] = (a, k)
        angle_lookup[(atom3, atom2, atom1)] = (a, k)
    equilibria = []
    force_constants = []
    for angle in angles:
        if (angle[0], angle[1], angle[2]) in angle_lookup:
            eq, k = angle_lookup[(angle[0], angle[1], angle[2])]
            equilibria.append(eq)
            force_constants.append(k)
        else:
            warnings.warn(f"Angle {angle} not found in force field.", UserWarning)
            equilibria.append(0.5)
            force_constants.append(1e-5)
    return np.array(equilibria), np.array(force_constants)


N_MAX_TORSION_TERMS = 6


def lookup_torsions(
        system,
        torsions,
        temperature
):
    """
    Parameters
    ----------
    system
    torsions
    temperature

    Returns
    -------

    """
    torsion_lookup = dict()
    nb_pairs = []
    for p1, p2, p3, p4, period, phase, k in torsion_parameters(system):
        torsion = (p1, p2, p3, p4)
        if torsion in torsion_lookup:
            torsion_lookup[torsion].append([period, phase, k])
        else:
            torsion_lookup[torsion] = [[period, phase, k]]
        nb_pairs.append([p1, p4])

    # cos-type interactions
    periodicities = np.zeros((len(torsions), N_MAX_TORSION_TERMS), dtype=int)
    phases = np.zeros((len(torsions), N_MAX_TORSION_TERMS))
    force_constants = np.zeros((len(torsions), N_MAX_TORSION_TERMS))

    thermodynamic_beta = 1.0 / (unit.constants.MOLAR_GAS_CONSTANT_R * temperature * unit.kelvin)
    for i, torsion in enumerate(torsions):
        parameters = torsion_lookup.get(tuple(torsion), [])
        for j, (period, phase, k) in enumerate(parameters):
            periodicities[i, j] = period
            phases[i, j] = phase.value_in_unit(unit.radian)
            force_constants[i, j] = thermodynamic_beta * k

    # 1-4 interactions
    chargeprods = np.zeros(len(torsions))
    sigmas = np.zeros(len(torsions))
    epsilons = np.zeros(len(torsions))
    ONE_4PI_EPS0 = 138.935456
    #vacuum_permittivity = (8.8541878128e-12 * unit.farad / unit.meter)

    nonbonds = nonbonded_parameters(system, [[p1, p4] for p1, _, _, p4, *_ in torsions])
    for i, (qq, sigma, epsilon) in enumerate(nonbonds):
        #qq = thermodynamic_beta * qq / (4 * np.pi * vacuum_permittivity) * unit.AVOGADRO_CONSTANT_NA
        qq = ONE_4PI_EPS0 * qq.value_in_unit(unit.elementary_charge**2) * unit.kilojoule_per_mole
        chargeprods[i] = thermodynamic_beta * qq
        sigmas[i] = sigma.value_in_unit(unit.nanometer)
        epsilons[i] = (thermodynamic_beta * epsilon)

    return periodicities, phases, force_constants, chargeprods, sigmas, epsilons

