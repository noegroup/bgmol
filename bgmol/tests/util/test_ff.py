

import pytest
import numpy as np
import torch
from bgmol.util.ff import (
    bond_parameters, bond_constraints, bond_marginal_estimate, angle_marginal_estimate,
    torsion_marginal_cdf_estimate, N_MAX_TORSION_TERMS, lookup_torsions,
    torsion_energies_from_ff_parameters, torsion_energies_from_scan, _ic_system
)
from bgmol.systems.ala2 import DEFAULT_GLOBAL_Z_MATRIX, DEFAULT_Z_MATRIX, DEFAULT_RIGID_BLOCK
from bgmol.systems.ala2 import AlanineDipeptideTSF
from bgflow import GlobalInternalCoordinateTransformation, MixedCoordinateTransformation
from bgmol.util.importing import import_openmm
mm, unit, _ = import_openmm()


def test_bond_parameters():
    assert len(list(bond_parameters(AlanineDipeptideTSF().system))) == 9


def test_bond_constraints():
    ala2 = AlanineDipeptideTSF()
    crd_trafo = GlobalInternalCoordinateTransformation(DEFAULT_GLOBAL_Z_MATRIX)
    indices, lengths = bond_constraints(ala2.system, crd_trafo)
    hydrogens = ala2.mdtraj_topology.select("element H")
    bonds_with_hydrogens = np.where(np.logical_or(
        np.isin(crd_trafo.bond_indices[:, 0], hydrogens),
        np.isin(crd_trafo.bond_indices[:, 1], hydrogens)
    ))[0]
    assert np.allclose(bonds_with_hydrogens, indices)
    assert lengths == pytest.approx(0.109 * np.ones_like(lengths), abs=1e-2)


def test_bond_marginal_estimate(ala2dataset, ctx):
    crd_trafo = GlobalInternalCoordinateTransformation(DEFAULT_GLOBAL_Z_MATRIX).to(**ctx)
    constrained_indices, constrained_lengths = bond_constraints(ala2dataset.system.system, crd_trafo)
    estimate = bond_marginal_estimate(ala2dataset.system.system, crd_trafo, 1000.).to(**ctx)
    samples = estimate.sample(10)
    mean = samples.mean(dim=0)
    sigma = samples.std(dim=0)
    unconstrained_indices = np.setdiff1d(np.arange(crd_trafo.dim_bonds), constrained_indices)
    constrained_lengths = torch.tensor(constrained_lengths, **ctx)

    bonds, *_ = crd_trafo.forward(torch.tensor(ala2dataset.xyz, **ctx).reshape(-1, 66))
    assert torch.allclose(bonds[:, unconstrained_indices].mean(dim=0), mean, atol=0.01)
    assert torch.allclose(bonds[:, unconstrained_indices].std(dim=0), sigma, atol=0.01)
    assert torch.allclose(bonds[:, constrained_indices].mean(dim=0),
                          constrained_lengths, atol=1e-5)
    assert torch.allclose(bonds[:, constrained_indices].std(dim=0),
                          torch.zeros_like(constrained_lengths), atol=1e-5)


@pytest.mark.parametrize("normalize", [True, False])
def test_angle_marginal_estimate(ala2dataset, ctx, normalize):
    crd_trafo = GlobalInternalCoordinateTransformation(
        DEFAULT_GLOBAL_Z_MATRIX, normalize_angles=normalize
    ).to(**ctx)
    estimate = angle_marginal_estimate(ala2dataset.system.system, crd_trafo, 1000.).to(**ctx)
    samples = estimate.sample(100)
    mean = samples.mean(dim=0)
    sigma = samples.std(dim=0)

    _, angles, *_ = crd_trafo.forward(torch.tensor(ala2dataset.xyz, **ctx).reshape(-1, 66))
    tol = 0.1 if normalize else np.pi * 0.1
    assert torch.allclose(angles.mean(dim=0), mean, atol=tol)
    assert torch.allclose(angles.std(dim=0), sigma, atol=tol)


def make_basic_system(n_particles=4, bond_k=1.0, sigma=0.01, epsilon=0.01):
    system = mm.System()
    for i in range(n_particles):
        system.addParticle(1.0)

    nonbonded_force = mm.NonbondedForce()
    bond_force = mm.HarmonicBondForce()
    angle_force = mm.HarmonicAngleForce()
    torsion_force = mm.PeriodicTorsionForce()
    for i in range(n_particles):
        nonbonded_force.addParticle(0.1, sigma, epsilon)
    for i in range(n_particles-1):
        bond_force.addBond(i, i+1, 0.1, bond_k)
    for i in range(n_particles-2):
        angle_force.addAngle(i, i+1, i+2, 1.0, 1.0)
    for i in range(n_particles-3):
        torsion_force.addTorsion(i, i+1, i+2, i+3, i+1, 0.3, 1.0)
    torsion_force.setForceGroup(10)
    nonbonded_force.setForceGroup(11)
    system.addForce(nonbonded_force)
    system.addForce(bond_force)
    system.addForce(angle_force)
    system.addForce(torsion_force)
    return system


def ground_truth_coulomb():
    s = make_basic_system(n_particles=2, bond_k=0., epsilon=0.0)
    c = mm.Context(s, mm.LangevinIntegrator(300, 0.001, 1.0))
    c.setPositions([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    e = c.getState(getEnergy=True).getPotentialEnergy()
    kbT = unit.MOLAR_GAS_CONSTANT_R * 300. * unit.kelvin
    e /= kbT
    return e


def test_torsion_energies():
    basic_system = make_basic_system()
    context = mm.Context(basic_system, mm.LangevinIntegrator(300, 0.001, 1.0))

    periods, phases, ks, qs, sigmas, epsilons = lookup_torsions(basic_system, [[0,1,2,3]], 300.)
    assert periods == pytest.approx(np.array([[1,] + [0]*(N_MAX_TORSION_TERMS-1)]))
    assert phases == pytest.approx(np.array([[0.3,] + [0]*(N_MAX_TORSION_TERMS-1)]))
    kbT = (unit.MOLAR_GAS_CONSTANT_R * 300. * unit.kelvin).value_in_unit(unit.kilojoule_per_mole)
    assert ks == pytest.approx(np.array([[1.0 / kbT,] + [0]*(N_MAX_TORSION_TERMS-1)]))
    assert qs == pytest.approx(np.array([ground_truth_coulomb(),]))
    assert sigmas == pytest.approx(np.array([0.01]))
    assert epsilons == pytest.approx([0.01 / kbT ])

    ictrafo = GlobalInternalCoordinateTransformation(
        z_matrix=np.array(
            [[0, -1, -1, -1], [1, 0, -1, -1],
             [2, 1, 0, -1], [3, 2, 1, 0],
             ], dtype=int),
        normalize_angles=False
    )

    def eval_mm(torsion_angle, context=context, ictrafo=ictrafo, groups={0,10,11}):
        with torch.no_grad():
            r, _ = ictrafo.forward(
                bonds=torch.tensor([[0.1]*ictrafo.dim_bonds]),
                angles=torch.tensor([[1.0]*ictrafo.dim_angles]),
                torsions=torch.tensor([[torsion_angle]*ictrafo.dim_torsions]),
                x0=torch.zeros((1, 1, 3)),
                R=torch.zeros((1, 3)),
                inverse=True
            )
            r = r.numpy()
        n_particles = basic_system.getNumParticles()
        context.setPositions(
            r.reshape(n_particles, 3)
        )
        beta = 1. / (300. * unit.kelvin * unit.MOLAR_GAS_CONSTANT_R)
        energy = beta * context.getState(getEnergy=True, groups=groups).getPotentialEnergy()
        return energy

    n_bins = 4
    discrete_torsions = np.linspace(-np.pi, np.pi, n_bins)
    #print(np.array([eval_mm(t, groups={10}) for t in discrete_torsions]))
    #print(np.array([eval_mm(t, groups={11}) for t in discrete_torsions]))
    mm_energies = np.array([eval_mm(t) for t in discrete_torsions])
    mm_energies -= mm_energies.min()
    energies = torsion_energies_from_ff_parameters(basic_system, ictrafo, 300., discrete_torsions[None, :])
    energies -= energies.min()
    assert energies.shape == (1, n_bins)
    assert np.allclose(energies[0], mm_energies)


@pytest.mark.parametrize("normalize", [True, False])
@pytest.mark.parametrize("method", ["scan", "ff"])
def test_torsion_marginal_estimate(method, ctx, normalize):
   # pytest.importorskip("bgflow.nn.flow.spline.PeriodicTabulatedTransform")
    ala2 = AlanineDipeptideTSF()
    positions = torch.tensor(ala2.positions, **ctx)
    crd_trafo = GlobalInternalCoordinateTransformation(
        DEFAULT_GLOBAL_Z_MATRIX, normalize_angles=normalize
    ).to(**ctx)
    estimate = torsion_marginal_cdf_estimate(ala2.system, crd_trafo, 300., method=method, coordinates=positions)
    estimate = estimate.to(**ctx)

    _, _, torsions, *_ = crd_trafo.forward(positions.reshape(-1, 66))
    tol = 0.1 if normalize else np.pi * 0.1
    uniform = estimate.forward(torsions)[0]
    assert torch.allclose(uniform.mean(dim=0), 0.5 * torch.ones_like(uniform[0]), atol=0.35)
    assert torch.all(uniform.min(dim=0).values >= 0.)
    assert torch.all(uniform.max(dim=0).values <= 1.)


def test_ic_system():
    ala2 = AlanineDipeptideTSF()
    icsystem = _ic_system(ala2.system)
    # make sure the nonbonded forces have been deleted
    for f in icsystem.getForces():
        assert not isinstance(f, mm.NonbondedForce)
        assert not isinstance(f, mm.CustomGBForce)
    # make sure the 1-4 force has been created
    assert any(isinstance(f, mm.CustomBondForce) for f in icsystem.getForces())


@pytest.mark.parametrize("normalize_angles", (False, True))
def test_torsion_scan(normalize_angles, ala2dataset, **ctx):
    # TODO: bugfix and enable for normalized_angles
    ala2 = ala2dataset.system
    data = torch.tensor(ala2dataset.xyz, **ctx).reshape(-1, 22)
    num_bins = 24
    trafo = MixedCoordinateTransformation(
        data=torch.tensor(ala2dataset.xyz, **ctx),
        z_matrix=DEFAULT_Z_MATRIX,
        fixed_atoms=DEFAULT_RIGID_BLOCK,
        normalize_angles=normalize_angles
    )
    energies, bins = torsion_energies_from_scan(ala2.system, trafo, data[0], torchnum_bins=num_bins)
    assert energies.shape == (len(DEFAULT_Z_MATRIX), num_bins + 1)
    # assert continuity around the periodic boundary
    assert torch.allclose(energies[:, 0], energies[:, -1])