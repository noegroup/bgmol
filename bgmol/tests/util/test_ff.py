

import pytest
import numpy as np
import torch
from bgmol.util.ff import bond_parameters, bond_constraints, bond_marginal_estimate, angle_marginal_estimate
from bgmol.systems.ala2 import DEFAULT_GLOBAL_Z_MATRIX
from bgmol.systems.ala2 import AlanineDipeptideTSF
from bgflow import GlobalInternalCoordinateTransformation


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
