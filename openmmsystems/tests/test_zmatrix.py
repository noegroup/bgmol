import pytest
import warnings
import numpy as np
import torch
from simtk import unit
from openmmsystems.systems import ImplicitBPTI, AlanineDipeptideImplicit
from openmmsystems.zmatrix import ZMatrixFactory, build_fake_topology
from bgtorch import (
    RelativeInternalCoordinateTransformation,
    MixedCoordinateTransformation,
    GlobalInternalCoordinateTransformation
)


def _check_trafo_complete(trafo, system):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # we don't care if this is numerically sound for now
        xyz = torch.tensor(np.array(system.positions._value)).reshape(1,-1)
        trafo.to(xyz)
        *out, dlogp = trafo.forward(xyz)
        xyz2, dlogp2 = trafo.forward(*out, inverse=True)
        assert xyz.shape == xyz2.shape


@pytest.mark.parametrize("system", [
    #AlanineDipeptideImplicit(),
    ImplicitBPTI()
])
def test_z_matrix_mixed_trafo(system):
    top = system.mdtraj_topology
    factory = ZMatrixFactory(top, "backbone")
    zmatrix, fixed = factory.build_with_templates()
    positions = torch.tensor(system.positions.value_in_unit(unit.nanometer)).reshape(1,-1)
    data = positions.repeat(100, 1) + torch.randn(100, len(positions))
    trafo = MixedCoordinateTransformation(data, zmatrix, fixed, keepdims=10)
    _check_trafo_complete(trafo, system)


@pytest.mark.parametrize("system", [
    #AlanineDipeptideImplicit(),
    ImplicitBPTI()
])
def test_z_matrix_global_trafo(system):
    top = system.mdtraj_topology
    factory = ZMatrixFactory(top)
    zmatrix, _ = factory.build_with_templates()
    trafo = GlobalInternalCoordinateTransformation(zmatrix)
    _check_trafo_complete(trafo, system)


def test_z_factory_global():
    system = AlanineDipeptideImplicit()
    top = system.mdtraj_topology
    factory = ZMatrixFactory(top)
    z, _ = factory.build_naive()
    assert len(z) == top.n_atoms
    assert (np.sort(z[:,0]) == np.arange(top.n_atoms)).all()
    trafo = GlobalInternalCoordinateTransformation(z)
    _check_trafo_complete(trafo, system)


def test_z_factory_relative():
    system = AlanineDipeptideImplicit()
    top = system.mdtraj_topology
    factory = ZMatrixFactory(top, cartesian="element C")
    z, fixed = factory.build_naive()
    assert len(z) == top.n_atoms
    assert (np.sort(z[:,0]) == np.arange(top.n_atoms)).all()
    trafo = RelativeInternalCoordinateTransformation(z, fixed)
    _check_trafo_complete(trafo, system)


def test_z_factory_with_fake_topology():
    top, _ = build_fake_topology(20)
    factory = ZMatrixFactory(top, cartesian=[5, 10, 15])
    z, fixed = factory.build_naive()
    assert len(z) == top.n_atoms
    RelativeInternalCoordinateTransformation(z, fixed)