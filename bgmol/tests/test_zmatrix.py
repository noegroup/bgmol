import pytest
import warnings
import numpy as np
import torch
from bgmol.util.importing import import_openmm
_, unit, _ = import_openmm()
from bgmol.systems import ImplicitBPTI, AlanineDipeptideImplicit, ChignolinC22Implicit
from bgmol.zmatrix import ZMatrixFactory, build_fake_topology
from bgflow import (
    RelativeInternalCoordinateTransformation,
    MixedCoordinateTransformation,
    GlobalInternalCoordinateTransformation
)


def _check_trafo_complete(trafo, system):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # we don't care if this is numerically sound for now
        xyz = torch.tensor(np.array(system.positions)).reshape(1, -1)
        trafo.to(xyz)
        *out, dlogp = trafo.forward(xyz)
        xyz2, dlogp2 = trafo.forward(*out, inverse=True)
        assert xyz.shape == xyz2.shape


@pytest.mark.parametrize("system", [
    AlanineDipeptideImplicit(),
    ImplicitBPTI()
])
def test_z_matrix_mixed_trafo(system):
    top = system.mdtraj_topology
    factory = ZMatrixFactory(top, "backbone")
    zmatrix, fixed = factory.build_with_templates()
    positions = torch.tensor(system.positions).reshape(1, -1)
    data = positions.repeat(100, 1) + torch.randn(100, len(positions))
    trafo = MixedCoordinateTransformation(data, zmatrix, fixed, keepdims=10)
    _check_trafo_complete(trafo, system)


@pytest.mark.parametrize("system", [
    AlanineDipeptideImplicit(),
    ChignolinC22Implicit(),
    ImplicitBPTI(),
])
def test_z_matrix_global_trafo(system):
    top = system.mdtraj_topology
    factory = ZMatrixFactory(top)
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        zmatrix, _ = factory.build_with_templates()
    trafo = GlobalInternalCoordinateTransformation(zmatrix)
    _check_trafo_complete(trafo, system)


def test_z_factory_global():
    system = AlanineDipeptideImplicit()
    top = system.mdtraj_topology
    factory = ZMatrixFactory(top)
    z, _ = factory.build_naive()
    assert factory.is_independent(z)
    assert len(z) == top.n_atoms
    trafo = GlobalInternalCoordinateTransformation(z)
    _check_trafo_complete(trafo, system)


def test_z_factory_relative():
    system = AlanineDipeptideImplicit()
    top = system.mdtraj_topology
    factory = ZMatrixFactory(top, cartesian="element C")
    z, fixed = factory.build_naive()
    assert len(z) + len(fixed) == top.n_atoms
    #assert (np.sort(z[:,0]) == np.arange(top.n_atoms)).all()
    trafo = RelativeInternalCoordinateTransformation(z, fixed)
    _check_trafo_complete(trafo, system)


def test_z_factory_with_fake_topology():
    top, _ = build_fake_topology(20)
    factory = ZMatrixFactory(top, cartesian=[5, 10, 15])
    z, fixed = factory.build_naive()
    assert len(z) + len(fixed) == top.n_atoms
    RelativeInternalCoordinateTransformation(z, fixed)


def test_z_factory_naive_cgn(chignolin):
    factory = ZMatrixFactory(chignolin.mdtraj_topology, cartesian="name == CA")
    z, fixed = factory.build_naive()
    assert fixed.shape == (10, )
    assert z.shape == (165, 4)

