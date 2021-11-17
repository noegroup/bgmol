import pytest
import warnings
from bgmol.zmatrix import ZMatrixFactory
from bgflow import GlobalInternalCoordinateTransformation
from ..test_zmatrix import _check_trafo_complete, _check_bonds_physical


def test_fastfolder_systems(fastfolder_system):
    n_atoms = fastfolder_system.mdtraj_topology.n_atoms
    assert n_atoms == len(fastfolder_system.positions)


def test_fastfolder_energy(fastfolder_system):
    torch = pytest.importorskip("torch")
    n_atoms = fastfolder_system.mdtraj_topology.n_atoms
    pos = torch.tensor(fastfolder_system.positions.reshape(1, n_atoms * 3))
    fastfolder_system.energy_model.energy(pos)


def test_z_factory_fastfolders(fastfolder_system):
    if fastfolder_system.protein == "villin":
        pytest.skip("NLE not implemented")
    if fastfolder_system.solvated:
        pytest.skip("Not available for solvated systems")
    factory = ZMatrixFactory(fastfolder_system.mdtraj_topology)
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        zmatrix, _ = factory.build_with_templates()

    # make sure that all bonds in the zmatrix are physical bonds
    _check_bonds_physical(zmatrix, fastfolder_system.mdtraj_topology)

    # check trafo
    trafo = GlobalInternalCoordinateTransformation(zmatrix)
    _check_trafo_complete(trafo, fastfolder_system)
