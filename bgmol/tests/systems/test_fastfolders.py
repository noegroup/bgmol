import pytest
from itertools import product

from bgmol.systems.fastfolders import FastFolder, FAST_FOLDER_NAMES


@pytest.fixture(scope="session", params=product(
    FAST_FOLDER_NAMES,
    (False, True),
    ("charmm36m", ["amber99sbildn.xml", "tip3p.xml"])
))
def fastfolder_system(request, tmpdir_factory):
    protein, solvated, forcefield = request.param
    tmpdir = tmpdir_factory.mktemp(protein)
    if forcefield != "charmm36m" and protein in ["ntl9", "villin"]:
        pytest.skip("Nonstandard residues.")
    yield FastFolder(protein, download=True, solvated=solvated, forcefield=forcefield, root=str(tmpdir))


def test_fastfolder_systems(fastfolder_system):
    n_atoms = fastfolder_system.mdtraj_topology.n_atoms
    assert n_atoms == len(fastfolder_system.positions)


def test_fastfolder_energy(fastfolder_system):
    torch = pytest.importorskip("torch")
    n_atoms = fastfolder_system.mdtraj_topology.n_atoms
    pos = torch.tensor(fastfolder_system.positions.reshape(1, n_atoms * 3))
    fastfolder_system.reinitialize_energy_model(n_workers=1)
    fastfolder_system.energy_model.energy(pos)
