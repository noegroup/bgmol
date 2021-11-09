import pytest

from bgmol.systems.fastfolders import FastFolder, FASTFOLDER_NAMES


@pytest.fixture(params=list(FASTFOLDER_NAMES.keys()), scope="session")
def fastfolder_system(request, tmpdir_factory):
    protein = request.param
    print(protein)
    tmpdir = tmpdir_factory.mktemp(protein)
    return FastFolder(protein, download=True, root=str(tmpdir))


def test_fastfolder_systems(fastfolder_system):
    n_atoms = fastfolder_system.mdtraj_topology.n_atoms
    assert n_atoms == len(fastfolder_system.positions)


@pytest.mark.slow
def test_fastfolder_energy(fastfolder_system):
    torch = pytest.importorskip("torch")
    n_atoms = fastfolder_system.mdtraj_topology.n_atoms
    pos = torch.tensor(fastfolder_system.positions.reshape(1, n_atoms * 3))
    fastfolder_system.energy_model.energy(pos)
