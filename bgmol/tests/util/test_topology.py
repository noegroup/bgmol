import pytest
from bgmol.systems.ala2 import DEFAULT_GLOBAL_Z_MATRIX
from bgmol.systems import MiniPeptide
from bgmol.util import rewire_chiral_torsions, find_rings, is_proper_torsion
from bgmol.zmatrix import ZMatrixFactory


@pytest.mark.parametrize("aminoacid, n_rings, n_ring_atoms", [
    ("A", 0, {}),
    ("C", 0, {}),
    ("D", 0, {}),
    ("E", 0, {}),
    ("F", 1, {6}),
    ("G", 0, {}),
    ("H", 1, {5}),
    ("I", 0, {}),
    ("K", 0, {}),
    ("L", 0, {}),
    ("M", 0, {}),
    ("N", 0, {}),
    ("P", 1, {5}),
    ("Q", 0, {}),
    ("R", 0, {}),
    ("S", 0, {}),
    ("T", 0, {}),
    ("V", 0, {}),
    ("W", 2, {5,6}),
    ("Y", 1, {6})
])
def test_find_rings(aminoacid, n_rings, n_ring_atoms):
    rings = find_rings(MiniPeptide(aminoacid).mdtraj_topology)
    assert len(rings) == n_rings
    if n_rings > 0:
        assert {len(ring) for ring in rings} == n_ring_atoms
        for ring in rings:
            assert isinstance(ring, list)
            for index in ring:
                assert isinstance(index, int)


def test_proper_torsions(ala2dataset):
    top = ala2dataset.system.mdtraj_topology
    torsions = (
        ["CB", "CA", "N", "C"],  # improper
        ["HA", "CA", "N", "CB"],  # improper
        ["O", "C", "CA", "CB"],  # proper
        ["H", "N", "CA", "HA"],  # proper
    )
    torsion_indices = [
        [top.select(f"resname ALA and name {name}")[0] for name in torsion]
        for torsion in torsions
    ]
    is_proper = is_proper_torsion(torsion_indices, top)
    assert is_proper.tolist() == [False, False, True, True]


@pytest.mark.parametrize("build_naive", [True, False])
def test_rewire_torsions(build_naive):
    top = MiniPeptide("A").mdtraj_topology
    zmatrix = ZMatrixFactory(top).build_naive()[0] if build_naive else DEFAULT_GLOBAL_Z_MATRIX
    zmat, indices = rewire_chiral_torsions(zmatrix, top)
    assert len(indices) == 1
    for atom, name in zip(zmat[indices[0]], ("CB", "CA", "N", "C")):
        assert top.atom(atom).name == name
