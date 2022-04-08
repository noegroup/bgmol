import pytest
import numpy as np
from bgmol.systems.ala2 import DEFAULT_GLOBAL_Z_MATRIX
from bgmol.systems import MiniPeptide, ChignolinC22Implicit
from bgmol.util import (
    rewire_chiral_torsions, find_rings, is_proper_torsion, is_chiral_torsion, is_ring_torsion, is_methyl_torsion,
    is_ramachandran_torsion, is_type_torsion
)
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
        ["HA", "CA", "N", "C"],  # improper
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
    z = ZMatrixFactory(top).build_naive(rewire_chiral=True)[0] if build_naive else DEFAULT_GLOBAL_Z_MATRIX
    z = rewire_chiral_torsions(z, top, verbose=False)
    indices = np.where(is_chiral_torsion(z, top))[0]
    assert len(indices) == 1
    for atom, name in zip(z[indices[0]], ("HA", "CA", "N", "C")):
        assert top.atom(atom).name == name


def test_rewire_torsions_valid_after():
    bgflow = pytest.importorskip("bgflow")
    top = ChignolinC22Implicit().mdtraj_topology
    z = ZMatrixFactory(top).build_naive(rewire_chiral=True)[0]
    #z = rewire_chiral_torsions(z, top, verbose=False)
    trafo = bgflow.GlobalInternalCoordinateTransformation(z)
    assert is_chiral_torsion(trafo.torsion_indices, top).sum() == 9  # one less for GLY


def test_is_ring_torsion():
    top = MiniPeptide("Y").mdtraj_topology
    z = ZMatrixFactory(top).build_with_templates()[0]
    is_ring = is_ring_torsion(z, top)
    assert is_ring.sum() == 5
    atoms = z[is_ring, 0]
    atomnames = [top.atom(i).name for i in atoms]
    assert set(atomnames) == {"CE2", "CE1", "CD1", "CD2", "CZ"}

    top = MiniPeptide("W").mdtraj_topology
    z = ZMatrixFactory(top).build_with_templates()[0]
    is_ring = is_ring_torsion(z, top)
    assert is_ring.sum() == 8

    top = ChignolinC22Implicit().mdtraj_topology
    z = ZMatrixFactory(top).build_with_templates()[0]
    is_ring = is_ring_torsion(z, top)
    assert is_ring.sum() == 4 + 3 * 5 + 8  # 1 PRO, 3 TYR, 1 TRP


def test_is_methyl_torsion():
    top = MiniPeptide("A").mdtraj_topology
    z = ZMatrixFactory(top).build_with_templates()[0]
    # there are three methyl groups in ala2
    assert is_methyl_torsion(z, top).sum() == 3


def test_is_ramachandran_torsion():
    top = MiniPeptide("A").mdtraj_topology
    z = ZMatrixFactory(top).build_with_templates()[0]
    # there are two ramachandran angles in ala2
    ramachandans = is_ramachandran_torsion(z, top)
    assert ramachandans.sum() == 2


@pytest.mark.parametrize("amino", ["A", "Y"])
def test_is_type_torsion(amino):
    top = MiniPeptide(amino).mdtraj_topology
    z = ZMatrixFactory(top).build_with_templates()[0]

    assert np.array_equal(is_type_torsion("ramachandran", z, top), is_ramachandran_torsion(z, top))
    assert np.array_equal(is_type_torsion("proper", z, top), is_proper_torsion(z, top))
    assert np.array_equal(is_type_torsion("methyl", z, top), is_methyl_torsion(z, top))
    assert np.array_equal(is_type_torsion("ring", z, top), is_ring_torsion(z, top))
    assert np.array_equal(is_type_torsion("chiral", z, top), is_chiral_torsion(z, top))
    assert is_type_torsion("phi", z, top).sum() == 1
    assert is_type_torsion("psi", z, top).sum() == 1
    assert is_type_torsion("omega", z, top).sum() == 0
    assert is_type_torsion('chi4', z, top).sum() == 0
