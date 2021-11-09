
from bgmol.systems.ala2 import DEFAULT_GLOBAL_Z_MATRIX
from bgmol.util import rewire_chiral_torsions


def test_find_chiral_torsions(ala2dataset):
    from bgflow import GlobalInternalCoordinateTransformation
    top = ala2dataset.system.mdtraj_topology
    zmatrix = DEFAULT_GLOBAL_Z_MATRIX
    torsion_indices = rewire_chiral_torsions(zmatrix, top)
    print(torsion_indices)
    print(zmatrix[17])
    assert len(torsion_indices) == 1
    assert zmatrix[torsion_indices, 0] == top.select("name HA")
    crd_trafo = GlobalInternalCoordinateTransformation(zmatrix)
