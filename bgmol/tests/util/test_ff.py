
import pytest
import numpy as np
from bgmol.util.ff import torsions, bond_parameters, bond_constraints
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
