
from bgmol.systems import TwoMiniPeptides, MiniPeptide
import numpy as np


def test_two_minipeptides():
    pep1 = MiniPeptide("A", solvated=True)
    pep2 = MiniPeptide("YG", solvated=True)
    twopep = TwoMiniPeptides("A", "YG")
    assert twopep.n_atoms == pep1.n_atoms + pep2.n_atoms
    assert isinstance(twopep.positions, np.ndarray)
    assert twopep.positions.shape == (twopep.n_atoms, 3)
    a,b,c = twopep.system.getDefaultPeriodicBoxVectors()
    assert a[0] > 1.8 * b[1]
