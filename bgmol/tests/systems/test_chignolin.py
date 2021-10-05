import pytest

from bgmol.util.importing import import_openmm
mm, _, app = import_openmm()


def test_chignolin_tics(chignolin):
    tics = chignolin.to_tics(chignolin.positions[None, ...], eigs_kept=10)
    assert tics.shape == (1, 10)


def test_chignolin_top(chignolin):
    assert chignolin.topology.getNumAtoms() == 175
    assert chignolin.topology.getNumBonds() >= 175
