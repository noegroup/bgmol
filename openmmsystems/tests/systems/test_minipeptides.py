import pytest

from openmmsystems.systems import MiniPeptide, AMINO_ACIDS


@pytest.mark.slow
@pytest.mark.parametrize("aminoacid", AMINO_ACIDS)
@pytest.mark.parametrize("solvated", [True, False])
def test_monopeptide_systems(aminoacid, solvated):
    MiniPeptide(aminoacid, solvated=solvated)


@pytest.mark.slow
@pytest.mark.parametrize("aminoacid1", AMINO_ACIDS)
@pytest.mark.parametrize("aminoacid2", AMINO_ACIDS)
@pytest.mark.parametrize("solvated", [True, False])
def test_dipeptide_systems(aminoacid1, aminoacid2, solvated):
    MiniPeptide(f"{aminoacid1}{aminoacid2}", solvated=solvated)

