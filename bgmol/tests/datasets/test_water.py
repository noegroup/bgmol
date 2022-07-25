
from bgmol.datasets.water import DimerFlexibleTIP3P


def test_water_dimer(tmpdir):
    dimer = DimerFlexibleTIP3P(root=str(tmpdir), read=True, download=True)
    assert dimer.coordinates.shape == (50000, 6, 3)
    assert dimer.forces.shape == (50000, 6, 3)
