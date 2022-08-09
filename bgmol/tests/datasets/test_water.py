import pytest

from bgmol.datasets.water import WaterDimerFlexibleTIP3P, WaterDimerRigidTIP3P


@pytest.mark.parametrize("Dataset", [WaterDimerRigidTIP3P, WaterDimerFlexibleTIP3P])
def test_water_dimer(tmpdir, Dataset):
    dimer = Dataset(root=str(tmpdir), read=True, download=True)
    assert dimer.coordinates.shape == (50000, 6, 3)
    assert dimer.forces.shape == (50000, 6, 3)
