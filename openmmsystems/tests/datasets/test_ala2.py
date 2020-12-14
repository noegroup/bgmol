
import pytest
from openmmsystems.datasets import Ala2Implicit1000, Ala2TSF300, Ala2TSF600, Ala2TSF1000


def test_ala2_300(tmpdir):
    dataset = Ala2Implicit1000(root=tmpdir, download=True, read=True)
    assert dataset.coordinates.shape == (len(dataset), 22, 3)
    assert dataset.forces.shape == (len(dataset), 22, 3)
    assert dataset.energies.shape == (len(dataset), )
    print(dataset.system)


@pytest.mark.slow
@pytest.mark.parametrize("Dataset", [Ala2TSF300, Ala2TSF600, Ala2TSF1000])
def test_ala2_tsf(tmpdir, Dataset):
    dataset = Dataset(root=tmpdir, download=True, read=True)
    assert dataset.coordinates.shape == (len(dataset), 22, 3)

