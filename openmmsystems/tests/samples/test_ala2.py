
from openmmsystems.samples import Ala2Implicit300


def test_dataset(tmpdir):
    dataset = Ala2Implicit300(root=tmpdir, download=True, read=True)
    assert dataset.coordinates.shape == (len(dataset), 22, 3)
    assert dataset.forces.shape == (len(dataset), 22, 3)
    assert dataset.energies.shape == (len(dataset), )

