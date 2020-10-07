
from openmmsystems.datasets import Ala2Implicit1000


def test_ala2_300(tmpdir):
    dataset = Ala2Implicit1000(root=tmpdir, download=True, read=True)
    assert dataset.coordinates.shape == (len(dataset), 22, 3)
    assert dataset.forces.shape == (len(dataset), 22, 3)
    assert dataset.energies.shape == (len(dataset), )
    print(dataset.system)


