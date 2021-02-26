
import pytest

import mdtraj as md
from openmmsystems.datasets import Ala2Implicit1000, Ala2TSF300, Ala2TSF600, Ala2TSF1000
from openmmsystems.systems.ala2 import DEFAULT_Z_MATRIX, DEFAULT_RIGID_BLOCK


def test_ala2_1000(tmpdir):
    dataset = Ala2Implicit1000(root=tmpdir, download=True, read=True)
    assert dataset.coordinates.shape == (len(dataset), 22, 3)
    assert dataset.forces.shape == (len(dataset), 22, 3)
    assert dataset.energies.shape == (len(dataset), )
    # test api
    assert isinstance(dataset.trajectory, md.Trajectory)
    dataset.system.compute_phi_psi(dataset.trajectory)
    assert (dataset.system.z_matrix == DEFAULT_Z_MATRIX).all()
    assert (dataset.system.rigid_block == DEFAULT_RIGID_BLOCK).all()


@pytest.mark.slow
@pytest.mark.parametrize("Dataset", [Ala2TSF300, Ala2TSF600, Ala2TSF1000])
def test_ala2_tsf(tmpdir, Dataset):
    dataset = Dataset(root=tmpdir, download=True, read=True)
    assert dataset.coordinates.shape == (len(dataset), 22, 3)
    # test api
    assert isinstance(dataset.trajectory, md.Trajectory)
    dataset.system.compute_phi_psi(dataset.trajectory)
    assert (dataset.system.z_matrix == DEFAULT_Z_MATRIX).all()
    assert (dataset.system.rigid_block == DEFAULT_RIGID_BLOCK).all()


