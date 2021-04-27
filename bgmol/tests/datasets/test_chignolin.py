import pytest
import shutil
import numpy as np
from bgmol.datasets.chignolin import ChignolinOBC2PT


@pytest.mark.slow
def test_chignolin_obc_pt(tmpdir):
    """Try donwloading and reading at all temperatures."""
    ChignolinOBC2PT(root=str(tmpdir), read=True, download=True)
    with pytest.raises(ValueError):
        ChignolinOBC2PT(root=str(tmpdir), read=True, download=False, temperature=10000.)
    for temperature in ChignolinOBC2PT.temperatures:
        print(temperature)
        dataset = ChignolinOBC2PT(root=str(tmpdir), read=True, download=False, temperature=temperature)
        assert dataset.xyz.shape == (160000, 175, 3)
        assert np.isclose(dataset.temperature, temperature)
    shutil.rmtree(tmpdir)