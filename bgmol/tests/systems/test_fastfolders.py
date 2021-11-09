import pytest

from bgmol.systems.fastfolders import FastFolder, FASTFOLDER_NAMES


@pytest.mark.parametrize("protein", list(FASTFOLDER_NAMES.keys()))
def test_fastfolder_systems(protein, tmpdir):
    FastFolder(protein, download=True, root=str(tmpdir))
