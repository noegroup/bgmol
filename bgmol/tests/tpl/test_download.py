
import os
from bgmol.tpl.download import download_url, download_and_extract_archive


def test_download(tmpdir):
    url = "ftp://ftp.mi.fu-berlin.de/pub/cmb-data/bgmol/README.md"
    download_url(url, tmpdir)
    assert os.path.isfile(tmpdir/"README.md")


def test_download_and_extract_archive(tmpdir):
    url = "ftp://ftp.mi.fu-berlin.de/pub/cmb-data/bgmol/test.tgz"
    download_and_extract_archive(url, tmpdir, tmpdir/"extracted")
    assert os.path.isfile(tmpdir/"extracted/README.md")
