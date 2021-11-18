import os
import tempfile
from ..systems.base import OpenMMSystem
from torchvision.datasets.utils import download_url

from ..util.importing import import_openmm
mm, unit, app = import_openmm()


__all__ = ["ArrestinActive"]


class ArrestinActive(OpenMMSystem):
    """
    """
    url = "http://ftp.mi.fu-berlin.de/pub/cmb-data/bgmol/systems/arrestin/"

    def __init__(self,  root=tempfile.gettempdir(), download=True):
        super(ArrestinActive, self).__init__()

        if download:
            filename = self._download(root)
        assert os.path.isfile(filename)

        pdb = app.PDBFile(filename)
        forcefield = app.ForceField("amber99sbildn.xml", "amber96_obc.xml")

        self._system = forcefield.createSystem(
            pdb.topology,
            nonbondedMethod=app.PME,
            nonbondedCutoff=1 * unit.nanometer,
            constraints=app.HBonds
        )

        self._positions = pdb.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        self._topology = pdb.getTopology()

    def _download(self, root):
        md5 = "c428fe87c0e888bc19e444eef05a9b3d"
        filename = "arr2-active_start.pdb"
        download_url(self.url+filename, root, filename, md5)
        local_file = os.path.join(root, filename)
        return local_file

