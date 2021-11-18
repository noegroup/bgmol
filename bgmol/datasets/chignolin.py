
import os
import numpy as np

from .base import DataSet
from ..systems.chignolin import ChignolinC22Implicit

__all__ = ["ChignolinOBC2PT"]


class ChignolinOBC2PT(DataSet):
    """Chignolin miniprotein; parallel tempering in OBC2 implicit solvent.
    1600 ns with five temperatures
    `[250., 274.64013583, 301.70881683, 331.44540173, 364.11284061, 400.] == np.geomspace(250, 400, 6)`.
    Exchanges attempted once per ps.
    Langevin dynamics in single precision with 1/ps friction coefficient and 4 fs time step.
    Samples are spaced in 10 ps intervals. The dataset contains positions only.


    """
    url = "http://ftp.mi.fu-berlin.de/pub/cmb-data/bgmol/datasets/chignolin/ChignolinOBC2PT.tgz"
    md5 = "d1d5cd96414a5ab8915113a71a4a2325"
    num_frames = 160000
    size = 1826236
    selection = "all"
    openmm_version = "7.4.2"
    date = "2021/02/10"
    author = "Yaoyi Chen"
    temperatures = [250., 274.64013583, 301.70881683, 331.44540173, 364.11284061, 400.]

    def __init__(self, root=os.getcwd(), download: bool = False, read: bool = False,  temperature: float = 301.70881683):
        self._temperature_index = self._find_temperature_index(temperature)
        super().__init__(root=root, download=download, read=read)
        self._system = ChignolinC22Implicit()
        self._temperature = self.temperatures[self._temperature_index]

    def read(self, n_frames=None, stride=None, atom_indices=None):
        xyz = []
        for sequence in range(16):
            xyz.append(np.load(os.path.join(self.root, f"ChignolinOBC2PT/chi_gbsa_pt_100ns_{sequence}.npy")))
        xyz = np.concatenate(xyz, axis=1)
        self._xyz = xyz[self._temperature_index]

    def _find_temperature_index(self, temperature):
        isclose = np.isclose(temperature, self.temperatures, rtol=0.0, atol=1.0)
        if not np.any(isclose):
            raise ValueError(f"temperature has to be close to one of {ChignolinOBC2PT.temperatures}")
        return np.where(isclose)[0][0]