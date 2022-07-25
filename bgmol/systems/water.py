

from .base import OpenMMSystem
from ..util.importing import import_openmm
from .base import OpenMMToolsTestSystem

_, unit, _ = import_openmm()


class WaterDimer(OpenMMSystem):
    def __init__(
            self,
            harmonic_force_constant=3 * unit.kilojoule_per_mole / unit.nanometer ** 2,
            constrained=True
    ):

        super().__init__()



