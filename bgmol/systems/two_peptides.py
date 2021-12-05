
import numpy as np

from .base import OpenMMSystem
from .minipeptides import MiniPeptide
from ..util.importing import import_openmm

_, unit, app = import_openmm()


class TwoMiniPeptides(OpenMMSystem):
    def __init__(
            self,
            aminoacids1,
            aminoacids2,
            forcefield=["amber99sbildn.xml", "tip3p.xml"],
            constraints=app.HBonds,
            solvated=True,
            hydrogen_mass=None,
            nonbonded_cutoff=0.9 * unit.nanometer,
            switch_distance=0.75 * unit.nanometer,
            **kwargs
    ):
        super().__init__()
        # register parameters
        self.forcefield = self.system_parameter(
            "forcefield", forcefield, default=["amber99sbildn.xml", "tip3p.xml"]
        )
        self.constraints = self.system_parameter(
            "constraints", constraints, default=app.HBonds
        )
        self.solvated = self.system_parameter(
            "solvated", solvated, default=True
        )
        self.hydrogen_mass = self.system_parameter(
            "hydrogen_mass", hydrogen_mass, default=None
        )
        self.nonbonded_cutoff = self.system_parameter(
            "nonbonded_cutoff", nonbonded_cutoff, default=0.9 * unit.nanometer
        )
        self.switch_distance = self.system_parameter(
            "switch_distance", switch_distance, default=0.75 * unit.nanometer
        )

        system1 = MiniPeptide(aminoacids1, solvated=solvated)
        system2 = MiniPeptide(aminoacids2, solvated=solvated)
        offset = system1.positions[:, 0].max() - system2.positions[:, 0].min()
        positions1 = self._center(system1.positions, mean=np.array([-offset / 2, 0.0, 0.0]))
        positions2 = self._center(system2.positions, mean=np.array([offset / 2, 0.0, 0.0]))
        modeller = app.Modeller(system1.topology, positions1 * unit.nanometer)
        modeller.add(system2.topology, positions2 * unit.nanometer)

        forcefield = app.ForceField(*forcefield)
        nonbonded_method = app.PME if solvated else app.CutoffNonPeriodic
        self._system = forcefield.createSystem(
            modeller.getTopology(),
            removeCMMotion=True,
            nonbondedMethod=nonbonded_method,
            nonbondedCutoff=self.nonbonded_cutoff,
            switchDistance=self.switch_distance,
            constraints=self.constraints,
            hydrogenMass=self.hydrogen_mass,
            rigidWater=True
        )
        self._positions = np.row_stack(modeller.getPositions().value_in_unit(unit.nanometers))
        self._topology = modeller.getTopology()

    @staticmethod
    def _center(positions, mean=np.array([0.0, 0.0, 0.0])):
        return positions - positions.mean() + mean
