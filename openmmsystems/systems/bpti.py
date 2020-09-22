from simtk import unit
from simtk.openmm import app
from openmmsystems.systems.base import OpenMMSystem
from openmmsystems.util import get_data_file


__all__ = ["ImplicitBPTI"]


class ImplicitBPTI(OpenMMSystem):
    """Aprotinin in implicit solvent."""
    def __init__(self, forcefield=['amber10.xml', 'amber10_obc.xml'], constraints=app.HBonds):

        # call parent constructor
        super(ImplicitBPTI, self).__init__()

        # register parameters
        self.forcefield = self.system_parameter(
            "forcefield", forcefield, default=['amber10.xml', 'amber10_obc.xml']
        )
        self.constraints = self.system_parameter(
            "constraints", constraints, default=app.HBonds
        )

        # create system, positions, and topology
        pdb = app.PDBFile(get_data_file('../data/bpti_top.pdb'))
        forcefield = app.ForceField(*forcefield)
        self._system = forcefield.createSystem(
            pdb.topology, removeCMMotion=False,
            nonbondedMethod=app.CutoffNonPeriodic,
            nonbondedCutoff=1.0*unit.nanometer,
            constraints=constraints, rigidWater=True
        )
        self._positions = pdb.positions
        self._topology = pdb.topology

