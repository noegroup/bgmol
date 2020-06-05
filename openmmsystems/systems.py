"""
Implementations of OpenMMSystems.
"""

from simtk import unit
from simtk.openmm import app
from openmmsystems.base import OpenMMSystem
from openmmsystems.util import get_data_file


class ImplicitBPTI(OpenMMSystem):
    """Aprotinin in implicit solvent."""
    def __init__(self, forcefield=['amber10.xml', 'amber10_obc.xml']):

        # call parent constructor
        super(ImplicitBPTI, self).__init__()

        # register parameters
        self.forcefield = self.system_parameter(
            "forcefield", forcefield, default=['amber10.xml', 'amber10_obc.xml']
        )

        # create system, positions, and topology
        pdb = app.PDBFile(get_data_file('data/bpti_top.pdb'))
        forcefield = app.ForceField(*forcefield)
        self._system = forcefield.createSystem(
            pdb.topology, removeCMMotion=False,
            nonbondedMethod=app.CutoffNonPeriodic,
            nonbondedCutoff=1.0*unit.nanometer,
            constraints=None, rigidWater=True
        )
        self._positions = pdb.positions
        self._topology = pdb.topology

