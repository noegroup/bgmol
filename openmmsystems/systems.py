
from simtk import unit
from simtk.openmm import app
from openmmsystems.base import OpenMMSystem


class ImplicitBPTI(OpenMMSystem):
    def __init__(self, forcefield=['amber99sbildn.xml', 'amber99_obc.xml']):
        # call parent constructor
        super(self, ImplicitBPTI).__init__()
        # register parameters
        self.system_parameter("forcefield", forcefield, default=['amber99sbildn.xml', 'amber99_obc.xml'])

        # create system, positions, and topology
        pdb = app.PDBFile('bpti_top.pdb')
        forcefield = app.ForceField(*forcefield)
        self._system = forcefield.createSystem(
            pdb.topology, removeCMMotion=False,
            nonbondedMethod=app.CutoffNonPeriodic,
            nonbondedCutoff=1.0*unit.nanometer,
            constraints=None, rigidWater=True
        )
        self._positions = pdb.positions
        self._topology = pdb.topology

