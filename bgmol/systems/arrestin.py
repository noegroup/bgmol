import os
import tempfile
import numpy as np
from simtk.openmm import app
import mdtraj as md
from ..systems.base import OpenMMToolsTestSystem, OpenMMSystem
from torchvision.datasets.utils import download_url
import pickle5 as pickle


__all__ = ["ArrestinActive", "ArrestinInactive"]


class ArrestinActive(OpenMMSystem):
    """
    """
    def __init__(self,  root=tempfile.gettempdir(), download=True):
        super(ArrestinActive, self).__init__()

        filename = "openmm-arrestin/arr2-active_start.pdb"

        pdb = app.PDBFile(filename)
        ff = app.ForceField("amber99sbildn.xml", "amber96_obc.xml")

        name = 'openmm-arrestin/arr2-active_system'
        with open(name + '.pkl', 'rb') as f:
            system = pickle.load(f)
            self._system = system

        self._positions = pdb.getPositions(asNumpy=True)
        self._topology = pdb.getTopology()
        self.z_matrix = None
        self.rigid_block = None


class ArrestinInactive(OpenMMSystem):
    """
    """
    def __init__(self,  root=tempfile.gettempdir(), download=True):
        super(ArrestinInactive, self).__init__()

        filename = "openmm-arrestin/arr2-inactive_system"
        full_filename = os.path.join(root, filename)

        pdb = app.PDBFile(filename)
        ff = app.ForceField("amber99sbildn.xml", "amber96_obc.xml")

        name = 'arrestin_system_inactive'#file doesn't exist yet
        with open(name + '.pkl', 'rb') as f:
            system = pickle.load(f)
            self._system = system
        
        self._positions = pdb.getPositions(asNumpy=True)
        self._topology = pdb.getTopology()
        self.z_matrix = None
        self.rigid_block = None


