

import os
import tempfile
from collections import OrderedDict
from simtk.openmm import app
from simtk import unit
from bgmol.systems import OpenMMSystem
from bgmol.util import get_data_file
from torchvision.datasets.utils import download_url

__all__ = ["MiniPeptide", "AMINO_ACIDS"]


class MiniPeptide(OpenMMSystem):
    """Small peptides (one or two amino acids, solvated or not). Created by Yaoyi Chen.

    Attributes
    ----------
    aminoacids : str
        A string of one-letter codes for aminoacids. All combinations of one or two aminoacids are supported.
    forcefield : list of str
        The builtin openmm force field files to be used
    constraints : app.internal.singleton.Singleton
        Constraint types
    solvated : bool
        Whether to add explicit water
    hydrogen_mass : unit.Quantity or None
        If None, don't repartition hydrogen mass. Else, assign the specified mass to hydrogen atoms.
    nonbonded_cutoff : unit.Quantity
        The cutoff for nonbonded forces.
    switch_distance : unit.Quantity or None
        Switch distance for nonbonded (LJ) interactions. If None, don't use a switch distance.
    root : str
        The root directory to which to download the files.
    download : bool
        Whether files should be downloaded.

    Notes
    -----
    Requires an internet connection to download the initial structure.
    """
    url = "http://ftp.mi.fu-berlin.de/pub/cmb-data/bgmol/systems/minipeptides/"

    def __init__(
            self,
            aminoacids="AA",
            forcefield=["amber99sbildn.xml", "tip3p.xml", "amber99_obc.xml"],
            constraints=app.HBonds,
            solvated=False,
            hydrogen_mass=None,
            nonbonded_cutoff=0.9*unit.nanometer,
            switch_distance=None,
            root=tempfile.gettempdir(),
            download=True
    ):
        super(MiniPeptide, self).__init__()

        # register parameters
        self.aminoacid = self.system_parameter(
            "aminoacids", aminoacids, default="AA"
        )
        self.forcefield = self.system_parameter(
            "forcefield", forcefield, default=["amber99sbildn.xml", "tip3p.xml", "amber99_obc.xml"]
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
            "nonbonded_cutoff", nonbonded_cutoff, default=0.9*unit.nanometer
        )
        self.switch_distance = self.system_parameter(
            "switch_distance", switch_distance, default=None
        )

        # create system
        if solvated:
            nonbonded_method = app.PME
            filename = f"{aminoacids}_solvated.pdb"
            forcefield = [fname for fname in forcefield if fname != "amber99_obc.xml"]
        else:
            nonbonded_method = app.NoCutoff
            filename = f"{aminoacids}.pdb"
            forcefield = [fname for fname in forcefield if fname != "tip3p.xml"]

        # download file
        full_filename = os.path.join(root, filename)
        if download:
            self._download(filename, root=root)
        else:
            assert os.path.isfile(full_filename)

        pdb = app.PDBFile(full_filename)
        ff = app.ForceField(*forcefield)
        self._system = ff.createSystem(
            pdb.topology,
            removeCMMotion=True,
            nonbondedMethod=nonbonded_method,
            nonbondedCutoff=nonbonded_cutoff,
            switchDistance=switch_distance,
            constraints=constraints,
            hydrogenMass=hydrogen_mass,
            rigidWater=True
        )
        self._positions = pdb.positions
        self._topology = pdb.topology

    def _download(self, filename, root):
        # get checksum
        md5 = None
        with open(get_data_file("../data/md5sums.txt"), "r") as fp:
            for line in fp:
                if line.split()[1] == filename:
                    md5 = line.split()[0]
        assert md5 is not None
        # download
        download_url(self.url+filename, root, filename, md5)
        assert os.path.isfile(os.path.join(root, filename))


AMINO_ACIDS = OrderedDict([
    ("A", "Alanine"),
    ("C", "Cysteine"),
    ("D", "AsparticAcid"),
    ("E", "GlutamicAcid"),
    ("F", "Phenylalanine"),
    ("G", "Glycine"),
    ("H", "Histidine"),
    ("I", "Isoleucine"),
    ("K", "Lysine"),
    ("L", "Leucine"),
    ("M", "Methionine"),
    ("N", "Asparagine"),
    ("P", "Proline"),
    ("Q", "Glutamine"),
    ("R", "Arginine"),
    ("S", "Serine"),
    ("T", "Threonine"),
    ("V", "Valine"),
    ("W", "Tryptophan"),
    ("Y", "Tyrosine")
])