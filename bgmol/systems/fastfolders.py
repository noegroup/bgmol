

import os
import tempfile
from bgmol.util.importing import import_openmm
_, unit, app = import_openmm()
from bgmol.systems import OpenMMSystem
from bgmol.util import get_data_file
from torchvision.datasets.utils import download_url

__all__ = ["FASTFOLDER_NAMES", "FastFolder"]


class FastFolder(OpenMMSystem):
    """Small, fast-folding proteins.

    This is a remake of the proteins from the paper Created by Yaoyi Chen.

    Attributes
    ----------
    protein : str
        The name of the protein. See the dictionary `FASTFOLDER_NAMES` in this module for valid options.
    forcefield : list of str
        The builtin openmm force field files to be used
    constraints : app.internal.singleton.Singleton
        Constraint types
    solvated : bool
        Whether to add explicit water and ions.
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

    References
    ----------
    [1] K. Lindorff-Larsen, S. Piana, R.O. Dror, D.E. Shaw,
        How fast-folding proteins fold,
        Science (80-. ). 334 (2011) 517–520.
        doi:10.1126/science.1208351.

    Notes
    -----
    Requires an internet connection to download the topology and initial structure.
    """
    url = "http://ftp.mi.fu-berlin.de/pub/cmb-data/bgmol/systems/fastfolders/"

    def __init__(
            self,
            protein="chignolin",
            forcefield=["amber99sbildn.xml", "tip3p.xml", "amber99_obc.xml"],
            constraints=app.HBonds,
            solvated=False,
            hydrogen_mass=None,
            nonbonded_cutoff=0.9 * unit.nanometer,
            switch_distance=None,
            root=tempfile.gettempdir(),
            download=True
    ):
        super().__init__()

        if not protein in FASTFOLDER_NAMES:
            raise ValueError(f"protein has to be one of {list(FASTFOLDER_NAMES.keys())}")

        # register parameters
        self.protein = self.system_parameter(
            protein, protein, default="chignolin"
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
            "nonbonded_cutoff", nonbonded_cutoff, default=0.9 * unit.nanometer
        )
        self.switch_distance = self.system_parameter(
            "switch_distance", switch_distance, default=None
        )

        # create system
        if solvated:
            raise NotImplementedError("Solvated structures have not been added, yet.")
            #nonbonded_method = app.PME
            #filename = f"{protein}_solvated.pdb"
            #forcefield = [fname for fname in forcefield if fname != "amber99_obc.xml"]
        else:
            nonbonded_method = app.NoCutoff
            pdb_filename = f"{FASTFOLDER_NAMES[protein]}.pdb"
            psf_filename = f"{FASTFOLDER_NAMES[protein]}.psf"
            forcefield = [fname for fname in forcefield if fname != "tip3p.xml"]

        # download file
        full_pdb_filename = os.path.join(root, pdb_filename)
        full_psf_filename = os.path.join(root, psf_filename)
        if download:
            self._download(pdb_filename, root=root)
            self._download(psf_filename, root=root)
        else:
            assert os.path.isfile(full_pdb_filename)
            assert os.path.isfile(full_psf_filename)

        pdb = app.PDBFile(full_pdb_filename)
        psf = app.CharmmPsfFile(full_psf_filename)
        ff = app.ForceField(*forcefield)
        self._system = ff.createSystem(
            psf.topology,
            removeCMMotion=True,
            nonbondedMethod=nonbonded_method,
            nonbondedCutoff=nonbonded_cutoff,
            switchDistance=switch_distance,
            constraints=constraints,
            hydrogenMass=hydrogen_mass,
            rigidWater=True
        )
        self._positions = pdb.positions
        self._topology = psf.topology

    def _download(self, filename, root):
        # get checksum
        md5 = None
        with open(get_data_file("../data/md5sums.txt"), "r") as fp:
            for line in fp:
                if line.split()[1] == filename:
                    md5 = line.split()[0]
        assert md5 is not None
        # download
        download_url(self.url + filename, root, filename, md5)
        assert os.path.isfile(os.path.join(root, filename))


FASTFOLDER_NAMES = {
    "bba": "1FME",
    "villin": "2F4K",
    "trpcage": "2JOF",
    "bbl": "2WAV",
    "alpha3d": "A3D",
    "chignolin": "CLN025",
    "wwdomain": "GTT",
    "lambdarepressor": "lambda",
    "ntl9": "NTL9",
    "proteinG": "NuG2",
    "proteinB": "PRB",
    "homeodomain": "UVF",
}