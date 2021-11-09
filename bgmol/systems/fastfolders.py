

import os
import glob
import tempfile
from bgmol.util.importing import import_openmm
_, unit, app = import_openmm()
from bgmol.systems import OpenMMSystem
from bgmol.util import get_data_file
from torchvision.datasets.utils import download_url, download_and_extract_archive

__all__ = ["FASTFOLDER_NAMES", "FastFolder"]


class FastFolder(OpenMMSystem):
    """Small, fast-folding proteins.

    This is a remake of the proteins from reference [1] in charmm-gui.
    Created by Andreas Krämer.

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
    url = "http://ftp.mi.fu-berlin.de/pub/cmb-data/bgmol/systems/fastfolders/FastFolder/"

    def __init__(
            self,
            protein="chignolin",
            forcefield="charmm36m",
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
        if not forcefield in self.ALLOWED_FORCE_FIELDS:
            raise ValueError(f"forcefield has to be in {self.ALLOWED_FORCE_FIELDS}")

        # register parameters
        self.protein = self.system_parameter(
            protein, protein, default="chignolin"
        )
        self.forcefield = self.system_parameter(
            "forcefield", forcefield, default="charmm36m"
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

        # download files
        full_pdb_filename = os.path.join(root, pdb_filename)
        full_psf_filename = os.path.join(root, psf_filename)
        if download:
            self._download(pdb_filename, psf_filename, root=root)

        assert os.path.isfile(full_pdb_filename)
        assert os.path.isfile(full_psf_filename)
        for sourcefile in self.CHARMM_FILES:
            assert os.path.isfile(os.path.join(root, sourcefile))

        # Load the CHARMM files
        charmm_files = [
            os.path.join(root, "toppar", fname)
            for fname in [
                "top_all36_prot.rtf",
                "par_all36m_prot.prm",
                "toppar_all36_prot_c36m_d_aminoacids.str",
                "toppar_water_ions.str",
            ]
        ]
        params = app.CharmmParameterSet(
            *charmm_files
        )

        # create system
        psf = app.CharmmPsfFile(full_psf_filename)
        pdb = app.PDBFile(full_pdb_filename)
        self._system = psf.createSystem(
            params,
            nonbondedMethod=app.NoCutoff,
            constraints=constraints,
            hydrogenMass=hydrogen_mass,
            implicitSolvent=app.OBC2
        )

        self._positions = pdb.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        self._topology = psf.topology

    def _download(self, *files, root):
        # get checksum
        md5 = None
        for sourcefile in files:
            with open(get_data_file("../data/md5sums.txt"), "r") as fp:
                for line in fp:
                    if line.split()[1] == sourcefile:
                        md5 = line.split()[0]
            assert md5 is not None
            # download
            download_url(self.url + sourcefile, root, sourcefile, md5)
            assert os.path.isfile(os.path.join(root, sourcefile))

        for sourcefile in self.CHARMM_FILES:
            md5 = self.CHARMM_FILES[sourcefile]
            download_and_extract_archive(
                url=self.url+sourcefile,
                download_root=root,
                extract_root=root,
                md5=md5,
                remove_finished=True
            )
            assert os.path.isfile(os.path.join(root, sourcefile))

    CHARMM_FILES = {
        "c36.tgz": "75596f1993ce815bb1629b14361c8d32",
    }

    ALLOWED_FORCE_FIELDS = ["charmm36m"]


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