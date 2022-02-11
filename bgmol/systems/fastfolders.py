

import os
import tempfile
from bgmol.util.importing import import_openmm
_, unit, app = import_openmm()
from bgmol.systems import OpenMMSystem
from ..tpl.download import download_url
from ..util.pdbpatch import fixed_atom_names

__all__ = ["FAST_FOLDER_NAMES", "FastFolder"]


class FastFolder(OpenMMSystem):
    """Small, fast-folding proteins.

    This is a remake of the proteins from reference [1] in charmm-gui.
    Created by Andreas Krämer.

    Attributes
    ----------
    protein : str
        The name of the protein. See the list `FAST_FOLDER_NAMES` in this module for valid options.
    forcefield : "charmm36m" or list[str], optional
        Either use the charmm36m force field or a list of builtin openmm force field files.
        For openmm force field files not all proteins might be supported
        (e.g. ["amber99sbildn.xml"] does not support villin and ntl9 because of nonstandard residues,
        a norleucine (NLE) and a ALA-CONH2 C-terminus, respectively).
    implicit_solvent : app.internal.singleton.Singleton or None, optional
        Implicit solvent model to be used.
    constraints : app.internal.singleton.Singleton, optional
        Constraint types
    solvated : bool, optional
        Whether to add explicit water and ions.
    hydrogen_mass : unit.Quantity or None, optional
        If None, don't repartition hydrogen mass. Else, assign the specified mass to hydrogen atoms.
    nonbonded_cutoff : unit.Quantity, optional
        The cutoff for nonbonded forces.
    switch_distance : unit.Quantity or None, optional
        Switch distance for nonbonded (LJ) interactions. If None, don't use a switch distance.
    root : str, optional
        The root directory to which to download the files.
    download : bool, optional
        Whether files should be downloaded.

    References
    ----------
    [1] K. Lindorff-Larsen, S. Piana, R.O. Dror, D.E. Shaw,
        How fast-folding proteins fold,
        Science (80-. ). 334 (2011) 517–520.
        doi:10.1126/science.1208351.

    Notes
    -----
    - Requires an internet connection to download the topology and initial structure.
    - Initial solvated structures are not equilibrated
    """
    url = "http://ftp.mi.fu-berlin.de/pub/cmb-data/bgmol/systems/fastfolders/"

    def __init__(
            self,
            protein="chignolin",
            forcefield="charmm36m",
            implicit_solvent=app.OBC2,
            constraints=app.HBonds,
            solvated=False,
            hydrogen_mass=None,
            nonbonded_cutoff=0.9 * unit.nanometer,
            switch_distance=None,
            root=tempfile.gettempdir(),
            download=True
    ):
        super().__init__()

        if not protein in FAST_FOLDER_NAMES:
            raise ValueError(f"protein has to be one of {FAST_FOLDER_NAMES}")

        # register parameters
        self.protein = self.system_parameter(
            protein, protein, default="chignolin"
        )
        self.forcefield = self.system_parameter(
            "forcefield", forcefield, default="charmm36m"
        )
        self.implicit_solvent = self.system_parameter(
            "implicit_solvent", implicit_solvent, default=app.OBC2
        )
        self.constraints = self.system_parameter(
            "constraints", constraints, default=app.HBonds
        )
        self.solvated = self.system_parameter(
            "solvated", solvated, default=False
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

        # download files
        nonbonded_method = app.PME if self.solvated else app.NoCutoff
        if download:
            self._download(root=root)
        local_files = {key: os.path.join(root, value) for key, value in self.source_files.items()}
        for f in local_files.values():
            assert os.path.isfile(f)

        # Load the CHARMM files
        # (keep the CONH2 term of nlt9)
        with fixed_atom_names(ALA=["HT1", "HT2"]):
            psf = app.CharmmPsfFile(local_files.pop("psf"))
            pdb = app.PDBFile(local_files.pop("pdb"))

        # specify solvation properties
        if solvated:
            a = _INIT_BOX_SIZES[self.protein] * unit.nanometers
            psf.setBox(a, a, a, 90., 90., 90.)
            implicit_solvent = None
        else:
            implicit_solvent = self.implicit_solvent

        # create system
        create_system_kwargs = dict(
            removeCMMotion=True,
            nonbondedMethod=nonbonded_method,
            nonbondedCutoff=nonbonded_cutoff,
            switchDistance=switch_distance,
            constraints=constraints,
            hydrogenMass=hydrogen_mass,
            rigidWater=True,
            implicitSolvent=implicit_solvent
        )
        if self.forcefield == "charmm36m":
            params = app.CharmmParameterSet(
                *list(local_files.values())
            )
            self._system = psf.createSystem(params, **create_system_kwargs)
        elif isinstance(self.forcefield, tuple) or isinstance(self.forcefield, list):
            ff = app.ForceField(*self.forcefield)
            self._system = ff.createSystem(psf.topology, **create_system_kwargs)
        else:
            raise ValueError(f"forcefield {self.forcefield} is not supported for FastFolder")

        self._positions = pdb.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        self._topology = psf.topology

    def _download(self, root):
        md5list = _MD5SUMS.strip().split("\n")
        md5dict = {line.split()[1]: line.split()[0] for line in md5list}
        if not os.path.isdir(os.path.join(root, "FastFolder")):
            os.mkdir(os.path.join(root, "FastFolder"))
        for filename in self.source_files.values():
            # get checksum
            md5 = md5dict[filename]
            download_url(self.url + filename, root, filename, md5)

    @property
    def source_files(self):
        solvated_suffix = "_solvated" if self.solvated else ""
        files = {
            "pdb": f"FastFolder/{self.protein}{solvated_suffix}.pdb",
            "psf": f"FastFolder/{self.protein}{solvated_suffix}.psf"
        }
        if self.forcefield == "charmm36m":
            files["rtf"] = "FastFolder/top_all36_prot.rtf"
            files["prm"] = "FastFolder/par_all36m_prot.prm"
            files["daa"] = "FastFolder/toppar_all36_prot_c36m_d_aminoacids.str"
            files["str"] = "FastFolder/toppar_water_ions.str"
        return files

    def reinitialize_energy_model(self, temperature=300, **kwargs):
        """override default n_workers to 1"""
        n_workers = kwargs.pop("n_workers", 1)
        super().reinitialize_energy_model(temperature=temperature, n_workers=n_workers, **kwargs)


_INIT_BOX_SIZES = {
    "bba": 4.7,
    "lambdarepressor": 7.0,
    "proteinG": 5.5,
    "proteinB": 5.6,
    "alpha3d": 6.4,
    "villin": 5.4,
    "ntl9": 5.0,
    "trpcage": 3.7,
    "homeodomain": 5.3,
    "bbl": 5.7,
    "chignolin": 4.0,
    "wwdomain": 5.2,
}


FAST_FOLDER_NAMES = tuple(_INIT_BOX_SIZES.keys())


_MD5SUMS = """
67560228af0b5a5d352229de77554b80  FastFolder/alpha3d.crd
28e4ae3c3b172a051ff5eb30fc700d29  FastFolder/alpha3d.pdb
c3d9f19d9e5d25ab4c2dd207fd853dde  FastFolder/alpha3d.psf
7185f07b50944013d30e2e0447b435b3  FastFolder/alpha3d_solvated.crd
c25ff5227fd3c14683a4d17bb847f9fd  FastFolder/alpha3d_solvated.pdb
c98ce084c3d2c5c2f5d5c287e1e84131  FastFolder/alpha3d_solvated.psf
e88d1433b9925b857e3e5bb34036081f  FastFolder/bba.crd
621b92c09318082ca819cdb17a5ed0cd  FastFolder/bba.pdb
50986e441d304afe7a86cbed0a1c5b6d  FastFolder/bba.psf
37520eadebba18e2695e708941a4f414  FastFolder/bba_solvated.crd
f782409069244a9fc5869657f0ba9022  FastFolder/bba_solvated.pdb
1fe729bab232748af888d6a5cb919b57  FastFolder/bba_solvated.psf
37b5d677b1f30cb52174989c0639363e  FastFolder/bbl.crd
f78485777df3122e6f3f70b8aa1d27c1  FastFolder/bbl.pdb
018bfeb1949056def5a92042e0251750  FastFolder/bbl.psf
3e7fd42c004421bb85cf85f463d69ebb  FastFolder/bbl_solvated.crd
d55bf739c7a62a3aa5ff401987b573c0  FastFolder/bbl_solvated.pdb
4e12cacc2fe786b635b34fbdf5653df1  FastFolder/bbl_solvated.psf
1cc98afb3f016e66b032fc968144ccb5  FastFolder/chignolin.crd
353430297cdd6416e7445ddd90609809  FastFolder/chignolin.pdb
b0128e863e73c13e3dca02fd0727bb03  FastFolder/chignolin.psf
0a8fe81d078198819e223e71a49783ac  FastFolder/chignolin_solvated.crd
ac2a1055db280b4cfd0ed6232dc7c78b  FastFolder/chignolin_solvated.pdb
23a06e7e0b0dcbcb65c827f4f9b92aaf  FastFolder/chignolin_solvated.psf
68333d7289cabf75f8d59854d7d62c95  FastFolder/homeodomain.crd
42b9aa9d9f7045b0077d7e83772996e7  FastFolder/homeodomain.pdb
b8761ecf1140cc666ed828f8f8197ff7  FastFolder/homeodomain.psf
d21e6e393db71c755a740fd5462855c4  FastFolder/homeodomain_solvated.crd
bb7442699fd822828b42b49a35830981  FastFolder/homeodomain_solvated.pdb
9cd6983bc1fec833e529fb617acc41da  FastFolder/homeodomain_solvated.psf
a2b782664a3e3513de41f43541f24d1b  FastFolder/lambdarepressor.crd
fc874cc443e54e27566a95b032b43e0e  FastFolder/lambdarepressor.pdb
3b190eb6d5be589ccb7b03318fc7105a  FastFolder/lambdarepressor.psf
d24883eaf8f875f668c8a271394a4ca3  FastFolder/lambdarepressor_solvated.crd
f19bdfa199828bd4101e93316178141f  FastFolder/lambdarepressor_solvated.pdb
d4c684e13166793b141a620e58294164  FastFolder/lambdarepressor_solvated.psf
b31e0856849ec40fc188b3de4b32c1b2  FastFolder/ntl9.crd
5fe45ec6d01fef41c3b99e347bf5b609  FastFolder/ntl9.pdb
728257ce6df22a1ea8f748a774eb7429  FastFolder/ntl9.psf
51bd66c5a8264272178314f402395168  FastFolder/ntl9_solvated.crd
b0bfbd975c6bbcb9e477ed59b28ee33e  FastFolder/ntl9_solvated.pdb
b3118d4def6b7809ae7b254a0860ae05  FastFolder/ntl9_solvated.psf
74edc2fb5d4d0cefd2a712bb04c2da6f  FastFolder/par_all36m_prot.prm
5fd53602d3271c29715628690bd15dcc  FastFolder/proteinB.crd
1559e29257d0c375e4acdd0bf4427520  FastFolder/proteinB.pdb
1c9d169b674ebebee77f51d3e545e8bc  FastFolder/proteinB.psf
fcf3fdbecb4a2d2f592aa74e4a472fbc  FastFolder/proteinB_solvated.crd
c538f7e84d52bea059b3d71799e1f6c4  FastFolder/proteinB_solvated.pdb
c63631302a11b792bdb8518326f9cefa  FastFolder/proteinB_solvated.psf
4f54ef1a0399bbbf117458d879fcbbff  FastFolder/proteinG.crd
58f9ba85c536bf4721d5048e2f33717d  FastFolder/proteinG.pdb
97f663ad2668f4fbe66ddf65472401f4  FastFolder/proteinG.psf
bb82cc3dced4c984b74ef54a42a50d12  FastFolder/proteinG_solvated.crd
7d65cd36c9aefa88abf8da0a9dd6191e  FastFolder/proteinG_solvated.pdb
eeaf53232ea29761c160f262b9387513  FastFolder/proteinG_solvated.psf
34a1946f85db4679b151f598e36e74f4  FastFolder/top_all36_prot.rtf
53e2842122d49e1f3b70e9e4286e61d8  FastFolder/toppar_all36_prot_c36m_d_aminoacids.str
789b30edc593da50bfcf7c7765cdc26b  FastFolder/toppar_water_ions.str
69821232eedc2e20093fc6486eca455d  FastFolder/trpcage.crd
fff8860471ecac54efd9699e2f8f2d04  FastFolder/trpcage.pdb
46cfc37c46d2c1cdc2bc99e84431ff24  FastFolder/trpcage.psf
0dc1abc76c24a9cf7ae23bbb6e96d67f  FastFolder/trpcage_solvated.crd
849e8e0fa565efa205917503753a4e97  FastFolder/trpcage_solvated.pdb
911c2ddc9dc033d38408f3ef555188d1  FastFolder/trpcage_solvated.psf
93ffee76616bca1b95899f94542c86fc  FastFolder/villin.crd
fb4c4eb4032b27d072f143c462322fa1  FastFolder/villin.pdb
a785841ed12fc5beaaf9f1a29eda20ab  FastFolder/villin.psf
b03e832024bb4ce6bfb3041cdf2c0dd1  FastFolder/villin_solvated.crd
f010e0a9f6d7183630c31cd7d439bc64  FastFolder/villin_solvated.pdb
d2a7fc02cb992a9eaaf58b8c3a9e5b9f  FastFolder/villin_solvated.psf
87f90218414ae226c26337f7e394bd7a  FastFolder/wwdomain.crd
f3ce7eb0ff97645c7be86f026a9dbb5f  FastFolder/wwdomain.pdb
035bfbdd851f6c77a757c00cb9441379  FastFolder/wwdomain.psf
8a6eeb26a1eb4d20b193f58f2c3f0e15  FastFolder/wwdomain_solvated.crd
fd289a5f7056447fd783acd5830b390e  FastFolder/wwdomain_solvated.pdb
0bb4c746f8269c8bc02b1c451eca5f28  FastFolder/wwdomain_solvated.psf
"""
