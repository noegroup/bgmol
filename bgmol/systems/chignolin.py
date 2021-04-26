import os
import warnings
import tempfile
from simtk import unit
from simtk.openmm import app
from ..systems.base import OpenMMSystem
from ..tpl.download import download_url

__all__ = ["ChignolinC22Implicit"]


class ChignolinC22Implicit(OpenMMSystem):
    """Chignolin miniprotein with CHARMM22* force field in implicit solvent.

    Attributes
    ----------
    constraints : app.internal.singleton.Singleton
        Constraint types
    hydrogen_mass : unit.Quantity or None
        If None, don't repartition hydrogen mass. Else, assign the specified mass to hydrogen atoms.
    implicit_solvent : app.internal.singleton.Singleton or None
        Implicit solvent model to be used.
    root : str
        The root directory to which to download the files.
    download : bool
        Whether files should be downloaded.

    Notes
    -----
    Requires an internet connection to download the initial structure.
    """
    url = "ftp://ftp.mi.fu-berlin.de/pub/cmb-data/bgmol/systems/chignolin/chignolin_c22_implicit/"

    def __init__(
            self,
            constraints=app.HBonds,
            hydrogen_mass=4.0*unit.amu,
            implicit_solvent=app.OBC2,
            root=tempfile.gettempdir(),
            download=True
    ):
        super(ChignolinC22Implicit, self).__init__()

        self.constraints = self.system_parameter(
            "constraints", constraints, default=app.HBonds
        )
        self.hydrogen_mass = self.system_parameter(
            "hydrogen_mass", hydrogen_mass, default=4.0*unit.amu
        )
        self.implicit_solvent = self.system_parameter(
            "implicit_solvent", implicit_solvent, default=app.OBC2
        )

        # download the system files
        if download:
            self._download(root)
        for sourcefile in self.FILES:
            assert os.path.isfile(os.path.join(root, sourcefile))

        # Load the CHARMM files
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # I don't want to see warnings about mixin force field
            params = app.CharmmParameterSet(
                os.path.join(root, "top_all22star_prot.rtf"),
                os.path.join(root, "top_water_ions.rtf"),
                os.path.join(root, "parameters_ak.prm")
            )

        # create system
        psf = app.CharmmPsfFile(os.path.join(root, "structure.psf"))
        crds = app.PDBFile(os.path.join(root, "structure.pdb"))
        self._system = psf.createSystem(
            params,
            nonbondedMethod=app.NoCutoff,
            constraints=constraints,
            hydrogenMass=hydrogen_mass,
            implicitSolvent=implicit_solvent
        )
        self._positions = crds.positions
        self._topology = crds.topology

    def _download(self, root):
        # download
        for sourcefile in self.FILES:
            md5 = self.FILES[sourcefile]
            download_url(self.url+sourcefile, root, sourcefile, md5)
            assert os.path.isfile(os.path.join(root, sourcefile))

    FILES = {
        "parameters_ak.prm": "d953daf4925e4146a5fcece875ee4e57",
        "structure.pdb": "be19629a75e0ee4e1cc3c72a9ebc63c6",
        "structure.psf": "944b26edb992c7dbdaa441675b9e42c5",
        "top_all22star_prot.rtf": "d046c9a998369be142a6470fd5bb3de1",
        "top_water_ions.rtf": "ade085f88e869de304c814bf2d0e57fe"
    }
