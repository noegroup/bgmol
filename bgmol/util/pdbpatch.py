from contextlib import contextmanager
from .importing import import_openmm
_, _, app = import_openmm()

__all__ = ["fixed_atom_names"]


@contextmanager
def fixed_atom_names(**residues):
    """Suppress replacing of selected atom names when creating PSF or PDB files.

    Examples
    --------
    >>> with fixed_atom_names(ALA=["NT"]):
    >>>     psf = app.CharmmPsfFile("some.psf")
    """
    def modify(original):
        @staticmethod
        def modified_load_tables():
            original()
            for residue in residues:
                for atom in residues[residue]:
                    if atom in app.PDBFile._atomNameReplacements[residue]:
                        del app.PDBFile._atomNameReplacements[residue][atom]

        return modified_load_tables

    app.PDBFile._loadNameReplacementTables, original = (
        modify(app.PDBFile._loadNameReplacementTables),
        app.PDBFile._loadNameReplacementTables
    )
    yield
    app.PDBFile._loadNameReplacementTables = original
