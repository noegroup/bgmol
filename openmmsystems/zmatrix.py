
import os
import warnings

import numpy as np
import yaml
from .util import get_data_file


__all__ = ["make_protein_z_matrix", "ZMatrixFactory"]


basis_Zs = {}


basis_Zs['ALA'] = [
    ["H", "N", "CA", "C"],
    ["O", "C", "CA", "N"],
    ["HA", "CA", "N", "C"],
    ["CB", "CA", "N", "C"], # ... improper
    ["HB1", "CB", "CA", "N"],
    ["HB2", "CB", "CA", "HB1"],
    ["HB3", "CB", "CA", "HB2"]
    ]#

basis_Zs['LEU'] = [
    ["H", "N", "CA", "C"], # improper ...
    ["O", "C", "CA", "N"],
    ["HA", "CA", "N", "C"],
    ["CB", "CA", "N", "C"], # ... improper
    ["CG", "CB", "CA", "N"], # torsion
    ["CD1", "CG", "CB", "CA"], # torsion
    ["CD2", "CG", "CD1", "CB"], # improper
    ["HB2", "CB", "CA", "CG"], # improper
    ["HB3", "CB", "CA", "CG"], # improper
    ["HG", "CG", "CD1", "CD2"], # improper
    ["HD11", "CD1", "CG", "CB"], # torsion methyl rotation
    ["HD21", "CD2", "CG", "CB"], # torsion methyl rotation
    ["HD12", "CD1", "CG", "HD11"], # improper
    ["HD13", "CD1", "CG", "HD12"], # improper
    ["HD22", "CD2", "CG", "HD21"], # improper
    ["HD23", "CD2", "CG", "HD22"], # improper
    ]

basis_Zs['ILE'] = [
    ["H", "N", "CA", "C"],
    ["O", "C", "CA", "N"],
    ["HA", "CA", "N", "C"],
    ["CB", "CA", "N", "C"], # ... improper
    ["CG1", "CB", "CA", "N"], # torsion
    ["CG2", "CB", "CA", "CG1"], # improper
    ["CD1", "CG1", "CB", "CA"], # torsion
    ["HB", "CB", "CA", "CG1"], # improper
    ["HG12", "CG1", "CB", "CD1"], # improper
    ["HG13", "CG1", "CB", "CD1"], # improper
    ["HD11", "CD1", "CG1", "CB"], # torsion methyl rotation
    ["HD12", "CD1", "CG1", "HD11"], # improper
    ["HD13", "CD1", "CG1", "HD11"], # improper
    ["HG21", "CG2", "CB", "CA"], # methyl torsion
    ["HG22", "CG2", "CB", "HG21"], # improper
    ["HG23", "CG2", "CB", "HG21"], # improper
    ]

basis_Zs['CYS'] = [
    ["H", "N", "CA", "C"],
    ["O", "C", "CA", "N"],
    ["HA", "CA", "N", "C"],
    ["CB", "CA", "N", "C"], # ... improper
    ["SG", "CB", "CA", "N"], # torsion
    ["HB2", "CB", "CA", "SG"], # improper
    ["HB3", "CB", "CA", "SG"], # improper
    ]

basis_Zs['HIS'] = [
    ["H", "N", "CA", "C"],
    ["O", "C", "CA", "N"],
    ["HA", "CA", "N", "C"],
    ]

basis_Zs['ASP'] = [
    ["H", "N", "CA", "C"],
    ["O", "C", "CA", "N"],
    ["HA", "CA", "N", "C"],
    ["CB", "CA", "N", "C"], # ... improper
    ["CG", "CB", "CA", "N"], # torsion
    ["OD1", "CG", "CB", "CA"], # torsion
    ["OD2", "CG", "CB", "OD1"], # torsion
    ["HB2", "CB", "CA", "CG"], # improper
    ["HB3", "CB", "CA", "CG"], # improper
    ]

basis_Zs['ASN'] = [
    ["H", "N", "CA", "C"],
    ["O", "C", "CA", "N"],
    ["HA", "CA", "N", "C"],
    ["CB", "CA", "N", "C"], # ... improper
    ["CG", "CB", "CA", "N"], # torsion
    ["OD1", "CG", "CB", "CA"], # torsion
    ["ND2", "CG", "CB", "OD1"], # improper
    ["HB2", "CB", "CA", "CG"], # improper
    ["HB3", "CB", "CA", "CG"], # improper
    ["HD21", "ND2", "CG", "CB"], # torsion NH2
    ["HD22", "ND2", "CG", "HD21"], # improper
    ]

basis_Zs['GLN'] = [
    ["H", "N", "CA", "C"],
    ["O", "C", "CA", "N"],
    ["HA", "CA", "N", "C"],
    ["CB", "CA", "N", "C"], # ... improper
    ["CG", "CB", "CA", "N"], # torsion
    ["CD", "CG", "CB", "CA"], # torsion
    ["OE1", "CD", "CG", "CB"], # torsion
    ["NE2", "CD", "CG", "OE1"], # torsion
    ["HB2", "CB", "CA", "CG"], # improper
    ["HB3", "CB", "CA", "CG"], # improper
    ["HG2", "CG", "CB", "CD"], # improper
    ["HG3", "CG", "CB", "CD"], # improper
    ["HE21", "NE2", "CD", "CG"], # NH2 torsion
    ["HE22", "NE2", "CD", "HE21"], # improper
    ]

basis_Zs['GLU'] = [
    ["H", "N", "CA", "C"],
    ["O", "C", "CA", "N"],
    ["HA", "CA", "N", "C"],
    ["CB", "CA", "N", "C"], # ... improper
    ["CG", "CB", "CA", "N"], # torsion
    ["CD", "CG", "CB", "CA"], # torsion
    ["OE1", "CD", "CG", "CB"], # torsion
    ["OE2", "CD", "CG", "OE1"], # imppoper
    ["HB2", "CB", "CA", "CG"], # improper
    ["HB3", "CB", "CA", "CG"], # improper
    ["HG2", "CB", "CG", "CD"], # improper
    ["HG3", "CB", "CG", "CD"], # improper
    ]

basis_Zs['GLY'] = [
    ["H", "N", "CA", "C"],
    ["O", "C", "CA", "N"],
    ["HA2", "CA", "N", "C"],
    ["HA3", "CA", "C", "N"]
    ]

basis_Zs['TRP'] = [
    ["H", "N", "CA", "C"],
    ["O", "C", "CA", "N"],
    ["HA", "CA", "N", "C"],
    ]

basis_Zs['TYR'] = [
    ["H", "N", "CA", "C"],
    ["O", "C", "CA", "N"],
    ["HA", "CA", "N", "C"],
    ["CB", "CA", "N", "C"], # ... improper
    ["CG", "CB", "CA", "N"], # torsion
    ["CD1", "CG", "CB", "CA"], # torsion
    ["CD2", "CG", "CB", "CD1"], # improper
    ["CE1", "CD1", "CG", "CB"], # torsion but rigid
    ["CE2", "CD2", "CG", "CD1"], # torsion but rigid
    ["CZ", "CE1", "CE2", "CD1"], # improper
    ["OH", "CZ", "CE1", "CE2"], # improper
    ["HB2", "CB", "CA", "CG"], # improper
    ["HB3", "CB", "CA", "CG"], # improper
    ["HD1", "CD1", "CG", "CE1"], # improper
    ["HD2", "CD2", "CG", "CE2"], # improper
    ["HE1", "CE1", "CD1", "CZ"], # improper
    ["HE2", "CE2", "CD2", "CZ"], # improper
    ["HH", "OH", "CZ", "CE1"], # improper
    ]

basis_Zs['SER'] = [
    ["H", "N", "CA", "C"],
    ["O", "C", "CA", "N"],
    ["HA", "CA", "N", "C"],
    ["CB", "CA", "N", "C"], # ... improper
    ["OG", "CB", "CA", "N"], # torsion
    ["HB2", "CB", "CA", "OG"], # improper
    ["HB3", "CB", "CA", "OG"], # improper
    ["HG", "OG", "CB", "CA"], # torsion
    ]

basis_Zs['PRO'] = [
    ["O", "C", "CA", "N"],
    ["HA", "CA", "N", "C"],
    ["CB", "CA", "N", "C"], # ... improper
    ["CG", "CB", "CA", "N"], # torsion
    ["CD", "CG", "CB", "CA"], # torsion
    ["HB2", "CB", "CA", "CG"], # improper
    ["HB3", "CB", "CA", "CG"], # improper
    ["HG2", "CG", "CB", "CD"], # improper
    ["HG3", "CG", "CB", "CD"], # improper
    ["HD2", "CD", "CG", "N"], # improper
    ["HD3", "CD", "CG", "N"], # improper
    ]

basis_Zs['ARG'] = [
    ["H", "N", "CA", "C"],
    ["O", "C", "CA", "N"],
    ["HA", "CA", "N", "C"],
    ["CB", "CA", "N", "C"], # ... improper
    ["CG", "CB", "CA", "N"], # torsion
    ["CD", "CG", "CB", "CA"], # torsion
    ["NE", "CD", "CG", "CB"], # torsion
    ["CZ", "NE", "CD", "CG"], # torsion
    ["NH1", "CZ", "NE", "CD"], # improper
    ["NH2", "CZ", "NE", "NH1"], # improper
    ["HB2", "CB", "CA", "CG"], # improper
    ["HB3", "CB", "CA", "CG"], # improper
    ["HG2", "CG", "CB", "CD"], # improper
    ["HG3", "CG", "CB", "CD"], # improper
    ["HD2", "CD", "CG", "NE"], # improper
    ["HD3", "CD", "CG", "NE"], # improper
    ["HE", "NE", "CD", "CZ"], # improper
    ["HH11", "NH1", "CZ", "NE"], # improper
    ["HH12", "NH1", "CZ", "HH11"], # improper
    ["HH21", "NH2", "CZ", "NE"], # improper
    ["HH22", "NH2", "CZ", "HH21"], # improper
    ]

basis_Zs['LYS'] = [
    ["H", "N", "CA", "C"],
    ["O", "C", "CA", "N"],
    ["HA", "CA", "N", "C"],
    ["CB", "CA", "N", "C"], # ... improper
    ["CG", "CB", "CA", "N"], # torsion
    ["CD", "CG", "CB", "CA"], # torsion
    ["CE", "CD", "CG", "CB"], # torsion
    ["NZ", "CE", "CD", "CG"], # torsion
    ["HB2", "CB", "CA", "CG"], # improper
    ["HB3", "CB", "CA", "CG"], # improper
    ["HG2", "CG", "CB", "CD"], # improper
    ["HG3", "CG", "CB", "CD"], # improper
    ["HD2", "CD", "CG", "CE"], # improper
    ["HD3", "CD", "CG", "CE"], # improper
    ["HE2", "CE", "CD", "NZ"], # improper
    ["HE3", "CE", "CD", "NZ"], # improper
    ["HZ1", "NZ", "CE", "CD"], # NH3 torsion
    ["HZ2", "NZ", "CE", "HZ1"], # improper
    ["HZ3", "NZ", "CE", "HZ2"], # improper
    ]

basis_Zs['MET'] = [
    ["H", "N", "CA", "C"],
    ["O", "C", "CA", "N"],
    ["HA", "CA", "N", "C"],
    ["CB", "CA", "N", "C"], # ... improper
    ["CG", "CB", "CA", "N"], # torsion
    ["SD", "CG", "CB", "CA"], # torsion
    ["CE", "SD", "CG", "CB"], # torsion
    ["HB2", "CB", "CA", "CG"], # improper
    ["HB3", "CB", "CA", "CG"], # improper
    ["HG2", "CG", "CB", "SD"], # improper
    ["HG3", "CG", "CB", "SD"], # improper
    ["HE1", "CE", "SD", "CG"], # torsion
    ["HE2", "CE", "SD", "HE1"], # improper
    ["HE3", "CE", "SD", "HE2"], # improper
    ]

basis_Zs['THR'] = [
    ["H", "N", "CA", "C"],
    ["O", "C", "CA", "N"],
    ["HA", "CA", "N", "C"],
    ["CB", "CA", "N", "C"], # ... improper
    ["OG1", "CB", "CA", "N"], # torsion
    ["CG2", "CB", "CA", "OG1"], # improper
    ["HB", "CB", "CA", "OG1"], # improper
    ["HG1", "OG1", "CB", "CA"], # torsion
    ["HG21", "CG2", "CB", "CA"], # methyl torsion
    ["HG22", "CG2", "CB", "HG21"], # improper
    ["HG23", "CG2", "CB", "HG21"], # improper
    ]

basis_Zs['VAL'] = [
    ["H", "N", "CA", "C"],
    ["O", "C", "CA", "N"],
    ["HA", "CA", "N", "C"],
    ["CB", "CA", "N", "C"], # ... improper
    ["CG1", "CB", "CA", "N"], # torsion
    ["CG2", "CB", "CA", "CG1"], # improper
    ["HB", "CB", "CA", "CG1"], # improper
    ["HG11", "CG1", "CB", "CA"], # methyl torsion
    ["HG12", "CG1", "CB", "HG11"], # improper
    ["HG13", "CG1", "CB", "HG11"], # improper
    ["HG21", "CG2", "CB", "CA"], # methyl torsion
    ["HG22", "CG2", "CB", "HG21"], # improper
    ["HG23", "CG2", "CB", "HG21"], # improper
    ]

basis_Zs['PHE'] = [
    ["H", "N", "CA", "C"],
    ["O", "C", "CA", "N"],
    ["HA", "CA", "N", "C"],
    ["CB", "CA", "N", "C"], # ... improper
    ["CG", "CB", "CA", "N"], # torsion
    ["CD1", "CG", "CB", "CA"], # torsion
    ["CD2", "CG", "CB", "CD1"], # improper
    ["CE1", "CD1", "CG", "CB"], # torsion but rigid
    ["CE2", "CD2", "CG", "CD1"], # torsion but rigid
    ["CZ", "CE1", "CE2", "CD1"], # improper
    ["HB2", "CB", "CA", "CG"], # improper
    ["HB3", "CB", "CA", "CG"], # improper
    ["HD1", "CD1", "CG", "CE1"], # improper
    ["HD2", "CD2", "CG", "CE2"], # improper
    ["HE1", "CE1", "CD1", "CZ"], # improper
    ["HE2", "CE2", "CD2", "CZ"], # improper
    ["HZ", "CZ", "CE1", "CE2"], # improper
    ]


def make_protein_z_matrix(mdtraj_topology, cartesian="backbone"):
    """
        mdtraj_topology
        cartesian: if not None MDTraj selection string of atoms not to represent with internal coordinates
    """
    z_matrix = []

    not_ic = []
    if cartesian is not None:
        not_ic = mdtraj_topology.select(cartesian)

    for i, residue in enumerate(mdtraj_topology.residues):
        is_nterm = i == 0
        is_cterm = (i + 1) == mdtraj_topology.n_residues

        resatoms = {a.name: a.index for a in residue.atoms}
        resname = residue.name
        for entry in basis_Zs[resname]:  # template entry:
            if resatoms[entry[0]] not in not_ic:
                z_matrix.append([resatoms[_e] for _e in entry])

        if is_nterm:
            # set two additional N-term protons
            if resatoms["H2"] not in not_ic:
                z_matrix.append([resatoms["H2"], resatoms["N"], resatoms["CA"], resatoms["H"] ])
            if resatoms["H3"] not in not_ic:
                z_matrix.append([resatoms["H3"], resatoms["N"], resatoms["CA"], resatoms["H2"] ])
        elif is_cterm:
            # place OXT
            if resatoms["OXT"] not in not_ic:
                z_matrix.append([resatoms["OXT"], resatoms["C"], resatoms["CA"], resatoms["O"]])

    if cartesian is not None:
        return np.array(z_matrix), np.array(not_ic)
    else:
        return np.array(z_matrix)


class ZMatrixFactory:
    TEMPLATE_LOOKUP_DIR = get_data_file("../data/")

    @property
    def z_matrix(self):
        return np.array(self._z)

    @property
    def fixed(self):
        return np.array(self._cartesian)

    def __init__(self, mdtraj_topology, cartesian=()):
        self.top = mdtraj_topology
        self._cartesian = self._select(cartesian)
        self.graph = self.top.to_bondgraph()
        self._atoms = sorted(list(self.graph.nodes), key=lambda node: node.index)
        self._distances = self._build_distance_matrix(self.graph)
        self._z = []

    # === naive builder and helpers ===

    def build_naive(self, subset="all"):
        """Place atoms relative to the closest atoms that are already placed (wrt. bond topology)."""
        subset = self._select(subset)
        if len(list(self._placed_atoms())) == 0:
            self._z = [[0, -1, -1, -1]]
        for atom in self._placed_atoms():
            for neighbor in self._neighbors(atom):
                if neighbor in subset:
                    if not self._is_placed(neighbor):
                        closest = self._3closest_placed_atoms(neighbor)
                        # pad with -1
                        closest3 = np.pad(closest, (0, 3-len(closest)), constant_values=-1)
                        self._z.append([neighbor, *closest3])
        return self.z_matrix, self.fixed

    @staticmethod
    def _build_distance_matrix(graph):
        import networkx as nx
        distances = nx.all_pairs_shortest_path_length(graph)
        matrix = np.zeros((len(graph), len(graph)))
        for this, dist in distances:
            for other in dist:
                matrix[this.index, other.index] = dist[other]
        return matrix

    def _neighbors(self, i):
        for neighbor in self.graph.neighbors(self._atoms[i]):
            yield neighbor.index

    def _is_placed(self, i):
        return any(torsion[0] == i for torsion in self._z)

    def _placed_atoms(self):
        for f in self._cartesian:
            yield f
        for torsion in self._z:
            yield torsion[0]

    def _torsion_index(self, i):
        return list(self._placed_atoms()).index(i)

    def _3closest_placed_atoms(self, i):
        placed = np.array(list(self._placed_atoms()))
        argsort = np.argsort(self._distances[i, placed])
        closest_atoms = placed[argsort[:3]]
        return closest_atoms

    # === template builder and helpers ===

    def build_with_templates(self, yaml_file="z_protein.yaml", *yaml_files, build_protein_backbone=True, subset="all"):
        subset = self._select(subset)
        templates = self._load_templates(yaml_file, *yaml_files)

        # build_backbone
        if build_protein_backbone:
            self.build_naive(subset=np.intersect1d(self.top.select("backbone"), subset))

        # build residues
        for i, residue in enumerate(self.top.residues):
            is_nterm = (i == 0) and residue.is_protein
            is_cterm = ((i + 1) == self.top.n_residues) and residue.is_protein

            resatoms = {a.name: a.index for a in residue.atoms}
            resname = residue.name
            for entry in templates[resname]:  # template entry:
                if not self._is_placed(resatoms[entry[0]]):  # not in not_ic:
                    self._z.append([resatoms[_e] for _e in entry])

            if is_nterm:
                # set two additional N-term protons
                if not self._is_placed(resatoms["H2"]):  # not in not_ic:
                    self._z.append([resatoms["H2"], resatoms["N"], resatoms["CA"], resatoms["H"]])
                if not self._is_placed(resatoms["H3"]):  # not in not_ic:
                    self._z.append([resatoms["H3"], resatoms["N"], resatoms["CA"], resatoms["H2"]])
            elif is_cterm:
                # place OXT
                if not self._is_placed(resatoms["OXT"]):  # not in not_ic:
                    self._z.append([resatoms["OXT"], resatoms["C"], resatoms["CA"], resatoms["O"]])

        # append missing
        placed = np.array(list(self._placed_atoms()))
        if not len(placed) == self.top.n_atoms:
            warnings.warn(
                f"Not all atoms found in templates. Applying naive reconstruction for missing atoms: "
                f"{np.setdiff1d(np.arange(self.top.n_atoms), placed)}"
            )
            self.build_naive(subset)
        return self.z_matrix, self.fixed

    def _load_template(self, yaml_file):
        filename_in_pkg = os.path.join(self.TEMPLATE_LOOKUP_DIR, yaml_file)
        if os.path.isfile(filename_in_pkg):
            if os.path.isfile(yaml_file) and os.path.normpath(yaml_file) != os.path.normpath(filename_in_pkg):
                raise warnings.warn(
                    f"{yaml_file} exists locally and in the package templates. "
                    f"Taking the built-in one from the openmmsystems package.", UserWarning
                )
            yaml_file = filename_in_pkg
        with open(yaml_file, "r") as f:
            templates = yaml.load(f, yaml.SafeLoader)
            for residue in templates:
                if not isinstance(templates[residue], list):
                    raise IOError("File format is not acceptable.")
                for torsion in templates[residue]:
                    if not len(torsion) == 4:
                        raise IOError(f"Torsion {torsion} in residue {residue} does not have 4 atoms.")
                    if not all(isinstance(name, str) for name in torsion):
                        raise IOError(f"Torsion {torsion} in residue {residue} does not consist of atom names.")
            return templates

    def _load_templates(self, *yaml_files):
        templates = dict()
        for f in yaml_files:
            template = self._load_template(f)
            for key in template:
                if key in templates:
                    warnings.warn(f"Residue {key} found in multiple files. Updating with the definition from {f}.")
            templates.update(template)
        return templates

    def build_with_system(self, system):
        raise NotImplementedError()

    def _select(self, selection):
        return self.top.select(selection) if isinstance(selection, str) else selection



