
import os
import warnings

import numpy as np
import yaml
from .util import get_data_file


__all__ = ["ZMatrixFactory"]


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

    def _select(self, selection):
        return self.top.select(selection) if isinstance(selection, str) else selection

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




