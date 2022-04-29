
import os
import warnings
from copy import copy

import numpy as np
import mdtraj as md
import yaml
from .util import get_data_file, rewire_chiral_torsions
from .util.topology import _select_ha


__all__ = ["ZMatrixFactory", "build_fake_topology"]


class ZMatrixFactory:
    """Factory for internal coordinate representations.

    Parameters
    ----------
    mdtraj_topology : mdtraj.Topology
        The system's topology.
    cartesian : str or sequence of ints
        The ids of atoms that are not transformed into internal coordinates.
        This can either be a sequence of integers or a DSL selection string for mdtraj.

    Attributes
    ----------
    z_matrix : np.ndarray
        The internal coordinate definition of shape (n_atoms - n_cartesian, 4).
        For each line, the placing of atom `line[0]` is conditioned on the atoms `line[1:3]`.
        Negative entries (-1) in a line mean that the atom is a seed for the global transform.
    fixed : np.ndarray
        One-dimensional array of cartesian atom indices.
    """
    TEMPLATE_LOOKUP_DIR = get_data_file("../data/")

    @property
    def z_matrix(self):
        return np.array(self._z).astype(np.int64)

    @z_matrix.setter
    def z_matrix(self, value):
        if isinstance(value, list):
            self._z = value
        elif isinstance(value, np.ndarray):
            self._z = value.tolist()
        else:
            raise ValueError(f"Can't set z_matrix to {value}")

    @property
    def fixed(self):
        return np.array(self._cartesian).astype(np.int64)

    def __init__(self, mdtraj_topology, cartesian=()):
        self.top = mdtraj_topology
        self._cartesian = self._select(cartesian)
        self.graph = self.top.to_bondgraph()
        self._atoms = sorted(list(self.graph.nodes), key=lambda node: node.index)
        self._distances = self._build_distance_matrix(self.graph)
        self._z = []

    # === naive builder and helpers ===

    def build_naive(self, subset="all", render_independent=True, rewire_chiral=True, verbose=False):
        """Place atoms relative to the closest atoms that are already placed (wrt. bond topology).

        Parameters
        ----------
        subset : str or sequence of int
            A selection string or list of atoms. The z-matrix is only build for the subset.
        render_independent : bool
            Whether to make sure that no two positions depend on the same three other positions.
        rewire_chiral : bool
            Whether to make sure that all HA depend on (CA, N, C).

        Returns
        -------
        z_matrix : np.ndarray
        fixed_atoms : np.ndarray
        """
        subset = self._select(subset)
        current_atoms = set(self._placed_atoms())
        if len(current_atoms) < 3:
            self._z = self._seed_z(current_atoms, subset)
            for torsion in self._z:
                current_atoms.add(torsion[0])

        while any(not self._is_placed(atom) for atom in subset):
            # add neighbors
            next_atoms = copy(current_atoms)
            for atom in current_atoms:
                assert self._is_placed(atom)
                for neighbor in self._neighbors(atom):
                    if not self._is_placed(neighbor) and neighbor in subset:
                        next_atoms.add(neighbor)
                next_atoms.remove(atom)
            # if there weren't any neighbors, add one in the distance matrix
            if len(next_atoms) == 0:
                next_atoms = self._closest_in_bond_graph(list(current_atoms), list(self._remaining_atoms(subset)))

            current_atoms = next_atoms
            # build part of z matrix
            z = []
            for atom in current_atoms:
                closest = self._3closest_placed_atoms(atom, subset=subset)
                if len(closest) == 3:
                    z.append([atom, *closest])
                else:
                    print("!", atom, closest)
            self._z.extend(z)
        if rewire_chiral:
            self.z_matrix = rewire_chiral_torsions(self.z_matrix, self.top, verbose=verbose)
        if render_independent:
            self.render_independent(keep=_select_ha(self.top))
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

    def _seed_z(self, current, subset):
        current = list(current)
        if len(current) == 0:
            current = [subset[0]]
        closest = self._3closest_placed_atoms(current[0], subset, subset)
        seed = np.unique([*current, *closest])[:3]
        z = [
            [seed[0], -1, -1, -1],
            [seed[1], seed[0], -1, -1],
            [seed[2], seed[1], seed[0], -1],
        ]
        return z

    def _neighbors(self, i):
        for neighbor in self.graph.neighbors(self._atoms[i]):
            yield neighbor.index

    def _is_placed(self, i):
        return (
                any(torsion[0] == i for torsion in self._z)
                or any(f == i for f in self._cartesian)
        )

    def _placed_atoms(self):
        for f in self._cartesian:
            yield f
        for torsion in self._z:
            yield torsion[0]

    def _remaining_atoms(self, subset):
        for atom in subset:
            if not self._is_placed(atom):
                yield atom

    def _closest_in_bond_graph(self, placed, subset):
        min_distance = self._distances[placed][:, subset].min()
        is_in_closest = (self._distances == min_distance)[placed]
        is_in_closest = np.any(is_in_closest, axis=0)
        closest,  = np.where(is_in_closest)
        return set(np.intersect1d(closest, subset))

    def _torsion_index(self, i):
        return list(self._placed_atoms()).index(i)

    def _3closest_placed_atoms(self, i, placed=None, subset=None):
        if placed is None:
            placed = np.array(list(self._placed_atoms()))
        if subset is not None:
            placed = np.intersect1d(placed, subset)
        argsort = np.argsort(self._distances[i, placed])
        closest_atoms = placed[argsort[:3]]
        return closest_atoms

    def _select(self, selection):
        return self.top.select(selection) if isinstance(selection, str) else selection

    # === template builder and helpers ===

    def build_with_templates(
            self,
            *yaml_files,
            build_protein_backbone=True,
            subset="all"
    ):
        """Build ICs from template files.

        Parameters
        ----------
        *yaml_files : str
            filenames of any other template files; if non are passed, use the bundled "z_protein.yaml", "z_termini.yaml"
        build_protein_backbone : bool
            Whether to build the protein backbone first
        subset : str or sequence of int
            A selection string or list of atoms. The z-matrix is only build for the subset.

        Notes
        -----
        For the formatting of template files, see data/z_protein.yaml

        Returns
        -------
        z_matrix : np.ndarray
        fixed_atoms : np.ndarray
        """
        if len(yaml_files) == 0:
            yaml_files = ["z_protein.yaml", "z_termini.yaml"]
        subset = self._select(subset)
        templates = self._load_templates(*yaml_files)

        # build_backbone
        if build_protein_backbone:
            self.build_naive(subset=np.intersect1d(self.top.select("backbone and element != O"), subset))
        # build residues
        residues = list(self.top.residues)
        for i, residue in enumerate(residues):
            is_nterm = (i == 0) and residue.is_protein
            is_cterm = ((i + 1) == self.top.n_residues) and residue.is_protein
            resatoms = {a.name: a.index for a in residue.atoms}
            if not is_cterm:
                resatoms_neighbor = {f"+{a.name}": a.index for a in residues[i+1].atoms}
                resatoms.update(resatoms_neighbor)
            if not is_nterm:
                resatoms_neighbor = {f"-{a.name}": a.index for a in residues[i-1].atoms}
                resatoms.update(resatoms_neighbor)

            # add template definitions
            definitions = templates[residue.name]
            if is_nterm and "NTERM" in templates:
                definitions = definitions + templates["NTERM"]
            if is_cterm and "CTERM" in templates:
                definitions = definitions + templates["CTERM"]

            for entry in definitions:  # template entry:
                if any(e not in resatoms for e in entry):
                    continue  # skip torsions with non-matching atom names
                if resatoms[entry[0]] in subset and not self._is_placed(resatoms[entry[0]]):  # not in not_ic:
                    self._z.append([resatoms[_e] for _e in entry])

        # append missing
        placed = np.array(list(self._placed_atoms()))
        if not len(placed) == len(subset):
            missing = np.setdiff1d(np.arange(self.top.n_atoms), placed)
            warnings.warn(
                f"Not all atoms found in templates. Applying naive reconstruction for missing atoms: "
                f"{tuple(self._atoms[m] for m in missing)}"
            )
            self.build_naive(subset)
        return self.z_matrix, self.fixed

    def _load_template(self, yaml_file):
        filename_in_pkg = os.path.join(self.TEMPLATE_LOOKUP_DIR, yaml_file)
        if os.path.isfile(filename_in_pkg):
            if os.path.isfile(yaml_file) and os.path.normpath(yaml_file) != os.path.normpath(filename_in_pkg):
                raise warnings.warn(
                    f"{yaml_file} exists locally and in the package templates. "
                    f"Taking the built-in one from the bgmol package.", UserWarning
                )
            yaml_file = filename_in_pkg
        with open(yaml_file, "r") as f:
            templates = yaml.load(f, yaml.SafeLoader)
            for residue in templates:
                if not isinstance(templates[residue], list):
                    raise IOError("File format is not acceptable.")
                if "GENERAL" in templates:
                    templates[residue] = [*templates["GENERAL"], *templates[residue]]
                for torsion in templates[residue]:
                    if not len(torsion) == 4:
                        raise IOError(f"Torsion {torsion} in residue {residue} does not have 4 atoms.")
                    if not all(isinstance(name, str) for name in torsion):
                        raise IOError(f"Torsion {torsion} in residue {residue} does not consist of atom names.")
            if "GENERAL" in templates:
                del templates["GENERAL"]
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
        """
        Idea: build a lookup table of torsions and impropers; sort them by marginal entropy.
        Before doing the naive lookup, try to insert in the minimum-entropy torsions.
        Maybe: also do something regarding symmetries.
        """
        raise NotImplementedError()

    def render_independent(self, keep=None):
        """Rewire z matrix so that no two positions depend on the same bonded atom and angle/torsion

        Parameters
        ----------
        keep : Sequence[int]
            All atoms, whose placement should not be changed by any means.
            By default, don't rewire CB to be able to control the chirality.
        """
        keep = _select_ha() if keep is None else keep
        keep = self._select(keep)
        all234 = [(torsion[1], set(torsion[2:])) for torsion in self._z]
        for i in range(len(self._z)):
            torsion = self._z[i]
            this234 = all234[i]
            while this234 in all234[:i]:
                previous_index = all234.index(this234)
                previous_torsion = self._z[previous_index]
                if torsion[0] in keep:
                    assert not previous_torsion[0] in keep
                    # swap
                    self._z[i], self._z[previous_index] = previous_torsion, torsion
                    all234[i], all234[previous_index] = all234[previous_index], all234[i]
                    torsion = self._z[i]
                    this234 = all234[i]
                    continue
                # make sure there are no circular dependencies
                assert previous_torsion[0] not in torsion[:3]
                assert torsion[0] not in previous_torsion
                # rewire
                new_torsion = torsion[:3] + previous_torsion[:1]
                self._z[i] = new_torsion
                this234 = set(new_torsion[2:])
                all234[i] = this234
            assert not this234 in all234[:i]
            assert not self._z[i] in self._z[:i]
        if not self.is_independent(self._z):
            warnings.warn("Z-matrix torsions are not fully independent because of a constraint on HA.")
        return self._z

    @staticmethod
    def is_independent(z):
        all234 = [(torsion[1], set(torsion[2:])) for torsion in z]
        for i, other in enumerate(all234):
            if other in all234[:i]:
                return False
        return True


def build_fake_topology(n_atoms, bonds=None, atoms_by_residue=None, coordinates=None):
    """A stupid function to build an MDtraj topology with limited information.

    Returns
    -------
    topology : md.Topology
    trajectory : md.Trajectory or None
        None, if no coordinates were specified.
    """
    topology = md.Topology()

    if bonds is None:   # assume linear molecule
        bonds = np.column_stack([np.arange(n_atoms-1), np.arange(1, n_atoms)])
    if atoms_by_residue is None:   # 1 atom per residue
        atoms = set([atom for bond in bonds for atom in bond])
        atoms_by_residue = [[atom] for atom in atoms]

    chain = topology.add_chain()
    for res in atoms_by_residue:
        residue = topology.add_residue("C", chain)
        for _ in res:
            topology.add_atom("CA", chain, residue)
    atoms = list(topology.atoms)
    for bond in bonds:
        atom1 = atoms[bond[0]]
        atom2 = atoms[bond[1]]
        topology.add_bond(atom1, atom2)

    trajectory = None
    if coordinates is not None:
        coords = coordinates.reshape(-1, n_atoms, 3)
        trajectory = md.Trajectory(
            xyz=coords,
            topology=topology
        )
    return topology, trajectory
