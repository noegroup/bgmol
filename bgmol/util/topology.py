
from typing import Sequence
import numpy as np
import mdtraj as md


__all__ = ["rewire_chiral_torsions", "find_rings", "is_proper_torsion"]


def find_rings(mdtraj_topology: md.Topology):
    """Find rings in a molecule.

    Returns
    -------
    ring_atoms: List[List[int]]
        Atoms constituting each ring
    """
    import networkx as nx
    graph = mdtraj_topology.to_bondgraph()
    cycles = nx.cycle_basis(graph)
    return [[atom.index for atom in ring] for ring in cycles]


def is_proper_torsion(torsions: Sequence[Sequence[int]], mdtraj_topology: md.Topology):
    """Whether torsions are proper or improper dihedrals.

    Parameters
    ----------
    torsions : np.ndarray or Sequence[Sequence[int]]
        A list of torsions or a zmatrix.
    mdtraj_topology : md.Topology

    Returns
    -------
    is_proper : np.ndarray
        A boolean array which contains 1 for proper and 0 for improper torsions.
    """
    is_proper = np.zeros(len(torsions), dtype=bool)
    import networkx as nx
    graph = mdtraj_topology.to_bondgraph()
    for i, torsion in enumerate(torsions):
        if -1 in torsion:
            continue
        atoms = [mdtraj_topology.atom(torsion[i]) for i in range(4)]
        if (
                atoms[1] in graph.neighbors(atoms[0])
            and atoms[2] in graph.neighbors(atoms[1])
            and atoms[3] in graph.neighbors(atoms[2])
        ):
            is_proper[i] = True
    return is_proper


def rewire_chiral_torsions(z_matrix: np.ndarray, mdtraj_topology: md.Topology, verbose=True):
    """Redefine all torsions containing CB as CB-CA-N-C torsions.
    We use those to define chirality of aminoacids.

    Parameters
    ----------
    z_matrix : np.ndarray
        A matrix of atom indices defining the internal coordinate transform.
    mdtraj_topology : md.Topology

    Returns
    -------
    z_matrix : np.ndarray
        A modified z-matrix that has all CB positions conditioned on CA, N, C
    indices : np.ndarray
        Indices of the chirality-defining torsions in the z-matrix.
    """
    cbetas = mdtraj_topology.select("name CB")
    found_torsions = []
    for cb in cbetas:
        indices = []
        for i, torsion in enumerate(z_matrix):
            if cb == torsion[0]:
                if not -1 in torsion:
                    indices.append(i)
        if len(indices) != 1:
            warnings.warn(f"Coordinate transform's z-matrix has no unique torsion with CB (index: {cb})")
        found_torsions.append((cb, indices[0]))

    for cb, torsion_index in found_torsions:
        torsion = z_matrix[torsion_index]
        chiral_torsion = [cb]
        for name in ["CA", "N", "C"]:
            atom = mdtraj_topology.atom(cb)
            selection = mdtraj_topology.select(f"resid {atom.residue.index} and name {name}")
            assert len(selection) == 1
            chiral_torsion.append(selection[0])
        chiral_torsion = np.array(chiral_torsion)
        if (chiral_torsion != torsion).any():
            if verbose:
                print(f"rewire torsion {torsion} to {chiral_torsion}")
            z_matrix[torsion_index] = chiral_torsion
    return z_matrix, np.array([index for _, index in found_torsions])
