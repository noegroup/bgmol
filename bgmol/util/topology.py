
import warnings
from typing import Sequence
import numpy as np
import mdtraj as md


__all__ = ["rewire_chiral_torsions", "find_rings", "is_proper_torsion", "is_chiral_torsion", "is_ring_torsion"]


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


def is_ring_torsion(torsions: Sequence[Sequence[int]], mdtraj_topology: md.Topology):
    """Whether torsions are part of defining a ring.

    Returns
    -------
    is_ring : np.ndarray
        A boolean array which contains 1 for torsions of ring atoms and 0 for others
    """
    is_ring = np.zeros(len(torsions), dtype=bool)
    rings = find_rings(mdtraj_topology)
    ring_atoms = set(atom for ring in rings for atom in ring)
    for i, torsion in enumerate(torsions):
        if torsion[0] in ring_atoms and any(atom in ring_atoms for atom in torsion[1:]):
            is_ring[i] = True
    return is_ring


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

    Notes
    -----
    Requires networkx.
    """
    is_proper = np.zeros(len(torsions), dtype=bool)
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


def _select_ha(mdtraj_topology: md.Topology):
    halpha = mdtraj_topology.select("name HA")
    graph = mdtraj_topology.to_bondgraph()
    indices = []
    for ha in halpha:
        atom = mdtraj_topology.atom(ha)
        neighbors = list(graph.neighbors(atom))
        if len(neighbors) != 1:
            continue
        ca = neighbors[0]
        if ca.name != "CA":
            continue
        neighbors = list(graph.neighbors(ca))
        if all(neighbor.element.symbol != "C" for neighbor in neighbors):
            continue
        if all(neighbor.element.symbol != "N" for neighbor in neighbors):
            continue
        indices.append(ha)
    return np.array(indices)


def is_chiral_torsion(torsions: Sequence[Sequence[int]], mdtraj_topology: md.Topology):
    """True for all torsions that define HA positions."""
    halpha = _select_ha(mdtraj_topology)
    is_ha = np.zeros(len(torsions), dtype=bool)
    for i, torsion in enumerate(torsions):
        if torsion[0] in halpha:
            is_ha[i] = True
    return is_ha


def rewire_chiral_torsions(z_matrix: np.ndarray, mdtraj_topology: md.Topology, verbose=True):
    """Redefine all torsions containing HA as HA-CA-N-C torsions.
    We use those to define chirality of aminoacids.

    Parameters
    ----------
    z_matrix : np.ndarray
        A matrix of atom indices defining the internal coordinate transform.
    mdtraj_topology : md.Topology

    Returns
    -------
    z_matrix : np.ndarray
        A modified z-matrix that has all HA positions conditioned on CA, N, C

    Notes
    -----
    The indices do not correspond directly to indices
    """
    if len(z_matrix) == 0:
        return z_matrix
    halphas = _select_ha(mdtraj_topology)
    halphas = np.intersect1d(halphas, z_matrix[:, 0])

    found_torsions = []
    for ha in halphas:
        indices = []
        for i, torsion in enumerate(z_matrix):
            if ha == torsion[0]:
                if not -1 in torsion:
                    indices.append(i)
        if len(indices) != 1:
            warnings.warn(f"Coordinate transform's z-matrix has no unique torsion with HA (index: {ha})")
            continue
        found_torsions.append((ha, indices[0]))

    for ha, torsion_index in found_torsions:
        torsion = z_matrix[torsion_index]
        chiral_torsion = [ha]
        for name in ["CA", "N", "C"]:
            atom = mdtraj_topology.atom(ha)
            selection = mdtraj_topology.select(f"resid {atom.residue.index} and name {name}")
            assert len(selection) == 1
            chiral_torsion.append(selection[0])
        chiral_torsion = np.array(chiral_torsion)
        if (chiral_torsion == torsion).all():
            continue
        # replace torsion by chiral_torsion
        if verbose:
            print(f"replace torsion {torsion} by {chiral_torsion}")

        # remove circular dependencies
        for other_atom in chiral_torsion[1:]:
            if other_atom not in z_matrix[:, 0]:
                continue
            other_index = np.where(z_matrix[:, 0] == other_atom)[0][0]
            if ha in z_matrix[other_index, 1:]:
                replace_with_atoms_from_torsion = ha
                replaced = False
                while not replaced:
                    try:
                        replacement = np.setdiff1d(z_matrix[replace_with_atoms_from_torsion], z_matrix[other_index])[0]
                        replaced = True
                    except IndexError:
                        replace_with_atoms_from_torsion = z_matrix[replace_with_atoms_from_torsion][1]
                        continue
                    atom_index = np.where(z_matrix[other_index] == ha)[0][0]
                    if verbose:
                        print(f"- change {z_matrix[other_index]} to {replacement} at index {atom_index}")
                    z_matrix[other_index, atom_index] = replacement

        z_matrix[torsion_index] = chiral_torsion
    return z_matrix
