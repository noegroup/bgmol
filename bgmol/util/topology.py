
import warnings
from typing import Sequence
import numpy as np
import mdtraj as md


__all__ = [
    "rewire_chiral_torsions", "find_rings", "is_proper_torsion", "is_chiral_torsion",
    "is_ring_torsion", "is_methyl_torsion", "is_type_torsion", "is_ramachandran_torsion"
]


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
    for i, torsion in enumerate(torsions):
        # checking whether the atom that is placed by this torsion is part of any rings
        ring_indices = [j for j, ring in enumerate(rings) if torsion[0] in ring]
        if len(ring_indices) == 0:
            continue
        # checking whether any other atom in this torsion is part of the same ring
        for ring_index in ring_indices:
            if any(atom in rings[ring_index] for atom in torsion[1:]):
                is_ring[i] = True
                continue
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


def is_methyl_torsion(torsions: Sequence[Sequence[int]], mdtraj_topology: md.Topology):
    """Whether torsions are the first (proper) torsion of a methyl group.
    Methyl hydrogens are placed by one proper and two improper torsions. This function only indicates the former.

    Parameters
    ----------
    torsions : np.ndarray or Sequence[Sequence[int]]
        A list of torsions or a zmatrix.
    mdtraj_topology : md.Topology

    Returns
    -------
    is_methyl : np.ndarray
        A boolean array which contains 1 for proper and 0 for improper torsions.

    Notes
    -----
    Requires networkx.
    """
    is_methyl = np.zeros(len(torsions), dtype=bool)
    is_proper = is_proper_torsion(torsions, mdtraj_topology)
    graph = mdtraj_topology.to_bondgraph()
    for i, (torsion, proper) in enumerate(zip(torsions, is_proper)):
        if not proper:
            continue
        atom = mdtraj_topology.atom(torsion[0])
        if not atom.element.symbol == "H":
            continue
        neighbor = list(graph.neighbors(atom))[0]
        if not neighbor.element.symbol == "C":
            continue
        carbon_neighbors = graph.neighbors(neighbor)
        n_hydrogens = sum(n.element.symbol == "H" for n in carbon_neighbors)
        if n_hydrogens == 3:
            is_methyl[i] = True
    return is_methyl


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


def is_type_torsion(type_torsion: str, torsions: Sequence[Sequence[int]], mdtraj_topology: md.Topology):
    """Whether torsions are of the specified type.
    Types supported are: ramachandran, phi, psi, omega, chi1, chi2, chi3, chi4, ring, proper, methyl, chiral

    Parameters
    ----------
    type_torsion : str
        Which type of torsion among the supported ones.
    torsions : np.ndarray or Sequence[Sequence[int]]
        A list of torsions or a zmatrix.
    mdtraj_topology : md.Topology

    Returns
    -------
    is_type : np.ndarray
        A boolean array which contains True iff the torsion is of the specified type.
    """

    fake_traj = md.Trajectory(np.zeros((mdtraj_topology, 3)), mdtraj_topology)
    if type_torsion == 'ramachandran':
        indices = np.vstack((md.compute_phi(fake_traj)[0], md.compute_psi(fake_traj)[0]))
    if type_torsion == 'phi':
        indices = md.compute_phi(fake_traj)[0]
    elif type_torsion == 'psi':
        indices = md.compute_psi(fake_traj)[0]
    elif type_torsion == 'omega':
        indices = md.compute_omega(fake_traj)[0]
    elif type_torsion == 'chi1':
        indices = md.compute_chi1(fake_traj)[0]
    elif type_torsion == 'chi2':
        indices = md.compute_chi2(fake_traj)[0]
    elif type_torsion == 'chi3':
        indices = md.compute_chi3(fake_traj)[0]
    elif type_torsion == 'chi4':
        indices = md.compute_chi4(fake_traj)[0]
    elif type_torsion == 'ring':
        return is_ring_torsion(torsions, mdtraj_topology)
    elif type_torsion == 'proper':
        return is_proper_torsion(torsions, mdtraj_topology)
    elif type_torsion == 'methyl':
        return is_methyl_torsion(torsions, mdtraj_topology)
    elif type_torsion == 'chiral':
        return is_chiral_torsion(torsions, mdtraj_topology)
    else:
        raise ValueError("Supported torsion types are: "
          "'ramachandran', 'phi', 'psi', 'omega', 'chi1', 'chi2', "
          "'chi3', 'chi4', 'ring', 'proper', 'methyl', 'chiral'")

    ordered_torsions = np.sort(torsions, axis=1) #make sure index ordering is same as md
    is_type = np.array([np.any([np.all(j == ind) for j in indices]) for ind in ordered_torsions])

    return is_type


def is_ramachandran_torsion(torsions: Sequence[Sequence[int]], mdtraj_topology: md.Topology):
    """Whether torsions are Ramachandran angles or not.

    Parameters
    ----------
    torsions : np.ndarray or Sequence[Sequence[int]]
        A list of torsions or a zmatrix.
    mdtraj_topology : md.Topology

    Returns
    -------
    is_ramachandran : np.ndarray
        A boolean array which contains True iff the torsion is Ramachandran
    """
  return is_type_torsion('ramachandran', torsions, mdtraj_topology)


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
