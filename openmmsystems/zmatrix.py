

import numpy as np

__all__ = ["make_protein_z_matrix"]


def make_global_z_matrix():
    pass


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
