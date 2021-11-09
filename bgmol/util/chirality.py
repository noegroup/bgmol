
__all__ = ["rewire_chiral_torsions"]


def rewire_chiral_torsions(z_matrix, mdtraj_topology):
    ha = mdtraj_topology.select("name HA")
    all_indices = []
    for h in ha:
        indices = []
        for i, torsion in enumerate(z_matrix):
            print(h, torsion, i)
            if h == torsion[0]:
                indices.append(i)
        if len(indices) != 1:
            raise ValueError(f"Coordinate transform's z-matrix has no unique torsion with HA (index: {h})")
        all_indices.append(indices[0])
    # sanity check
    atoms = list(mdtraj_topology.atoms)
    for h, torsion_index in zip(ha, all_indices):
        torsion = z_matrix[torsion_index]
        # same residue
        for atom in torsion:
            assert atoms[atom].residue == atoms[h].residue
        # correct names
        assert atoms[torsion[1]].name == "CA"
        assert atoms[torsion[2]].name == "N"
        assert atoms[torsion[3]].name == "C"
    return all_indices
