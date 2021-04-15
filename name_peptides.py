
import os
from shutil import copy
from bgmol.util import AMINO_ACIDS

for i, aa in enumerate(AMINO_ACIDS):
    for j, aa2 in enumerate(AMINO_ACIDS):
        vac = f"bgmol/data/dipeptides/dpep_{i*len(AMINO_ACIDS)+j:04d}.pdb"
        sol = f"bgmol/data/dipeptides/dpep_eq_{i*len(AMINO_ACIDS)+j:04d}.pdb"
        assert os.path.exists(vac)
        assert os.path.exists(sol)
        with open(vac, "r") as f:
            lines = f.readlines()
            print(aa, lines[13].split()[3], aa2, lines[-10].split()[3])
        copy(vac, f"bgmol/data/minipeptides/{aa}{aa2}.pdb")
        copy(sol, f"bgmol/data/minipeptides/{aa}{aa2}_solvated.pdb")