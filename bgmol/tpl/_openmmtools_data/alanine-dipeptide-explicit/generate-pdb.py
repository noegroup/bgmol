"""
Generate PDB file containing periodic box data.

"""

from ...util.importing import import_openmm
mm, unit, app = import_openmm()

prmtop_filename = 'alanine-dipeptide.prmtop'
crd_filename = 'alanine-dipeptide.crd'
pdb_filename = 'alanine-dipeptide.pdb'

# Read topology and positions.
prmtop = app.AmberPrmtopFile(prmtop_filename)
inpcrd = app.AmberInpcrdFile(crd_filename)

# Write PDB.
outfile = open(pdb_filename, 'w')
app.PDBFile.writeFile(prmtop.topology, inpcrd.positions, file=outfile, keepIds=False)
outfile.close()

