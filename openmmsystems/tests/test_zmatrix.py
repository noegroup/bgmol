import pytest
import torch
from simtk import unit
from openmmsystems.systems import ImplicitBPTI, MiniPeptide, AlanineDipeptideImplicit
from openmmsystems.zmatrix import make_protein_z_matrix
from bgtorch import (
    RelativeInternalCoordinateTransformation,
    MixedCoordinateTransformation,
    GlobalInternalCoordinateTransformation
)


@pytest.mark.parametrize("system", [
    #AlanineDipeptideImplicit(),
    ImplicitBPTI()
])
def test_z_matrix_mixed_trafo(system):
    top = system.mdtraj_topology
    zmatrix, fixed = make_protein_z_matrix(top, "backbone")
    #print([list(top.atoms)[i] for i in top.select("backbone")])
    data = torch.tensor(system.positions.value_in_unit(unit.nanometer)).reshape(1,-1)
    MixedCoordinateTransformation(data, zmatrix, fixed)


@pytest.mark.parametrize("system", [
    #AlanineDipeptideImplicit(),
    ImplicitBPTI()
])
def test_z_matrix_global_trafo(system):
    top = system.mdtraj_topology
    zmatrix = make_protein_z_matrix(top, cartesian=None)
    #trafo = GlobalInternalCoordinateTransformation(zmatrix)