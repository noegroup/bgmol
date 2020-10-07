from openmmsystems.systems import OpenMMToolsTestSystem
from openmmsystems.systems.replicated import ReplicatedSystem
from openmmsystems.api import system_by_yaml

import pytest
import numpy as np
from simtk import unit
from simtk.openmm import Context, VerletIntegrator, Platform, app


@pytest.mark.parametrize(
    "base_system",
     [
         OpenMMToolsTestSystem("AlanineDipeptideVacuum"),
         #OpenMMToolsTestSystem("AlanineDipeptideImplicit")
     ]
)
@pytest.mark.parametrize("enable_energies", [True, False])
def test_replicated(base_system, enable_energies):

    n_replicas = 2
    s = ReplicatedSystem(base_system, n_replicas)
    assert s.system.getNumConstraints() == base_system.system.getNumConstraints() * n_replicas
    assert s.system.getNumParticles() == base_system.system.getNumParticles() * n_replicas

    s = ReplicatedSystem(base_system, n_replicas, enable_energies=enable_energies)

    context1 = Context(base_system.system, VerletIntegrator(0.001), Platform.getPlatformByName("Reference"))
    context1.setPositions(base_system.positions)
    ener1 = context1.getState(getEnergy=True).getPotentialEnergy().value_in_unit_system(unit.md_unit_system)

    positions2 = 1.01*base_system.positions
    context2 = Context(base_system.system, VerletIntegrator(0.001), Platform.getPlatformByName("Reference"))
    context2.setPositions(positions2)
    ener2 = context2.getState(getEnergy=True).getPotentialEnergy().value_in_unit_system(unit.md_unit_system)

    context3 = Context(s.system, VerletIntegrator(0.001), Platform.getPlatformByName("Reference"))
    context3.setPositions(np.row_stack([base_system.positions, positions2]))
    ener3 = context3.getState(getEnergy=True).getPotentialEnergy().value_in_unit_system(unit.md_unit_system)
    assert ener3 == pytest.approx(ener1 + ener2)

    if enable_energies:
        ener3a = context3.getState(getEnergy=True, groups={0}).getPotentialEnergy().value_in_unit_system(unit.md_unit_system)
        ener3b = context3.getState(getEnergy=True, groups={1}).getPotentialEnergy().value_in_unit_system(unit.md_unit_system)
        assert ener3a == pytest.approx(ener1)
        assert ener3b == pytest.approx(ener2)


def test_load_replicated():
    ala2 = OpenMMToolsTestSystem("AlanineDipeptideVacuum", constraints=app.HBonds)
    n_replicas = 2
    s = ReplicatedSystem(ala2, n_replicas)
    s2 = system_by_yaml(str(s))
    assert isinstance(s2._base_system, OpenMMToolsTestSystem)
    assert (str(s) == str(s2))
