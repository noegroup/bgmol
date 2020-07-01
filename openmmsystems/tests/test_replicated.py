from openmmsystems.base import OpenMMToolsTestSystem
from openmmsystems.replicated import ReplicatedSystem

import pytest
import numpy as np
from simtk import unit
from simtk.openmm import NonbondedForce, Context, VerletIntegrator, Platform


def test_replicated():
    base_system = OpenMMToolsTestSystem("AlanineDipeptideVacuum")
    base_system.system.getForces()

    n_replicas = 2
    s = ReplicatedSystem(base_system.system, n_replicas)
    assert s.system.getNumConstraints() == base_system.system.getNumConstraints() * n_replicas
    assert s.system.getNumParticles() == base_system.system.getNumParticles() * n_replicas

    s = ReplicatedSystem(base_system.system, n_replicas, enable_energies=False)

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