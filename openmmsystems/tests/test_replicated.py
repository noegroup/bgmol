from openmmsystems.base import OpenMMToolsTestSystem
from openmmsystems.replicated import ReplicatedSystem

import pytest
import numpy as np
from simtk import unit
from simtk.openmm import NonbondedForce, Context, VerletIntegrator, Platform


def test_replicated():
    base_system = OpenMMToolsTestSystem("AlanineDipeptideVacuum")

    n_replicas = 3
    s = ReplicatedSystem(base_system.system, n_replicas)
    assert s.system.getNumConstraints() == base_system.system.getNumConstraints() * n_replicas
    assert s.system.getNumParticles() == base_system.system.getNumParticles() * n_replicas

    s = ReplicatedSystem(base_system.system, n_replicas, enable_energies=True)

    context1 = Context(base_system.system, VerletIntegrator(0.001), Platform.getPlatformByName("Reference"))
    context1.setPositions(base_system.positions)
    ener1 = context1.getState(getEnergy=True).getPotentialEnergy().value_in_unit_system(unit.md_unit_system)

    context2 = Context(s.system, VerletIntegrator(0.001), Platform.getPlatformByName("Reference"))
    context2.setPositions(np.row_stack([base_system.positions]*n_replicas))
    ener2 = context2.getState(getEnergy=True).getPotentialEnergy().value_in_unit_system(unit.md_unit_system)
    assert ener2 == pytest.approx(ener1*n_replicas)