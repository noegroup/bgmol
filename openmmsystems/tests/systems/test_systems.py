"""Tests for systems in this package."""

import pytest

from openmmsystems.systems import *

from simtk.openmm import Context, VerletIntegrator, Platform
from simtk import unit
import numpy as np

# ===== general tests ======
# The following tests are run for all OpenMMSystem subclasses defined in openmmsystems.systems


def test_default_constructor(system_instance):
    """test construction of all systems"""
    pass


def test_num_particles(system_instance):
    """test arguments of all systems"""
    num_particles = system_instance.system.getNumParticles()
    assert num_particles > 0
    assert system_instance.mdtraj_topology.n_atoms == num_particles
    assert len(system_instance.positions) == num_particles


def test_energy_calculation(system_instance):
    """test calculation of initial energy"""
    platform = Platform.getPlatformByName("Reference")
    context = Context(system_instance.system, VerletIntegrator(0.001), platform)
    context.setPositions(system_instance.positions)
    state = context.getState(getEnergy=True)
    energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    assert not np.isnan(energy) and not np.isinf(energy)


# ====== specific tests =======
# tests for specific systems

@pytest.mark.parametrize("ff", [['amber99sbildn.xml', 'amber99_obc.xml']])
def test_bpti_implicit_parameters(ff):
    bpti = ImplicitBPTI(forcefield=ff)
