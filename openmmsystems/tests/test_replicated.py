from openmmsystems.systems import ImplicitBPTI
from openmmsystems.replicated import ReplicatedSystem

from simtk.openmm import NonbondedForce


def test_replicated():
    bpti = ImplicitBPTI()
    bpti.system.removeForce(bpti.system.getNumForces() - 1)  # retain only vacuum part
    bpti.system.getForce(bpti.system.getNumForces() - 1).setNonbondedMethod(NonbondedForce.NoCutoff)

    n_replicas = 2
    s = ReplicatedSystem(bpti.system, n_replicas)
    assert s.system.getNumConstraints() == bpti.system.getNumConstraints() * n_replicas
    assert s.system.getNumParticles() == bpti.system.getNumParticles() * n_replicas

    s = ReplicatedSystem(bpti.system, n_replicas, enable_energies=True)