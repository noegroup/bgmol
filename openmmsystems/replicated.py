
__all__ = ["ReplicatedSystem"]


from simtk.openmm import (
    System,
    LocalCoordinatesSite,
    OutOfPlaneSite,
    TwoParticleAverageSite,
    ThreeParticleAverageSite,
    HarmonicBondForce,
    HarmonicAngleForce,
    PeriodicTorsionForce,
    NonbondedForce,
    CustomBondForce,
    CustomNonbondedForce,
)


class ReplicatedSystem:
    """
    """

    def __init__(self, base_system, n_replicas, enable_energies=False):
        self._base_system = base_system
        self._n_replicas = n_replicas
        self._system = self._replicate_system(base_system, n_replicas, enable_energies)

    @property
    def system(self):
        return self._system

    @staticmethod
    def _replicate_system(base_system, n_replicas, enable_energies):
        system = System()
        n_particles = base_system.getNumParticles()
        # particles
        for j in range(n_replicas):
            for i in range(n_particles):
                system.addParticle(base_system.getParticleMass(i))
                if system.isVirtualSite(i):
                    vs = system.getVirtualSite(i)
                    vs_copy = self._replicate_virtual_site(vs, n_particles, j)
                    system.setVirtualSite(i + j * n_particles, vs_copy)
        # constraints
        for j in range(n_replicas):
            for i in range(base_system.getNumConstraints()):
                p1, p2, distance = base_system.getConstraintParameters(i)
                system.addConstraint(p1 + j * n_particles, p2 + j * n_particles, distance)
        # properties
        system.setDefaultPeriodicBoxVectors(*(base_system.getDefaultPeriodicBoxVectors()))
        # forces
        for force in base_system.getForces():
            forcename = force.__class__.__name__
            methodname = f"_replicate_{forcename}"
            assert hasattr(ReplicatedSystem, methodname), f"Replicating {forcename} not implemented."
            replicate_force_method = getattr(ReplicatedSystem, methodname)
            replicated_forces = replicate_force_method(force, n_particles, n_replicas, enable_energies)
            for f in replicated_forces:
                system.addForce(f)
        return system

    @staticmethod
    def _replicate_virtual_site(vs, n_particles, replica):
        if isinstance(vs, LocalCoordinatesSite):
            args = []
            for i in range(vs.getNumParticles()):
                particle_ids.append(vs.getParticle(i) + replica * n_particles)
            args.append(vs.getOriginWeights())
            args.append(vs.getXWeights())
            args.append(vs.getYWeights())
            args.append(vs.getLocalPosition())
            return LocalCoordinatesSite(*args)
        elif isinstance(vs, OutOfPlaneSite):
            args = []
            for i in range(vs.getNumParticles()):
                particle_ids.append(vs.getParticle(i) + replica * n_particles)
            args.append(vs.getWeight12())
            args.append(vs.getWeight13())
            args.append(vs.getWeightCross())
            return OutOfPlaneSite(*args)
        elif isinstance(vs, TwoParticleAverageSite):
            return TwoParticleAverageSite(
                vs.getParticle(0) + replica * n_particles,
                vs.getParticle(1) + replica * n_particles,
                vs.getWeight(0),
                vs.getWeight(1)
            )
        elif isinstance(vs, ThreeParticleAverageSite):
            return ThreeParticleAverageSite(
                vs.getParticle(0) + replica * n_particles,
                vs.getParticle(1) + replica * n_particles,
                vs.getParticle(2) + replica * n_particles,
                vs.getWeight(0),
                vs.getWeight(1),
                vs.getWeight(2)
            )
        else:
            raise f"Unknown virtual site type: {type(vs)}."

    @staticmethod
    def _replicate_HarmonicBondForce(force, n_particles, n_replicas, enable_energies):
        replicated_forces = []
        replicated_force = HarmonicBondForce()
        replicated_force.setUsesPeriodicBoundaryConditions(force.usesPeriodicBoundaryConditions())
        for j in range(n_replicas):
            for i in range(force.getNumBonds()):
                p1, p2, length, k = force.getBondParameters(i)
                replicated_force.addBond(p1 + j * n_particles, p2 + j * n_particles, length, k)
            if enable_energies:
                # create a new force object for each replicate
                replicated_force.setForceGroup(j)
                replicated_forces.append(replicated_force)
                replicated_force = HarmonicBondForce()
                replicated_force.setUsesPeriodicBoundaryConditions(force.usesPeriodicBoundaryConditions())
        return replicated_forces

    @staticmethod
    def _replicate_HarmonicAngleForce(force, n_particles, n_replicas, enable_energies):
        replicated_forces = []
        replicated_force = HarmonicAngleForce()
        replicated_force.setUsesPeriodicBoundaryConditions(force.usesPeriodicBoundaryConditions())
        for j in range(n_replicas):
            for i in range(force.getNumAngles()):
                p1, p2, p3, angle, k = force.getAngleParameters(i)
                replicated_force.addAngle(
                    p1 + j * n_particles,
                    p2 + j * n_particles,
                    p3 + j * n_particles,
                    angle, k)
            if enable_energies:
                # create a new force object for each replicate
                replicated_force.setForceGroup(j)
                replicated_forces.append(replicated_force)
                replicated_force = HarmonicAngleForce()
                replicated_force.setUsesPeriodicBoundaryConditions(force.usesPeriodicBoundaryConditions())
        return replicated_forces

    @staticmethod
    def _replicate_PeriodicTorsionForce(force, n_particles, n_replicas, enable_energies):
        replicated_forces = []
        replicated_force = PeriodicTorsionForce()
        replicated_force.setUsesPeriodicBoundaryConditions(force.usesPeriodicBoundaryConditions())
        for j in range(n_replicas):
            for i in range(force.getNumTorsions()):
                p1, p2, p3, p4, angle, mult, k = force.getTorsionParameters(i)
                replicated_force.addTorsion(
                    p1 + j * n_particles,
                    p2 + j * n_particles,
                    p3 + j * n_particles,
                    p4 + j * n_particles,
                    angle, mult, k)
            if enable_energies:
                # create a new force object for each replicate
                replicated_force.setForceGroup(j)
                replicated_forces.append(replicated_force)
                replicated_force = PeriodicTorsionForce()
                replicated_force.setUsesPeriodicBoundaryConditions(force.usesPeriodicBoundaryConditions())
        return replicated_forces

    @staticmethod
    def _replicate_NonbondedForce(force, n_particles, n_replicas, enable_energies):
        nonbonded_method = force.getNonbondedMethod()
        if nonbonded_method == NonbondedForce.NoCutoff:
            return ReplicatedSystem._replicate_nonbonded_as_custom_bond_force(
                force,
                n_particles,
                n_replicas,
                enable_energies
            )
        else:
            return NotImplemented

    @staticmethod
    def _replicate_nonbonded_as_custom_bond_force(force, n_particles, n_replicas, enable_energies):
        replicated_forces = []
        energy_string = ""
        replicated_force = CustomBondForce()
        force.usesPeriodicBoundaryConditions()
