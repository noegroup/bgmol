"""
Toy systems that do not require an MD engine.
"""


from bgflow.distribution.energy import Energy
from bgflow.utils.geometry import distance_vectors, distances_from_vectors


class DoubleWellEnergy(Energy):
    def __init__(self, dim, a=0, b=-4., c=1.):
        super().__init__(dim)
        self._a = a
        self._b = b
        self._c = c

    def _energy(self, x):
        d = x[:, [0]]
        v = x[:, 1:]
        e1 = self._a * d + self._b * d.pow(2) + self._c * d.pow(4)
        e2 = 0.5 * v.pow(2).sum(dim=-1, keepdim=True)
        return e1 + e2


class MultiDoubleWellPotential(Energy):
    def __init__(self, dim, n_particles, a, b, c, offset):
        super().__init__(dim)
        self._dim = dim
        self._n_particles = n_particles
        self._a = a
        self._b = b
        self._c = c
        self._offset = offset

    def _energy(self, x):
        n_batch = x.shape[0]
        x = x.view(n_batch, self._n_particles, -1)
        dists = distances_from_vectors(
            distance_vectors(x.view(n_batch, self._n_particles, -1))
        )
        dists = dists.view(-1, 1)
        dists = dists - self._offset
        energies = self._c * dists**4 + self._b * dists**2 + self._a
        return energies.view(n_batch, -1).sum(-1) / 2

