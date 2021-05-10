
import os
import pytest
import shutil

import numpy as np
from simtk.openmm import Context, VerletIntegrator, Platform
from bgmol import datasets
from bgmol.api import list_datasets


@pytest.fixture(params=list_datasets())
def all_datasets(request, tmpdir):
    subdir = os.path.join(tmpdir, request.param)
    os.mkdir(subdir)
    yield request.param, subdir
    shutil.rmtree(subdir)


@pytest.mark.slow
def test_all_datasets(all_datasets):
    """Try downloading all datasets and compare energies and forces for a random frame."""
    dataset_name, dataset_dir = all_datasets
    DatasetClass = getattr(datasets, dataset_name)
    dataset = DatasetClass(root=dataset_dir, read=True, download=True)
    system = dataset.system

    # check shapes
    atoms = system.select(dataset.selection)
    assert dataset.xyz.shape == (dataset.num_frames, len(atoms), 3)
    if dataset.forces is not None:
        assert dataset.forces.shape == (dataset.num_frames, len(atoms), 3)
    if dataset.energies is not None:
        assert dataset.energies.shape == (dataset.num_frames, )

    # check energies and forces
    if dataset.selection == "all":
        num_comparisons = 10
        for i in range(num_comparisons):
            random_frame = np.random.randint(0, dataset.num_frames)
            # check energy/force for a random frame
            context = Context(
                system.system,
                VerletIntegrator(0.001),
                Platform.getPlatformByName("Reference")
            )
            context.setPositions(dataset.xyz[random_frame])
            if dataset.unitcell_vectors is not None:
                context.setPeriodicBoxVectors(*dataset.unitcell_vectors[random_frame])
            state = context.getState(getEnergy=True, getForces=dataset.forces is not None)
            if dataset.energies is not None:
                assert dataset.energies[random_frame] == pytest.approx(state.getPotentialEnergy()._value, abs=3e-1, rel=0.0)
            if dataset.forces is not None:
                assert np.allclose(dataset.forces[random_frame], state.getForces(asNumpy=True)._value, atol=3e-1, rtol=0.0)

