bgmol
==============================
[//]: # (Badges)
[![Travis Build Status](https://travis-ci.com/REPLACE_WITH_OWNER_ACCOUNT/ProjectName.svg?branch=master)](https://travis-ci.com/REPLACE_WITH_OWNER_ACCOUNT/ProjectName)
[![AppVeyor Build status](https://ci.appveyor.com/api/projects/status/REPLACE_WITH_APPVEYOR_LINK/branch/master?svg=true)](https://ci.appveyor.com/project/REPLACE_WITH_OWNER_ACCOUNT/ProjectName/branch/master)
[![codecov](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/ProjectName/branch/master/graph/badge.svg)](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/ProjectName/branch/master)


A collection of OpenMM systems and data generated for those systems.

This package is an extension of the `openmmtools.testsystems` module with a twofold aim:
1) It contains a collection of OpenMM systems (with initial positions and topologies).
2) It provides infrastructure for sharing data (positions, forces, energies, velocities)
generated for these systems.


### Quickstart: Install
Make sure that you are in a conda environment that has **OpenMM** and **mdtraj** installed.

Clone the code
```
git clone git@github.com:noegroup/bgmol.git
```

And install
```
cd bgmol
python setup.py install
```


### Quickstart: Example

```python
from bgmol.datasets import Ala2Implicit300
dataset = Ala2Implicit300(download=True, read=True)
```

The dataset contains forces, energies and coordinates
it also holds a reference to the system that defines the potential energy function.
```python
openmmsystem = dataset.system
```

The system is an `OpenMMSystem` object, it provides access to the simtk.openmm.system instance,
the topology, and a set of initial coordinates. For example, we can run an OpenMM simulation
as follows
```python
from simtk.openmm.app import Simulation, LangevinIntegrator
integrator = LangevinIntegrator(dataset.temperature, 1, 0.001)
simulation = Simulation(openmmsystem.topology, openmmsystem.system, integrator)
simulation.context.setPositions(openmmsystem.positions)
simulation.step(10)
```

The dataset contains coordinates (xyz), forces, and energies.
```python
dataset.energies
```


# get data and OpenMMSystem instance
bpti_data = get_data("bpti_implicit_mini_example")
bpti = bpti_data.create_openmm_system()
print(bpti_data.info())

# download data
bpti_data.download()
bpti_data.read()
print(bpti_data.positions)
trajectory = bpti_data.to_mdtraj()

# evaluate energies
from simtk.openmm import LangevinIntegrator
energy_bridge = OpenMMEnergyBridge(bpti.system, LangevinIntegrator(300,1,0.0002))
energies, force = energy_bridge.evaluate(bpti_data.positions)
```

### Using Systems from this Repository

Each `OpenMMSystem` has a system object, a topology, a set of initial positions, and a unique identifier (the class 
or directory name).

Subclasses of `OpenMMSystem` can be created via their constructor:
```python
from bgmol import BPTIImplicit
bpti = BPTIImplicit()
```

All available systems can also be created via their identifier:

```python
from bgmol import get_system
ala = get_system("AlanineDipeptideImplicit")  # this is a system defined in openmmtools_testsystems
print(ala.info())
```

Identifiers can be:
1) Class names in `bgmol.systems` (subclasses of `bgmol.OpenMMSystem`)
2) Directory names in `bgmol/systems` (each subdirectory defines a system, as specified below)
3) Class names in `bgmol/openmmtools_testsystems` (subclasses of `openmmtools_testsystems.TestSystem`)

`OpenMMSystem` instances may have keyword arguments, for example:

```python
from simtk import unit
from bgmol import system_by_name

system_by_name(
    "HarmonicOscillator", 
    K=100.0*unit.kilocalories_per_mole / unit.angstroms**2, 
    mass=39.948*unit.amu
)
```


### Using Data from this Repository

The samples created for the different systems do not live in the repository 
but in some path (webserver, compute cluster, local network, local machine).
This means that not all samples are accessible for all users; they may require access to specific clusters, etc.

You can list all registered data (sets of positions, forces, energies, velocities) for a system:

```python
from bgmol import list_all_data
list_all_data("AlanineDipeptideImplicit")
```

You can also list only the data sets that are accessible for you from your machine:

```python
from bgmol import list_accessible_data
list_accessible_data("AlanineDipeptideImplicit")
```

or data that is available for a given set of constructor arguments

```python
from bgmol import list_accessible_data
ala_samples = list_accessible_data("HarmonicOscillator", mass=39.948*unit.amu)
```

Samples also have unique identifiers (the name of a yaml file, as specified below) and can be read as follows 
```python
bpti_data = get_data("bpti_implicit_mini_example")
print(bpti_data.info())
bpti_data.download()  # gets stored in $HOME/bgmol_data/datasets
bpti_data.read() # load data into memory
bpti_data.randomize() # random permutation

# Now, randomly reordered positions and forces are accessible as numpy arrays
bpti_data.positions
bpti_data.forces
```

For each data set, you can also get the system
```python
bpti_data = bpti_data.create_openmm_system()
```


### Adding Systems to the Repository

The preferred way to add systems is by adding a subclass of `OpenMMSystem` to the `bgmol.systems` module, 
see examples therein. 

### Adding Samples to the Repository


### Customization


### Evaluating Energies

TODO: The BGTorch openmm wrapper would have a natural place in this repository.

### Creating Samples

TODO: simple API to extend a data set.

### Copyright

Copyright (c) 2020, noegroup


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.2.
