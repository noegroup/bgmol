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

### Adding Systems to the Repository

The preferred way to add systems is by adding a subclass of `OpenMMSystem` to the `bgmol.systems` module, 
see examples therein. 

### Adding Samples to the Repository
TBD


### Copyright

Copyright (c) 2020, noegroup


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.2.
