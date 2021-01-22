"""
High-level API
"""

import inspect
from openmmsystems import systems
from openmmsystems import datasets
from openmmsystems.tpl import _openmmtools_testsystems
from openmmsystems.systems import replicated
from openmmsystems.util import OpenMMSystemsException, yaml_load


__all__ = ["system_by_name", "system_by_yaml", "list_datasets", "list_openmmsystems", "list_openmmtools_systems",
           "list_systems", "list_toysystems"]


def system_by_name(name, **kwargs):
    """Create an OpenMMSystem object from its name and some keyword arguments."""
    if name in list_openmmsystems():
        return getattr(systems, name)(**kwargs)
    elif name in list_openmmtools_systems():
        return systems.base.OpenMMToolsTestSystem(name, **kwargs)
    elif name == "ReplicatedSystem":
        base_system_name = kwargs.pop("base_system_name")
        n_replicas = kwargs.pop("n_replicas")
        enable_energies = kwargs.pop("enable_energies")
        base_system = system_by_name(base_system_name, **kwargs)
        return replicated.ReplicatedSystem(base_system, n_replicas, enable_energies)
    else:
        raise OpenMMSystemsException(f"System {name} not found.")


def system_by_yaml(stream):
    """Create an OpenMMSystem object from its yaml definition in an input stream."""
    system_dict = yaml_load(stream)
    assert "system" in system_dict
    if "parameters" in system_dict["system"]:
        parameters = system_dict["system"]["parameters"]
        if parameters is None:
            parameters = {}
    else:
        parameters = {}
    return system_by_name(name=system_dict["system"]["identifier"], **parameters)


def list_systems():
    """
    List names of all available systems.

    Returns
    -------
    names (list of str):
        The names of all available systems.
    """
    return list_toysystems() + list_openmmsystems() + list_openmmtools_systems()


def list_toysystems():
    return []


def list_openmmtools_systems():
    """
    List names of openmmtools testsystems.

    Returns
    -------
    names (list of str):
        The names all subclasses of openmmtools.TestSystem
    """
    names = []
    for key, value in _openmmtools_testsystems.__dict__.items():
        if (
                inspect.isclass(value)
                and issubclass(value, _openmmtools_testsystems.TestSystem)
                and value != _openmmtools_testsystems.TestSystem
        ):
            names.append(key)
    return names


def list_openmmsystems():
    """
    List names of OpenMMSystems in this package.

    Returns
    -------
    names (list of str):
        The names all subclasses of systems.OpenMMSystem
    """
    names = []
    for key, value in systems.__dict__.items():
        if (
                inspect.isclass(value)
                and issubclass(value, systems.OpenMMSystem)
                and key not in systems.base.__dict__.keys()  # not one of the base classes
        ):
            names.append(key)
    return names


def list_datasets():
    """
    List names of datasets in this package.

    Returns
    -------
    names (list of str):
        The names all subclasses of datasets
    """
    names = []
    for key, value in datasets.__dict__.items():
        if (
                inspect.isclass(value)
                and issubclass(value, datasets.DataSet)
                and key not in datasets.base.__dict__.keys()  # not one of the base classes
        ):
            names.append(key)
    return names