"""
High-level API
"""

from openmmsystems import _openmmtools_testsystems, systems, base, replicated
from openmmsystems.util import OpenMMSystemsException, yaml_load
import inspect


def get_system_by_name(name, **kwargs):
    """Create an OpenMMSystem object from its name and some keyword arguments."""
    if name in get_openmmsystems_names():
        return getattr(systems, name)(**kwargs)
    elif name in get_openmmtools_system_names():
        return base.OpenMMToolsTestSystem(name, **kwargs)
    elif name == "ReplicatedSystem":
        base_system_name = kwargs.pop("base_system_name")
        n_replicas = kwargs.pop("n_replicas")
        enable_energies = kwargs.pop("enable_energies")
        base_system = get_system_by_name(base_system_name, **kwargs)
        return replicated.ReplicatedSystem(base_system, n_replicas, enable_energies)
    else:
        raise OpenMMSystemsException(f"System {name} not found.")


def get_system_by_yaml(stream):
    """Create an OpenMMSystem object from its yaml definition in an input stream."""
    system_dict = yaml_load(stream)
    assert "system" in system_dict
    return get_system_by_name(name=system_dict["system"]["identifier"], **system_dict["system"]["parameters"])


def get_system_names():
    """
    List names of all available systems.

    Returns
    -------
    names (list of str):
        The names of all available systems.
    """
    return get_toysystem_names() + get_openmmsystems_names() + get_openmmtools_system_names()


def get_toysystem_names():
    return []


def get_openmmtools_system_names():
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


def get_openmmsystems_names():
    """
    List names of testsystems in this module.

    Returns
    -------
    names (list of str):
        The names all subclasses of openmmtools.TestSystem
    """
    names = []
    for key, value in systems.__dict__.items():
        if (
                inspect.isclass(value)
                and issubclass(value, systems.OpenMMSystem)
                and key not in base.__dict__.keys()  # not one of the base classes
        ):
            names.append(key)
    return names
