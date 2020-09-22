"""
High-level API
"""

from openmmsystems import systems
from openmmsystems.tpl import _openmmtools_testsystems
from openmmsystems.systems import base
from openmmsystems.util import OpenMMSystemsException, yaml_load
import inspect


def system_by_name(name, **kwargs):
    """Create an OpenMMSystem object from its name and some keyword arguments."""
    if name in list_openmmsystems():
        return getattr(systems, name)(**kwargs)
    elif name in list_openmmtools_systems():
        return base.OpenMMToolsTestSystem(name, **kwargs)
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
