"""
High-level API
"""

from openmmsystems import _openmmtools_testsystems, systems, base
import inspect


def get_system_names():
    """
    List names of all available systems.

    Returns
    -------
    names (list of str):
        The names of all available systems.
    """
    return get_openmmtools_system_names() + get_openmmsystems_names()


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
