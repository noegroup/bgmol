"""
Utility module
"""

import os
from pkg_resources import resource_filename

import yaml

from simtk import unit
from simtk.openmm import app
from simtk.openmm.app.internal.singleton import Singleton

from mdtraj.utils.unit import _str_to_unit

# ==========================================================================================
# YAML utilities
# ==========================================================================================


def parse_quantity(string):
    """
    Parameters
    ----------
    string : str
        A string object representing a quantity.

    Returns
    -------
    quantity : unit.Quantity
        The parsed quantity.
    """
    return _str_to_unit(string, simtk=True)


def quantity_representer(dumper, data):
    """Converts quantity object into a yaml string."""
    if data.unit == unit.dimensionless:
        return dumper.represent_scalar('!quantity', f'{data._value}')
    else:
        operator = '' if data.unit._name.startswith('/') else '*'
        u = data.unit._name.replace(' ','_')
        return dumper.represent_scalar('!quantity', f'{data._value} {operator} {u}')


def quantity_constructor(loader, node):
    """Parser for quantity objects from yaml file. Works only for float and int values."""
    string = loader.construct_scalar(node)
    return parse_quantity(string)


_OPENMM_SINGLETONS = [
    object for object in app.__dict__.values()
    if isinstance(object, Singleton)
]


def singleton_representer(dumper, data):
    """Converts simtk.openmm.app.internal.Singleton object into a yaml string."""
    return dumper.represent_scalar('!openmm', data.__class__.__name__)


def singleton_constructor(loader, node):
    """Parser for simtk.openmm.app.internal.Singleton objects from yaml file."""
    string = loader.construct_scalar(node)
    return getattr(app, string)


def yaml_dump(data, stream):
    """A version of yaml.dump with custom representers."""
    yaml.add_representer(unit.Quantity, quantity_representer)
    yaml.add_multi_representer(Singleton, singleton_representer)
    yaml.dump(data, stream)


def yaml_load(stream, Loader=yaml.SafeLoader):
    """A version of yaml.load with custom constructors."""
    yaml.add_constructor('!quantity', quantity_constructor, Loader=Loader)
    yaml.add_constructor('!openmm', singleton_constructor, Loader=Loader)
    return yaml.load(stream, Loader=Loader)

# ==========================================================================================
# Exceptions
# ==========================================================================================


class OpenMMSystemsException(Exception):
    """Base Exception for Package"""
    pass


# ==========================================================================================
# Paths
# ==========================================================================================

def get_data_file(relative_path):
    """Get the full path to one of the reference files.

    In the source distribution, these files are in ``openmmtools/_openmmtools_data/*/``,
    but on installation, they're moved to somewhere in the user's python
    site-packages directory.

    Parameters
    ----------
    name : str
        Name of the file to load (with respect to the repex folder).

    """
    fn = resource_filename('openmmsystems.systems', relative_path)

    if not os.path.exists(fn):
        raise ValueError("Sorry! %s does not exist. If you just added it, you'll have to re-install" % fn)

    return fn
