"""
YAML utilities
"""

import yaml
from mdtraj.utils.unit import _str_to_unit

from bgmol.util.importing import import_openmm
mm, unit, app = import_openmm()


__all__ = [
    "parse_quantity",
    "quantity_representer",
    "quantity_constructor",
    "singleton_constructor",
    "singleton_representer",
    "yaml_load",
    "yaml_dump"
]


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
    if isinstance(object, app.internal.singleton.Singleton)
]


def singleton_representer(dumper, data):
    """Converts openmm.app.internal.Singleton object into a yaml string."""
    return dumper.represent_scalar('!openmm', data.__class__.__name__)


def singleton_constructor(loader, node):
    """Parser for openmm.app.internal.Singleton objects from yaml file."""
    string = loader.construct_scalar(node)
    return getattr(app, string)


def yaml_dump(data, stream):
    """A version of yaml.dump with custom representers."""
    yaml.add_representer(unit.Quantity, quantity_representer)
    yaml.add_multi_representer(app.internal.singleton.Singleton, singleton_representer)
    yaml.dump(data, stream)


def yaml_load(stream, Loader=yaml.SafeLoader):
    """A version of yaml.load with custom constructors."""
    yaml.add_constructor('!quantity', quantity_constructor, Loader=Loader)
    yaml.add_constructor('!openmm', singleton_constructor, Loader=Loader)
    return yaml.load(stream, Loader=Loader)
