"""
Utility module
"""

from simtk import unit
import yaml


def quantity_representer(dumper, data):
    """Converts quantity object into a yaml string."""
    if data.unit == unit.dimensionless:
        return dumper.represent_scalar('!quantity', f'{data._value}')
    else:
        operator = '' if data.unit._name.startswith('/') else '*'
        u = data.unit._name.replace(' ','_')
        return dumper.represent_scalar('!quantity', f'{data._value} {operator} {u}')


def quantity_constructor(loader, node):
    """Parser for quantity objects from yaml file. Works only for float values."""
    string = loader.construct_scalar(node)
    split = string.replace(' ','').split('/')
    numerator_units = split[0].split('*')
    denominator_units = [] if len(split) == 1 else split[1].replace('(','').replace(')','').split('*')
    v = float(numerator_units[0])
    u = unit.dimensionless
    for nu in numerator_units[1:]:
        assert hasattr(unit, nu)
        u = u * getattr(unit, nu)
    for du in denominator_units:
        u = u / getattr(unit, du)
    return unit.Quantity(v, u)


def yaml_dump_with_quantity(data, stream):
    """A version of yaml.dump with quantity objects."""
    yaml.add_representer(unit.Quantity, quantity_representer)
    yaml.dump(data, stream)


def yaml_load_with_quantity(stream, Loader=yaml.SafeLoader):
    """A version of yaml.load with quantity objects."""
    yaml.add_constructor('!quantity', quantity_constructor, Loader=Loader)
    return yaml.load(stream, Loader=Loader)
