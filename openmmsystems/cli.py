"""
openmmsystems.py
Collection of OpenMM systems

Handles the primary functions
"""

import os

import click
from openmmsystems.api import get_openmmsystems_names, get_openmmtools_system_names, get_toysystem_names


@click.group()
@click.version_option()
def main():
    pass


@main.command()
def systems():
    """Print a list of all available systems"""
    print("\n- ".join(
        ["Toy Systems\n-------------"] +
        get_toysystem_names()
    ))
    print("\n- ".join(
        ["\n\nOpenMMSystems\n-------------"] +
        get_openmmsystems_names()
    ))
    print("\n- ".join(
        ["\n\nopenmmtools\n-------------"] +
        get_openmmtools_system_names()
    ))


@main.command()
def data():
    """Print a list of all available data sets"""
    pass


@main.command()
@click.argument("identifier")
def show(identifier):
    """Show details about a data set or system"""
    print("Show")

