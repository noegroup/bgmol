"""
openmmsystems.py
Collection of OpenMM systems

Handles the primary functions
"""

import os

import click
from openmmsystems.api import list_openmmsystems, list_openmmtools_systems, list_toysystems


@click.group()
@click.version_option()
def main():
    pass


@main.command()
def systems():
    """Print a list of all available systems"""
    print("\n- ".join(
        ["Toy Systems\n-------------"] +
        list_toysystems()
    ))
    print("\n- ".join(
        ["\n\nOpenMMSystems\n-------------"] +
        list_openmmsystems()
    ))
    print("\n- ".join(
        ["\n\nopenmmtools\n-------------"] +
        list_openmmtools_systems()
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

