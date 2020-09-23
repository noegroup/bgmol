"""
ProjectName
Collection of OpenMM systems
"""


# package-wide paths
import os
SYSTEMPATH = os.path.join(os.path.dirname(os.path.normpath(__file__)), "systems")
#SAMPLEPATH = os.path.join(os.path.dirname(os.path.normpath(__file__)), "datasets")
DATAPATH = os.path.join(os.path.dirname(os.path.normpath(__file__)), "data")


# Add imports here
from .api import *
from . import datasets
from . import systems

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions

