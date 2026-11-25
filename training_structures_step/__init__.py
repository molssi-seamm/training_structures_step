# -*- coding: utf-8 -*-

"""
training_structures_step
A SEAMM plug-in for Training Structures
"""

# Bring up the classes so that they appear to be directly in
# the training_structures_step package.

from .training_structures import TrainingStructures  # noqa: F401, E501
from .training_structures_parameters import TrainingStructuresParameters  # noqa: F401, E501
from .training_structures_step import TrainingStructuresStep  # noqa: F401, E501
from .tk_training_structures import TkTrainingStructures  # noqa: F401, E501

from .metadata import metadata  # noqa: F401

# Handle versioneer
from ._version import get_versions

__author__ = "Paul Saxe"
__email__ = "psaxe@molssi.org"
versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions
