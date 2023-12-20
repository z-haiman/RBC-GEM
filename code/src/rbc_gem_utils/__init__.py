"""Contains main functions and constants to facilitate working with the RBC-GEM reconstruction."""

from .io import read_rbc_model, write_rbc_model
from .util import (
    REPO_PATH, GEM_NAME, COBRA_CONFIGURATION, show_versions
)


# FIXME Find a way to single source version with version.txt outside repository, only if tools are to be included as same version as model
__version__ = "0.0.1"

