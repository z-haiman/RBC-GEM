import logging
from pathlib import Path
from cobra import Configuration

from depinfo import print_dependencies

LOGGER = logging.getLogger(__name__)

# Constants
COBRA_CONFIGURATION = Configuration()
# TODO Probably a better way to determine the root directory, considering its the name of the project. 
# FIXME Right now this only works with an editable install
REPO_PATH = str(Path(__file__).resolve().parent.parent.parent.parent)
GEM_NAME = "RBC-GEM"


def show_versions():
    """Print dependency information."""
    print_dependencies("rbc_gem_utils")


def split_string(string_to_split, sep=";"): 
    """Split a string using the given seperator and remove any white space."""
    return set([x.strip() for x in string_to_split.split(sep)])


def build_string(set_of_components, sep=";"): 
    """Build a string using the given seperator after sorting the given set of components."""
    return sep.join(sorted(set(set_of_components)))