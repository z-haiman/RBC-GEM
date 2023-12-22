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
