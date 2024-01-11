"""
Contains various miscellaneous utility functions to help work with the RBC-GEM repository
"""
import logging
from pathlib import Path
from cobra import Configuration

from depinfo import print_dependencies

LOGGER = logging.getLogger(__name__)

# Constants
COBRA_CONFIGURATION = Configuration()
RAW_GH_URL = "https://raw.githubusercontent.com"
# TODO Probably a better way to determine the root directory, considering its the name of the project. 
# FIXME Right now this only works with an editable install './RBC-GEM'
ROOT_PATH = str(Path(__file__).resolve().parent.parent.parent.parent)
# Relevant relative paths
MODEL_PATH = "/model"
RAW_PATH = "/data/raw"
INTERIM_PATH = "/data/interim"
PROCESSED_PATH = "/data/processed"
EXTERNAL_PATH = "/data/external"
ANNOTATION_PATH = "/data/annotation"
DATABASE_PATH = f"{EXTERNAL_PATH}/database"


GEM_NAME = "RBC-GEM"


def show_versions():
    """Print dependency information."""
    print_dependencies("rbc_gem_utils")


def split_string(string_to_split, sep=";"): 
    """Split a string using the given seperator, removing any white space on each side of the item.
    
    Identical duplicates are removed, however the order of items is preserved in the returned list.
    """
    # Let any ordering occur before this method, e.g., first item is active acession and others are secondary accessions
    # Using dict constructor removes duplicates and preserves order based on entry into dict.
    return list(dict.fromkeys([x.strip() for x in string_to_split.split(sep)]))


def build_string(set_of_components, sep=";"): 
    """Build a string using the given seperator after sorting the given set of components.
    
    Identical duplicates are removed, however the order of items is preserved in the returned string.
    """
    # Using dict constructor removes duplicates and preserves order based on entry into dict.
    # Let any ordering occur before this
    return sep.join(list(dict.fromkeys(ensure_iterable(set_of_components))))

def ensure_iterable(to_check):
    """Ensure the item is an iterable before returning."""
    if isinstance(to_check, str) or not hasattr(to_check, "__iter__"):
        return [to_check]
    
    return to_check 


def check_if_valid(to_check, valid_values, msg):
    """Check if provided values are valid and raise a ValueError for unrecongized values.
    
    Parameters
    ----------
    to_check : iterable, str
        The value(s) to check. 
    valid_values : iterable, str
        The value(s) that are considered to be valid. 

    Returns
    -------
    iterable
        Valid values that were checked. If the original data was a string, it will be returned within a list.
    """
    to_check = ensure_iterable(to_check)
    valid_values = ensure_iterable(valid_values)
    invalid = set(to_check).difference(valid_values)
    if invalid:
        raise ValueError(f"{msg}: {invalid}")
    
    return to_check


def strip_plural(string):
    if string.endswith("ies"):
        return f'{string.rstrip("ies")}y'
    elif string.endswith("s"):
        return string.rstrip("s")
    else:
        return string
