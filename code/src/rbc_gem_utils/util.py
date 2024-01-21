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
DATA_PATH = "/data"
RAW_PATH = f"{DATA_PATH}/raw"
INTERIM_PATH = f"{DATA_PATH}/interim"
PROCESSED_PATH = f"{DATA_PATH}/processed"
EXTERNAL_PATH = f"{DATA_PATH}/external"
ANNOTATION_PATH = f"{DATA_PATH}/annotation"
CURATION_PATH = f"{DATA_PATH}/curation"
DATABASE_PATH = f"{EXTERNAL_PATH}/database"
RESULTS_PATH = f"{DATA_PATH}/results"


GEM_NAME = "RBC-GEM"
GEM_URL = f"{RAW_GH_URL}/z-haiman/{GEM_NAME}"
GEM_MODEL_FILETYPES = {"mat", "json", "xml", "yml"}


def show_versions():
    """Print dependency information."""
    print_dependencies("rbc_gem_utils")


def split_string(string_to_split, sep=None, raise_error=False):
    """Split a string using the given seperator, removing any white space on each side of the item.

    Identical duplicates are removed, however the order of items is preserved in the returned list.
    """
    if sep is None:
        sep = ";"
    # Let any ordering occur before this method, e.g., first item is active acession and others are secondary accessions
    # Using dict constructor removes duplicates and preserves order based on entry into dict.
    if isinstance(string_to_split, str):
        return list(dict.fromkeys([x.strip() for x in string_to_split.split(sep)]))

    elif raise_error:
        raise TypeError("Must be a string")
    else:
        return string_to_split


def build_string(set_of_components, sep=None, raise_error=False):
    """Build a string using the given seperator after sorting the given set of components.

    Identical duplicates are removed, however the order of items is preserved in the returned string.
    """
    # Using dict constructor removes duplicates and preserves order based on entry into dict.
    # Let any ordering occur before this
    if sep is None:
        sep = ";"
    if not isinstance(set_of_components, float):
        return sep.join(list(dict.fromkeys(ensure_iterable(set_of_components))))

    elif raise_error:
        raise TypeError("Cannot be a float")
    else:
        return set_of_components


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


def explode_column(df, name, sep=None):
    df = df.copy()
    name = ensure_iterable(name)
    for key in name:
        df[key] = df[key].apply(lambda x: split_string(x, sep=sep))
    return df.explode(name)


def has_value_type(element):
    return bool(element.text is not None and bool(element.text.strip()))
