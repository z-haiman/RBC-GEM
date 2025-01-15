"""
Contains various miscellaneous utility functions to help work with the RBC-GEM repository
"""

import logging
from pathlib import Path

import numpy as np
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
MAP_PATH = "/map"
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
PARAMETERIZATION_PATH = f"{DATA_PATH}/parameterization"

GEM_NAME = "RBC-GEM"
GEM_URL = f"{RAW_GH_URL}/z-haiman/{GEM_NAME}"
GEM_MODEL_FILETYPES = {"mat", "json", "xml", "yml"}
MAP_NAMES = {f"{GEM_NAME}.full.map"}


AVOGADRO_NUMBER = 6.02214076e23
DEFAULT_DRY_MASS_PER_CELL = 30  # pg / living RBC
# 31.4 pg/cell         PMID: 14919571
# 27 to 32 pg/cell     PMID: 23005682
# 27.3 +/- 5.3 pg/cell PMID: 22934287
# 33.4 +/- 4.1 pg/cell PMID: 26720876
DEFAULT_VOLUME_PER_CELL = 90  # fL / living RBC, fL = um^3
# bionumbers:101722;bionumbers:101723;bionumers:101724
DEFAULT_WATER_PER_VOLUME = 0.717
# 87 +/- 7 fL/cell PMID: 21250103
# 100.6 +/- 4 fL/cell PMID: 26720876


def show_versions():
    """Print dependency information."""
    print_dependencies("rbc_gem_utils")


def split_string(string_to_split, sep=";", raise_error=False):
    """Split a string using the given seperator, removing any white space on each side of the item.

    Identical duplicates are removed, however the order of items is preserved in the returned list.
    """
    # Let any ordering occur before this method, e.g., first item is active acession and others are secondary accessions
    # Using dict constructor removes duplicates and preserves order based on entry into dict.
    if isinstance(string_to_split, str):
        return list(dict.fromkeys([x.strip() for x in string_to_split.split(sep)]))

    elif raise_error:
        raise TypeError("Must be a string")
    else:
        return string_to_split


def build_string(set_of_components, sep=";", raise_error=False):
    """Build a string using the given seperator after sorting the given set of components.

    Identical duplicates are removed, however the order of items is preserved in the returned string.
    """
    # Using dict constructor removes duplicates and preserves order based on entry into dict.
    # Let any ordering occur before this
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


def strip_plural(s):
    """TODO DOCSTRING."""
    if s.endswith("ies"):
        return f"{s[:-3]}y"
    elif s.endswith("xes"):
        return f"{s[:-2]}"
    elif s.endswith("s"):
        return f"{s[:-1]}"
    else:
        return s


def explode_column(df, name, sep=";"):
    """TODO DOCSTRING."""
    df = df.copy()
    name = ensure_iterable(name)
    for key in name:
        df[key] = df[key].apply(lambda x: split_string(x, sep=sep))
    return df.explode(name)


def has_value_type(element):
    """TODO DOCSTRING."""
    return bool(element.text is not None and bool(element.text.strip()))


def format_summary(header, body):
    """Format the summary for pretty printing."""
    LOGGER.debug("Formatting summary")
    max_len = max([len(line) for line in body.split("\n")] + [len(header)])
    header = "".join(
        (
            " " * int((max_len - len(header)) / 2),
            header,
            " " * int((max_len - len(header)) / 2),
        )
    )
    summary = "\n".join(("=" * max_len, header, "-" * max_len, body, "=" * max_len))
    return summary


def log_msg(logger, lvl, msg, *args, print_lvl=0):
    """Log the message at the request level, and print to console if desired."""
    logger.log(lvl, msg, *args)
    if print_lvl and 0 <= lvl - print_lvl:
        print(msg % args)


def convert_gDW_to_L(
    value,
    pgDW_per_cell=DEFAULT_DRY_MASS_PER_CELL,
    fL_per_cell=DEFAULT_VOLUME_PER_CELL,
    water_fraction=DEFAULT_WATER_PER_VOLUME,
):
    """Convert the value given in per dry weight to per liter"""
    # Conversion from Aq. volume
    # (pgDW / (fL RBC * (H2O / RBC))) --> (pgDW / 0.001 pL) --> 1000 * (gDW / L) conversion factor to value per L
    # Assuming value is given as mmol / gDW, then new value is returned as mmol / L
    return value * 1000 * (pgDW_per_cell / (fL_per_cell * water_fraction))


def convert_L_to_gDW(
    value,
    pgDW_per_cell=DEFAULT_DRY_MASS_PER_CELL,
    fL_per_cell=DEFAULT_VOLUME_PER_CELL,
    water_fraction=DEFAULT_WATER_PER_VOLUME,
):
    """Convert the value given in per liter to per dry weight"""
    # (fL RBC * (H2O / RBC)/ pgDW) --> (0.001 pL / pgDW) --> 0.001 * (L / gDW) conversion factor to value per L
    # Assuming value is given as mmol / L, then new value is returned as mmol / gDW
    return value * 0.001 * ((fL_per_cell * water_fraction) / pgDW_per_cell)
