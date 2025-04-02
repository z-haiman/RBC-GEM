"""Functions to extract relevant information from ENZYME

Notes
-----
* Main site: https://enzyme.expasy.org/
* Code written based on release of 27-Mar-24 (6759 active entries)
* Code last updated: May 2024

"""

import logging
import re
from pathlib import Path
from warnings import warn

import requests

from rbc_gem_utils.util import DATABASE_PATH, ROOT_PATH, check_if_valid


LOGGER = logging.getLogger(__name__)

EC_URL = "https://ftp.expasy.org/databases/enzyme/"
EC_RELEASE_RE = re.compile(r"Release:     (?P<release>\d+-\S+-\d+)")
EC_FILENAMES = {
    "enzyme.dat",
    "enzyme.rdf",
    "enzclass.txt",
    "enzuser.txt",
    "enzyme.get",
}
EC_DB_TAG = "EC"
EC_PATH = Path(EC_DB_TAG)
EC_VERSION_EXPECTED = "27-Nov-2024"


def get_version_EC():
    """Return the current version of ENZYME nomenclature database."""
    response = requests.get(f"{EC_URL}/enzclass.txt")
    response.raise_for_status()

    match = EC_RELEASE_RE.search(response.text)
    if match is None:
        warn("Could not retrieve version")
        return None

    return match.group("release")


def download_database_EC(filename=None, database_dirpath=None):
    """TODO DOCSTRING."""
    if filename is None:
        filename = EC_FILENAMES
    else:
        filename = check_if_valid(filename, EC_FILENAMES, "Invalid filenames: ")

    if database_dirpath is not None:
        # Ensure the path exists
        Path(database_dirpath).mkdir(parents=False, exist_ok=True)
    else:
        database_dirpath = ROOT_PATH / DATABASE_PATH / EC_PATH

    for fname in filename:
        fileurl = f"{EC_URL}{fname}"
        response = requests.get(fileurl)
        response.raise_for_status()

        fname = f"ec_{fname}"
        with open(database_dirpath / fname, "w") as file:
            file.write(response.text)

        LOGGER.info("`%s` saved at `%s`", fname, database_dirpath)
