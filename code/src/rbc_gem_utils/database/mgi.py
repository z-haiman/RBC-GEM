"""Functions to extract relevant information from the MGI database.

Notes
-----
* Main site: https://www.informatics.jax.org/
* Code written or updated based on MGI 6.24, released 05/06/2025
* Code last updated: May 2025
"""

import logging
from datetime import datetime
from pathlib import Path

import requests
from bs4 import BeautifulSoup

from rbc_gem_utils.util import DATABASE_PATH, ROOT_PATH, check_if_valid


LOGGER = logging.getLogger(__name__)

MGI_URL = "https://www.informatics.jax.org"
MGI_FILENAMES = {
    "HGNC_AllianceHomology.rpt",
    "MRK_SwissProt.rpt",
    "MRK_SwissProt_TrEMBL.rpt",
    "MRK_Sequence.rpt",
}
MGI_RELEASE_EXPECTED = "MGI 6.24, 08/05/2025"
MGI_DB_TAG = "MGI"
MGI_PATH = Path(MGI_DB_TAG)


def download_database_MGI(filename=None, database_dirpath=None):
    """TODO DOCSTRING."""
    if filename is None:
        filename = MGI_FILENAMES
    else:
        filename = check_if_valid(filename, MGI_FILENAMES, "Invalid filenames: ")

    if database_dirpath is not None:
        # Ensure the path exists
        Path(database_dirpath).mkdir(parents=False, exist_ok=True)
    else:
        database_dirpath = ROOT_PATH / DATABASE_PATH / MGI_PATH

    for fname in filename:
        fileurl = f"{MGI_URL}/downloads/reports/{fname}"
        # TODO handle non-.rpt files
        suffix = Path(fileurl).suffix
        if suffix == ".rpt":
            response = requests.get(fileurl)
            response.raise_for_status()
            with open(database_dirpath / fname.replace(suffix, ".tsv"), "w") as file:
                file.write(response.text)
        else:
            NotImplementedError("Handling of {suffix} files not yet implemented.")
        LOGGER.info("`%s.tsv` saved at `%s`", fname, database_dirpath)


def get_release_MGI():
    """Return the current version of MGI - Mouse Genomics Informatics database."""
    response = requests.get(f"{MGI_URL}/mgihome/projects/aboutmgi.shtml")
    response.raise_for_status()
    soup = BeautifulSoup(response.text, "html.parser")
    items = soup.find(
        "td", attrs=dict(align="center", style="font-size:9px;", width="34%")
    )
    items = [x.strip() for x in items.text.split("\n") if x.strip()]
    for item in items:
        try:
            item = datetime.strptime(item, "%m/%d/%Y")
        except ValueError:
            if "MGI " in item:
                release = item
            else:
                continue
        else:
            date = item.strftime("%m/%d/%Y")
    return f"{release}, {date}"
