import logging
from pathlib import Path

import requests

from rbc_gem_utils.util import DATABASE_PATH, ROOT_PATH, check_if_valid


LOGGER = logging.getLogger(__name__)


REACTOME_URL = "https://reactome.org"
REACTOME_FILENAMES = {
    "ComplexParticipantsPubMedIdentifiers_human",
}
REACTOME_DB_TAG = "Reactome"
REACTOME_PATH = Path(REACTOME_DB_TAG)


def download_database_Reactome(filename=None, database_dirpath=None):
    """TODO DOCSTRING."""
    if filename is None:
        filename = REACTOME_FILENAMES
    else:
        filename = check_if_valid(filename, REACTOME_FILENAMES, "Invalid filenames: ")

    if database_dirpath is not None:
        # Ensure the path exists
        Path(database_dirpath).mkdir(parents=False, exist_ok=True)
    else:
        database_dirpath = ROOT_PATH / DATABASE_PATH / REACTOME_PATH

    for fname in filename:
        fileurl = f"{REACTOME_URL}/download/current/{fname}"
        response = requests.get(fileurl)
        response.raise_for_status()

        with open(database_dirpath / f"{fname}.tsv", "w") as file:
            file.write(response.text)

        LOGGER.info("`%s.tsv` saved at `%s`", fname, database_dirpath)
