import logging
import pathlib

import requests

from rbc_gem_utils.util import DATABASE_PATH, ROOT_PATH, check_if_valid


LOGGER = logging.getLogger(__name__)


REACTOME_URL = "https://reactome.org"
REACTOME_FILENAMES = {
    "ComplexParticipantsPubMedIdentifiers_human",
}
REACTOME_PATH = "/Reactome"
REACTOME_DB_TAG = "Reactome"


def download_database_Reactome(filename=None, database_dirpath=None):
    """TODO DOCSTRING."""
    if filename is None:
        filename = REACTOME_FILENAMES
    else:
        filename = check_if_valid(filename, REACTOME_FILENAMES, "Invalid filenames: ")

    if database_dirpath is not None:
        # Ensure the path exists
        pathlib.Path(f"{database_dirpath}").mkdir(parents=False, exist_ok=True)
    else:
        database_dirpath = f"{ROOT_PATH}{DATABASE_PATH}{REACTOME_PATH}"

    for fname in filename:
        fileurl = f"{REACTOME_URL}/download/current/{fname}"
        response = requests.get(fileurl)
        response.raise_for_status()

        with open(f"{database_dirpath}/{fname}.tsv", "w") as file:
            file.write(response.text)

        LOGGER.info("`%s.tsv` saved at `%s`", fname, database_dirpath)
