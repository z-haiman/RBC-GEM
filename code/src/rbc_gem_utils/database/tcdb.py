import logging
from pathlib import Path

import requests

from rbc_gem_utils.util import DATABASE_PATH, ROOT_PATH, check_if_valid


LOGGER = logging.getLogger(__name__)


TCDB_URL = "https://www.tcdb.org"
TCDB_FILENAMES = {
    "getSubstrates",
    "families",
    "refseq",
    "listSuperfamilies",
    "acc2tcid",
    "go",
    "pdb",
    "pfam",
}
TCDB_MIRIAMS = ["tcdb", "uniprot", "pfam", "refseq", "pdb", "go", "chebi"]
TCDB_DB_TAG = "TCDB"
TCDB_PATH = Path(TCDB_DB_TAG)


def download_database_TCDB(filename=None, database_dirpath=None):
    """TODO DOCSTRING."""
    if filename is None:
        filename = TCDB_FILENAMES
    else:
        filename = check_if_valid(filename, TCDB_FILENAMES, "Invalid filenames: ")

    if database_dirpath is not None:
        # Ensure the path exists
        Path(database_dirpath).mkdir(parents=False, exist_ok=True)
    else:
        database_dirpath = ROOT_PATH / DATABASE_PATH / TCDB_PATH

    for fname in filename:
        if fname in {"getSubstrates", "listSuperfamilies"}:
            fileurl = f"{TCDB_URL}/cgi-bin/substrates/{fname}.py"
        else:
            fileurl = f"{TCDB_URL}/cgi-bin/projectv/public/{fname}.py"
        response = requests.get(fileurl)
        response.raise_for_status()

        with open(database_dirpath / f"{fname}.tsv", "w") as file:
            file.write(response.text)

        LOGGER.info("`%s.tsv` saved at `%s`", fname, database_dirpath)
