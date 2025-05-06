from pathlib import Path
from warnings import warn

from rbc_gem_utils.util import ANNOTATION_PATH, DATABASE_PATH, INTERIM_PATH, ROOT_PATH

from .complexportal import (
    COMPLEXPORTAL_DB_TAG,
    COMPLEXPORTAL_PATH,
    COMPLEXPORTAL_RELEASE_EXPECTED,
    get_release_ComplexPortal,
)
from .drugbank import (
    DRUGBANK_DB_TAG,
    DRUGBANK_PATH,
    DRUGBANK_RELEASE_EXPECTED,
    get_release_DrugBank,
)
from .ec import EC_DB_TAG, EC_PATH, EC_RELEASE_EXPECTED, get_release_EC
from .metatlas import (
    HUMANGEM_DB_TAG,
    HUMANGEM_PATH,
    HUMANGEM_RELEASE_EXPECTED,
    get_release_HumanGEM,
)
from .mgi import MGI_DB_TAG, MGI_PATH, MGI_RELEASE_EXPECTED, get_release_MGI
from .mim import MIM_DB_TAG, MIM_PATH, get_last_updated_dates_MIM
from .tcdb import TCDB_DB_TAG, TCDB_PATH
from .uniprot import (
    UNIPROT_DB_TAG,
    UNIPROT_PATH,
    UNIPROT_RELEASE_EXPECTED,
    get_release_UniProt,
)


CDCDB_DB_TAG = "CDCDB"
CDCDB_PATH = Path(CDCDB_DB_TAG)
DRUGCENTRAL_DB_TAG = "DrugCentral"
DRUGCENTRAL_PATH = Path(DRUGCENTRAL_DB_TAG)


def get_database_dirpath(database, use_interim=False):
    """Return the path variable corresponding to the database."""
    database_release_paths = {
        CDCDB_DB_TAG: CDCDB_PATH,
        COMPLEXPORTAL_DB_TAG: COMPLEXPORTAL_PATH,
        DRUGBANK_DB_TAG: DRUGBANK_PATH,
        DRUGCENTRAL_DB_TAG: DRUGCENTRAL_PATH,
        EC_DB_TAG: EC_PATH,
        HUMANGEM_DB_TAG: HUMANGEM_PATH,
        MGI_DB_TAG: MGI_PATH,
        MIM_DB_TAG: MIM_PATH,
        TCDB_DB_TAG: TCDB_PATH,
        UNIPROT_DB_TAG: UNIPROT_PATH,
    }
    try:
        database_path = database_release_paths[database]
    except KeyError as e:
        raise KeyError(
            f"Unrecognized database: {e}. "
            f"Must be one of the following: {list(database_release_paths)}"
        )
    if use_interim:
        return ROOT_PATH / INTERIM_PATH / database_path
    else:
        return ROOT_PATH / DATABASE_PATH / database_path


def get_annotation_dirpath(use_interim=False):
    """Return the path variable corresponding to the annotations."""
    if use_interim:
        return ROOT_PATH / INTERIM_PATH / ANNOTATION_PATH
    else:
        return ROOT_PATH / ANNOTATION_PATH


def check_database_release_online(database, expected=None, verbose=False, **kwargs):
    """Check the database release online against the expected release.

    Indicates whether the expected database release matches the release retrieved online.
    If False, the database may have changed since the last time the RBC-GEM codebase was updated.

    """
    database_release_functions = {
        COMPLEXPORTAL_DB_TAG: (
            get_release_ComplexPortal,
            COMPLEXPORTAL_RELEASE_EXPECTED,
        ),
        DRUGBANK_DB_TAG: (get_release_DrugBank, DRUGBANK_RELEASE_EXPECTED),
        EC_DB_TAG: (get_release_EC, EC_RELEASE_EXPECTED),
        HUMANGEM_DB_TAG: (get_release_HumanGEM, HUMANGEM_RELEASE_EXPECTED),
        MGI_DB_TAG: (get_release_MGI, MGI_RELEASE_EXPECTED),
        UNIPROT_DB_TAG: (get_release_UniProt, UNIPROT_RELEASE_EXPECTED),
    }
    try:
        get_release_func, default_expected = database_release_functions[database]
    except KeyError as e:
        raise KeyError(
            f"Unrecognized database: {e}. "
            f"Must be one of the following: {list(database_release_functions)}"
        )
    else:
        if expected is None:
            expected = default_expected
    if kwargs:
        return check_release(get_release_func(**kwargs), expected, verbose=verbose)
    else:
        return check_release(get_release_func(), expected, verbose=verbose)


def check_release(current, expected, verbose=False):
    """Check the current release against the expected release, returning ``True`` if they match."""
    if current == expected:
        if verbose:
            print(f"Current and expected releases match. Current release: {current}")
        return True
    else:
        if verbose:
            warn("Current and expected releases are not the same.")
            print(f"Current release: {current}.")
            print(f"Expected release: {expected}.")
        return False
