from warnings import warn

from rbc_gem_utils.util import ensure_iterable

from .complexportal import (
    COMPLEXPORTAL_DB_TAG,
    COMPLEXPORTAL_PATH,
    COMPLEXPORTAL_VERSION_EXPECTED,
    get_version_ComplexPortal,
)
from .drugbank import (
    DRUGBANK_DB_TAG,
    DRUGBANK_PATH,
    DRUGBANK_VERSION_EXPECTED,
    get_version_DrugBank,
)
from .ec import EC_DB_TAG, EC_VERSION_EXPECTED, get_version_EC
from .metatlas import (
    HUMANGEM_DB_TAG,
    HUMANGEM_PATH,
    HUMANGEM_VERSION_EXPECTED,
    get_version_HumanGEM,
)
from .mim import MIM_DB_TAG, MIM_PATH, get_last_updated_dates_MIM
from .tcdb import TCDB_DB_TAG, TCDB_PATH
from .uniprot import (
    UNIPROT_DB_TAG,
    UNIPROT_PATH,
    UNIPROT_VERSION_EXPECTED,
    get_version_UniProt,
)


CDCDB_DB_TAG = "CDCDB"
CDCDB_PATH = "/CDCDB"
DRUGCENTRAL_DB_TAG = "DrugCentral"
DRUGCENTRAL_PATH = "/DrugCentral"


def check_database_version_online(database, expected=None, verbose=False):
    """Check the database version online against the expected version.

    Indicates whether the expected database version matches the version retrieved online.
    If False, the database may have changed since the last time the RBC-GEM codebase was updated.

    """
    database_versfunc_expected_dict = {
        COMPLEXPORTAL_DB_TAG: (
            get_version_ComplexPortal,
            COMPLEXPORTAL_VERSION_EXPECTED,
        ),
        DRUGBANK_DB_TAG: (get_version_DrugBank, DRUGBANK_VERSION_EXPECTED),
        EC_DB_TAG: (get_version_EC, EC_VERSION_EXPECTED),
        HUMANGEM_DB_TAG: (get_version_HumanGEM, HUMANGEM_VERSION_EXPECTED),
        UNIPROT_DB_TAG: (get_version_UniProt, UNIPROT_VERSION_EXPECTED),
    }
    try:
        get_version_func, default_expected = database_versfunc_expected_dict[database]
    except KeyError as e:
        raise KeyError(
            f"Unrecognized database: {e}\n"
            f"Must be one of the following: {list(database_versfunc_expected_dict)}"
        )
    else:
        if expected is None:
            expected = default_expected
    return check_version(get_version_func(), expected, verbose=verbose)


def check_version(current, expected, verbose=False):
    """Check the current version against the expected version, returning ``True`` if they match."""
    if current == expected:
        if verbose:
            print("Current and expected versions match.")
        return True
    else:
        if verbose:
            warn("Current and expected versions are not the same.")
            print(f"Current version: {current}.")
            print(f"Expected version: {expected}.")
        return False
