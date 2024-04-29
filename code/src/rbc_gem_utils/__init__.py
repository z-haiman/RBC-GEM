"""Contains main functions and constants to facilitate working with the RBC-GEM reconstruction.

Currently houses all code associated with the RBC-GEM.
"""

from .annotation import get_annotation_df, get_id_annotation_mapping
from .database import check_database_version_online, check_version
from .io import read_cobra_model, read_rbc_model, write_cobra_model, write_rbc_model
from .qc import compare_tables
from .util import (
    ANNOTATION_PATH,
    COBRA_CONFIGURATION,
    CURATION_PATH,
    DATABASE_PATH,
    EXTERNAL_PATH,
    GEM_NAME,
    GEM_URL,
    INTERIM_PATH,
    MAP_NAMES,
    MAP_PATH,
    MODEL_PATH,
    PARAMETERIZATION_PATH,
    PROCESSED_PATH,
    RAW_PATH,
    RESULTS_PATH,
    ROOT_PATH,
    build_string,
    explode_column,
    show_versions,
    split_string,
)
from .visualization import visualize_comparison


# FIXME Find a way to single source version with version.txt outside repository, only if tools are to be included as same version as model
__version__ = "0.0.1"
