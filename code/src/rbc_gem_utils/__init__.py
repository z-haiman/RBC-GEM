"""Contains main functions and constants to facilitate working with the RBC-GEM reconstruction.

Currently houses all code associated with the RBC-GEM.
"""

from .annotation import get_annotation_df, get_id_annotation_mapping
from .database import DATABASE_DIRPATHS, check_database_release_online, check_release
from .io import read_cobra_model, write_cobra_model
from .qc import compare_tables
from .util import (
    ANALYSIS_PATH,
    ANNOTATION_PATH,
    COBRA_CONFIGURATION,
    CURATION_PATH,
    DATA_PATH,
    DATABASE_PATH,
    DEPRECATEDIDS_DIRPATH,
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
    check_if_valid,
    explode_column,
    show_versions,
    split_string,
)
from .visualization import visualize_comparison


def get_dirpath(*args, use_temp=None):
    if len(args) == 0:
        main_dir = ""
    if len(args) == 1:
        main_dir = args[0]
        sub_dir = None
    elif len(args) == 2:
        main_dir, sub_dir = args

    dirpath = ROOT_PATH

    temp_dirpath_dict = {
        "external": EXTERNAL_PATH,
        "raw": RAW_PATH,
        "interim": INTERIM_PATH,
        "processed": PROCESSED_PATH,
    }

    if use_temp is not None:
        try:
            dirpath /= temp_dirpath_dict[use_temp]
        except KeyError as e:
            raise KeyError(
                f"Unrecognized temporary directoruy: {e}. "
                f"Must be one of the following: {list(temp_dirpath_dict)}"
            )

    # Return path variable corresponding to model
    if main_dir == "model":
        dirpath /= MODEL_PATH
    # Return path variable corresponding to annotations
    elif main_dir == "annotation":
        dirpath /= ANNOTATION_PATH
    # Return path variable corresponding to annotations
    elif main_dir == "curation":
        dirpath /= CURATION_PATH
    elif main_dir == "analysis":
        dirpath /= ANALYSIS_PATH
    elif main_dir == "map":
        dirpath /= MAP_PATH
    # Return path variable corresponding to omics
    elif main_dir.endswith("omics"):
        valid = {"proteomics", "metabolomics", "genomics"}
        check_if_valid(
            main_dir,
            valid,
            msg=f"Must be one of the following {valid}. Invalid directory: ",
        )
        dirpath /= main_dir
        if sub_dir is not None:
            dirpath /= sub_dir

    elif main_dir == "database":
        dirpath /= DATABASE_PATH
        if sub_dir is not None:
            try:
                dirpath /= DATABASE_DIRPATHS[sub_dir]
            except KeyError as e:
                raise KeyError(
                    f"Unrecognized database: {e}. "
                    f"Must be one of the following: {list(DATABASE_DIRPATHS)}"
                )

    elif main_dir == "deprecatedIdentifiers":
        dirpath /= DEPRECATEDIDS_DIRPATH

    return dirpath


# FIXME Find a way to single source version with version.txt outside repository, only if tools are to be included as same version as model
__version__ = "0.0.3"
