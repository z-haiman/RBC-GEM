"""Contains main functions and constants to facilitate working with the RBC-GEM reconstruction."""

from .annotation import (
    get_annotation_df, get_id_annotation_mapping
)
from .io import (
    read_rbc_model, write_rbc_model, 
    read_cobra_model, write_cobra_model
)
from .qc import (
    compare_tables
)
from .util import (
    ROOT_PATH, 
    MODEL_PATH, 
    DATABASE_PATH, 
    ANNOTATION_PATH,
    INTERIM_PATH,
    GEM_NAME, 
    COBRA_CONFIGURATION, 
    show_versions,
)
from .visualization import (
    visualize_comparison
)

# FIXME Find a way to single source version with version.txt outside repository, only if tools are to be included as same version as model
__version__ = "0.0.1"

