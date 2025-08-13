"""
Functions for importing and exporting the RBC-GEM using COBRApy from anywhere in the repo.
"""

import io
import logging
from pathlib import Path
from warnings import warn


LOGGER = logging.getLogger(__name__)

from cobra.io import (
    load_json_model,
    load_matlab_model,
    load_yaml_model,
    read_sbml_model,
    save_json_model,
    save_matlab_model,
    save_yaml_model,
    write_sbml_model,
)

from .util import GEM_NAME, MODEL_PATH, ROOT_PATH


IO_FUNCTIONS_DICT = {
    "xml": {"read": read_sbml_model, "write": write_sbml_model},
    "mat": {"read": load_matlab_model, "write": save_matlab_model},
    "json": {"read": load_json_model, "write": save_json_model},
    "yml": {"read": load_yaml_model, "write": save_yaml_model},
}

ALTERNATE_FILETYPES_DICT = {
    "sbml": "xml",
    "yaml": "yml",
}


def write_cobra_model(model, filename, filetype=None, **kwargs):
    """Save/write the COBRA model file. Determines the function based on the file extension.

    Parameters
    ----------
    model : cobra.Model
        The COBRA model to represent.
    filename : str or file-like
        File path or descriptor to which the model is written.
    filetype : str
        Model filetype. Required for BytesIO and StringIO objects, otherwise attempts automatic determination from file extension.
    **kwargs
        Passed to the function utilized in saving the model.

    Returns
    -------
    cobra.Model
        The loaded COBRA model.

    See Also
    --------
    save_json_model
        Underlying function utilized for saving a COBRA model as a `.json` file.
    save_matlab_model
        Underlying function utilized for saving a COBRA model as a `.mat` file.
    save_yaml_model
        Underlying function utilized for saving a COBRA model as a `.yml` file or `.yaml` file.
    write_sbml_model
        Underlying function utilized for saving a COBRA model as a `.xml` file or `.sbml` file.

    """

    if isinstance(filename, (io.BytesIO, io.StringIO)) and filetype is None:
        raise ValueError("The `filetype` must be defined for IO objects.")
    elif not isinstance(filename, (io.BytesIO, io.StringIO)) and filetype is None:
        fname = filename.name if hasattr(filename, "name") else str(filename)
        filetype = fname.rsplit(".")[-1]
        if filetype in {"zip", "bz", "gz"}:
            filetype = fname.rsplit(".")[-2]

    filetype = ALTERNATE_FILETYPES_DICT.get(filetype, filetype)
    try:
        write_function = IO_FUNCTIONS_DICT[filetype]["write"]
    except KeyError as e:
        raise ValueError(f"Unrecognized file type {e}.")

    write_function(model, filename, **kwargs)


def read_cobra_model(filename, filetype=None, **kwargs):
    """Load/read the COBRA model file. Determines the function based on the file extension.

    Parameters
    ----------
    filename : str or file-like
        File path or descriptor that contains the model.
    filetype : str
        Model filetype. Required for BytesIO and StringIO objects, otherwise attempts automatic determination from file extension.
    **kwargs
        Passed to the function utilized in loading the model.

    Returns
    -------
    cobra.Model
        The loaded COBRA model.

    See Also
    --------
    load_json_model
        Underlying function utilized for loading a COBRA model from a `.json` file.
    load_matlab_model
        Underlying function utilized for loading a COBRA model from a `.mat` file.
    read_sbml_model
        Underlying function utilized for loading a COBRA model from a `.xml` file or `.sbml` file.
    load_yaml_model
        Underlying function utilized for loading a COBRA model from a `.yml` file or `.yaml` file.
    """
    if isinstance(filename, io.StringIO) and filetype is None:
        raise ValueError("The `filetype` must be defined for IO objects.")
    elif not isinstance(filename, io.StringIO) and filetype is None:
        fname = filename.name if hasattr(filename, "name") else str(filename)
        filetype = fname.rsplit(".")[-1]
        if filetype in {"zip", "bz", "gz"}:
            filetype = fname.rsplit(".")[-2]

    filetype = ALTERNATE_FILETYPES_DICT.get(filetype, filetype)
    try:
        read_function = IO_FUNCTIONS_DICT[filetype]["read"]
    except KeyError as e:
        raise ValueError(f"Unrecognized file type {e}.")

    model = read_function(filename, **kwargs)

    return model


# def write_rbc_model(model, filetype="xml", directory=None, **kwargs):
#     """Save/write the RBC-GEM as a COBRA model. The `filetype` determines which function to use.

#     Parameters
#     ----------
#     model : cobra.Model
#         The RBC-GEM model to represent.
#     filetype : str, iterable, optional
#         The type of model file. Default value is 'xml' for the SBML model file.
#         If ``filetype="all"``, all possible functions will be used.
#         Alternatively, an iterable of strings representing the filetypes can be provided.
#     directory : str, optional
#         The directory to save the model file.
#         If ``None``, defaults to the {ROOT_PATH}{MODEL_PATH}.
#     **kwargs
#         Passed to the underlying function utilized in saving/writing the model file.
#         If ``filetype="all"`` or an iterable is provided to `filetype`, a `dict` should be passed with
#         the keys representing filetypes and values representing dictionaries of ``kwargs``
#         for the `filetype` passed.

#     See Also
#     --------
#     save_json_model
#         Underlying function utilized for saving a COBRA model as a `.json` file.
#     save_matlab_model
#         Underlying function utilized for saving a COBRA model as a `.mat` file.
#     save_yaml_model
#         Underlying function utilized for saving a COBRA model as a `.yml` file or `.yaml` file.
#     write_sbml_model
#         Underlying function utilized for saving a COBRA model as a `.xml` file or `.sbml` file.

#     """
#     if isinstance(filetype, str):
#         if filetype == "all":
#             # Ensure no duplicates when it comes to SBML or yaml files.
#             filetype = sorted(
#                 set(IO_FUNCTIONS_DICT).difference(set(ALTERNATE_FILETYPES_DICT))
#             )
#         else:
#             kwargs = {filetype: kwargs}
#             filetype = [filetype]
#     if not directory:
#         directory = ROOT_PATH / MODEL_PATH
#     for ftype in set(filetype):
#         ftype_kwargs = kwargs.get(ftype, {})
#         if ftype in ALTERNATE_FILETYPES_DICT:
#             alt = ALTERNATE_FILETYPES_DICT[ftype]
#             msg = f'Extension "{ftype}" used instead of "{alt}"'
#             if ftype not in filetype:
#                 warn(f"{msg}, correcting file format for export")
#                 ftype = alt
#                 ftype_kwargs = kwargs.get(alt, {})
#             else:
#                 # No need to duplicate the export
#                 warn(f"{msg}, skipping to prevent duplicate export")
#                 continue

#         write_cobra_model(
#             model, Path(directory) / f"{GEM_NAME}.{ftype}", **ftype_kwargs
#         )
#         LOGGER.info("Model `%s` saved as a `.%s` file", model.id, ftype)


# def read_rbc_model(filetype="xml", directory=None, **kwargs):
#     """Load/read the RBC-GEM as a COBRA model. The `filetype` determines which function to use.

#     Parameters
#     ----------
#     filetype : str, optional
#         The type of model file. Default value is 'xml' for the SBML model file.
#     directory : str, optional
#         The directory where the model file is kept.
#         If ``None``, defaults to the {ROOT_PATH}{MODEL_PATH}.
#     **kwargs
#         Passed to the underlying function utilized in loading/reading the model file.

#     Returns
#     -------
#     cobra.Model
#         The RBC-GEM as a COBRA model.

#     See Also
#     --------
#     load_json_model
#         Underlying function utilized for loading a COBRA model from a `.json` file.
#     load_matlab_model
#         Underlying function utilized for loading a COBRA model from a `.mat` file.
#     load_yaml_model
#         Underlying function utilized for loading a COBRA model from a `.yml` file or `.yaml` file.
#     read_sbml_model
#         Underlying function utilized for loading a COBRA model from a `.xml` file or `.sbml` file.

#     """
#     if not directory:
#         directory = ROOT_PATH / MODEL_PATH
#     model = read_cobra_model(Path(directory) / f"{GEM_NAME}.{filetype}", **kwargs)
#     LOGGER.info("Model `%s` loaded from a `.%s` file", model.id, filetype)
#     return model
