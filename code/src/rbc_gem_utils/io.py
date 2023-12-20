"""
Functions for importing and exporting the RBC-GEM using COBRApy from anywhere in the repo.
"""
import logging

LOGGER = logging.getLogger(__name__)

from cobra.io import (
    load_json_model, load_matlab_model, load_yaml_model, read_sbml_model,
    save_json_model, save_matlab_model, save_yaml_model, write_sbml_model,
)

try:
    import simplejson as json
except ImportError:
    import json


from .util import GEM_NAME, REPO_PATH

def write_rbc_model(model, filetype="xml", **kwargs):
    """Save the RBC-GEM as a COBRA model. The `filetype` determines which function to use.

    Parameters
    ----------
    model : cobra.Model
        The RBC-GEM model to represent.
    filetype : str, optional
        The type of model file. Default value is 'xml' for the SBML model file.
    **kwargs
        Passed to the underlying function utilized in saving/writing the model file.

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
    functions_dict = {
        "sbml": write_sbml_model,
        "xml": write_sbml_model,
        "mat": save_matlab_model,
        "json": save_json_model,
        "yml": save_yaml_model,
        "yaml": save_yaml_model,
    }
    try:
        save_function = functions_dict[filetype]
    except KeyError as e:
        raise ValueError(f"Unrecognized file type {e}.")
    else:
        save_function(model, f"{REPO_PATH}/model/{GEM_NAME}.{filetype}", **kwargs)
        LOGGER.info("Model `%s` saved as a `.%s` file", model.id, filetype)


def read_rbc_model(filetype='xml', **kwargs):
    """Load the RBC-GEM as a COBRA model. The `filetype` determines which function to use.

    Parameters
    ----------
    filetype : str, optional
        The type of model file. Default value is 'xml' for the SBML model file.
    **kwargs
        Passed to the underlying function utilized in loading/reading the model file.

    Returns
    -------
    cobra.Model
        The RBC-GEM as a COBRA model.

    See Also
    --------
    load_json_model
        Underlying function utilized for loading a COBRA model from a `.json` file.
    load_matlab_model
        Underlying function utilized for loading a COBRA model from a `.mat` file.
    load_yaml_model
        Underlying function utilized for loading a COBRA model from a `.yml` file or `.yaml` file.
    read_sbml_model
        Underlying function utilized for loading a COBRA model from a `.xml` file or `.sbml` file.

    """
    functions_dict = {
        "sbml": read_sbml_model,
        "xml": read_sbml_model,
        "mat": load_matlab_model,
        "json": load_json_model,
        "yml": load_yaml_model,
        "yaml": load_yaml_model,
    }
    try:
        load_function = functions_dict[filetype]
    except KeyError as e:
        raise ValueError(f"Unrecognized file type {e}.")
    else:
        model = load_function(f"{REPO_PATH}/model/{GEM_NAME}.{filetype}", **kwargs)
        LOGGER.info("Model `%s` loaded from an `.%s` file", model.id, filetype)

    return model
