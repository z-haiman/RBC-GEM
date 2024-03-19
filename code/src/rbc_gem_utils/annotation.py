import logging

import pandas as pd
from cobra import DictList

from rbc_gem_utils.util import build_string, check_if_valid, ensure_iterable


LOGGER = logging.getLogger(__name__)


def set_sbo_default_annotations(model, annotation_type, verbose=False):
    """Set default SBO annotations on the model."""
    annotation_type = check_if_valid(
        annotation_type,
        ["reactions", "metabolites", "genes"],
        "Unrecognized annotations:",
    )
    for attr in annotation_type:
        for model_object in getattr(model, attr):
            if attr == "reactions":
                if model_object.boundary and model_object in model.exchanges:
                    sbo_term = "SBO:0000627"  # Exchange reaction
                elif model_object.boundary and model_object in model.demands:
                    sbo_term = "SBO:0000628"  # Demand reaction
                elif model_object.boundary and model_object in model.sinks:
                    sbo_term = "SBO:0000632"  # Sink reaction
                elif len(model_object.compartments) > 1:
                    sbo_term = "SBO:0000185"  # Transport reaction
                else:
                    sbo_term = "SBO:0000176"  # Biochemical reaction
            elif attr == "metabolites":
                sbo_term = "SBO:0000247"  # Metabolite
            elif attr == "genes":
                sbo_term = "SBO:0000243"  # Gene
            else:
                continue

            try:
                assert model_object.annotation["sbo"] == sbo_term
            except KeyError:
                msg = f"SBO term set for {model_object.id}"
                LOGGER.info(msg)
                if verbose:
                    print(msg)
                model_object.annotation["sbo"] = sbo_term
            except AssertionError:
                msg = f"SBO term changed for {model_object.id}: {model_object.annotation['sbo']} --> {sbo_term}"
                LOGGER.info(msg)
                if verbose:
                    print(msg)
                model_object.annotation["sbo"] = sbo_term

    return model


def get_id_annotation_mapping(objects, annotation_key, default=None):
    """Return a dictionary containing model identifiers mapped to annotation values.

    Parameters
    ----------
    objects :
    annotation_key : str, iterable
    default :

    Returns
    -------

    """
    objects = DictList(objects)
    objects = objects.query(
        lambda x: x.annotation.get(annotation_key, default) is not None
    )
    objects.sort()
    return dict(
        zip(
            objects.list_attr("id"),
            [x.annotation.get(annotation_key, default) for x in objects],
        )
    )


def get_annotation_df(objects, annotation_key, sep=";"):
    """Return a table containing object identifiers mapped to their annotations.

    Parameters
    ----------
    objects :
    annotation_key : str, iterable

    Returns
    -------

    """
    # Ensure annotation key is iterable
    annotation_key = ensure_iterable(annotation_key)
    # Get annotation data for each key in the annotation keys.
    annotation_data = {}
    for key in annotation_key:
        # Get values and convert any lists to strings by joinining them together
        values = get_id_annotation_mapping(
            objects,
            key,
            default="",  # Default is None so entries without values are not returned
        )
        # No float NaN values can exist before this method (empty strings can)
        annotation_data[key] = {k: build_string(v, sep) for k, v in values.items()}

    # Create a DataFrame, index will be object identifiers
    df = pd.DataFrame.from_dict(annotation_data, orient="columns")
    # Sort index to match original objects, then reset index to make IDs a regular column
    df.index.name = "id"
    df = df.loc[objects.list_attr("id")].reset_index()
    # Replace empty strings with NaN
    df = df.replace("", float("nan"))
    return df
