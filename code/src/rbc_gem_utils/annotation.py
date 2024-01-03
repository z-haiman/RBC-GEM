import pandas as pd
from cobra import DictList

from rbc_gem_utils.util import ensure_iterable, build_string

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
    objects = objects.query(lambda x: x.annotation.get(annotation_key, default) is not None)
    objects.sort()
    return dict(zip(objects.list_attr("id"), [x.annotation.get(annotation_key, default) for x in objects]))


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
            default="" # Default is None so entries without values are not returned
        )
        # No float NaN values can exist before this method (empty strings can)
        annotation_data[key] = {k: build_string(v, sep) for k, v in values.items()}

    # Create a DataFrame, index will be object identifiers
    df = pd.DataFrame.from_dict(
        annotation_data,
        orient="columns"
    )
    # Sort index to match original objects, then reset index to make IDs a regular column
    df.index.name = "id"
    df = df.loc[objects.list_attr("id")].reset_index()
    # Replace empty strings with NaN
    df = df.replace("", float("nan"))
    return df
