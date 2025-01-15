import itertools

import pandas as pd
from sympy import logic

from .util import COBRA_CONFIGURATION, build_string, split_string


def create_table_of_reaction_bounds(
    model, default_bounds=False, additional_cols=None, sort=True
):
    """Create a table for the reaction bounds of the model.

    Parameters
    ----------
    model : cobra.Model

    default_bounds : bool
        Whether to include reactions with default bounds.
        Default is False.

    additional_cols : dict, optional
        Additional columns to append to the table. Keys correspond to column names.
        Values must also be dictionaries, with reaction IDs mapped to entry values.

    Returns
    -------
    df
        table
    """
    reactions = model.reactions
    if not default_bounds:
        default_bounds = set(
            [
                COBRA_CONFIGURATION.bounds,  # Reversible
                (0.0, COBRA_CONFIGURATION.upper_bound),  # Irreversible forward
                (COBRA_CONFIGURATION.lower_bound, 0.0),  # Irreversible reverse
            ]
        )
        reactions = reactions.query(lambda x: x.bounds not in default_bounds)

    data = {}
    for idx, reaction in enumerate(reactions):
        if get_btype_by_prefix(reaction):
            reaction_type = get_btype_by_prefix(reaction)
        elif get_btype_by_prefix(reaction):
            reaction_type = get_btype_by_sbo(reaction)
        elif reaction.id.startswith("LOAD_"):
            reaction_type = "load"
        elif len(reaction.compartments) >= 2:
            reaction_type = "transport"
        else:
            reaction_type = "biochemical"

        data[idx] = {
            "reactions": reaction.id,
            "name": reaction.name,
            "type": reaction_type,
            "reaction": reaction.reaction,
            "subsystem": reaction.subsystem,
            "lower_bound": reaction.lower_bound,
            "upper_bound": reaction.upper_bound,
        }

    df = pd.DataFrame.from_dict(data, orient="index")
    # Set additional columns, these should be included in all tables.
    if additional_cols:
        for key, value in additional_cols.items():
            df = pd.merge(
                df,
                pd.DataFrame.from_dict(value, columns=[key], orient="index"),
                left_on="reactions",
                right_index=True,
                how="left",
            )
    if sort:
        df = sort_df_boundary_by_type_and_id(df)
    else:
        df = df.reset_index(drop=True)
    return df


def get_btype_by_prefix(boundary):
    if boundary.id.startswith("EX_"):
        return "exchange"
    elif boundary.id.startswith("DM_"):
        return "demand"
    elif boundary.id.startswith("SK_"):
        return "sink"
    else:
        return ""


def get_btype_by_sbo(boundary):
    sbo = boundary.annotation.get("sbo")
    if sbo == "SBO:0000627":
        return "exchange"
    elif sbo == "SBO:0000628":
        return "demand"
    elif sbo == "SBO:0000632":
        return "sink"
    else:
        return None


def sort_df_boundary_by_type_and_id(df_boundary):
    df_boundary = df_boundary.sort_values(by="reactions")
    df_boundary = pd.concat(
        (
            df_boundary[df_boundary["type"] == "exchange"],
            df_boundary[df_boundary["type"] == "sink"],
            df_boundary[df_boundary["type"] == "demand"],
            df_boundary[~df_boundary["type"].isin(["exchange", "sink", "demand"])],
        )
    ).reset_index(drop=True)
    return df_boundary


def parse_gprs_to_complexes(
    model,
    genes_to_proteins=None,
    cofactor_genes=None,
    additional_cols=None,
    replace_compartments=None,
):
    """Create a table for the complexes of the model.

    Parameters
    ----------
    model : cobra.Model
        The model to use. Complexes are extracted form the gene reaction rules.
    genes_to_proteins : dict
        Maps gene identifiers to protein identifiers.
    cofactor_genes : set, dict, optional
        Gene-encoded proteins that are necessary for a protein or complex to catalyze a reaction.
        Specifically, these proteins act more in a substrate-like role and may not be formally recognized as part of the complex itself.
        Examples can include Thioredoxin (TXN), Glutaredoxin (GLRX), Polyubiquitin-B (UBB), etc.
    additional_cols : dict, optional
        Additional columns to append to the table. Keys correspond to column names.
        Values must also be dictionaries, with reaction IDs mapped to entry values.
    replace_compartments : dict, optional
        Used to replaced compartment identifiers.
        Ideal for compartments merged together when a complex spans a membrane.

    Returns
    -------
    df
        table
    """
    if genes_to_proteins is None:
        genes_to_proteins = {}
    else:
        genes_to_proteins = {
            k: tuple(split_string(v)) for k, v in genes_to_proteins.items()
        }

    data = {}
    for idx, reaction in enumerate(
        model.reactions.query(lambda x: x.gene_reaction_rule)
    ):
        # Use sympy to parse GPR to complexes
        symbolic_gpr = reaction.gpr.as_symbolic()
        # Use DNF form
        symbolic_gpr = logic.to_dnf(symbolic_gpr)
        cplx_list = [
            x.strip("(").strip(")").strip() for x in str(symbolic_gpr).split(" | ")
        ]
        data[idx] = {
            "genes": cplx_list,
            "reactions": reaction.id,
            "compartment": "".join((sorted(reaction.compartments))),
        }

    df = pd.DataFrame.from_dict(data, orient="index")
    df = df.explode("genes")
    df["genes"] = df["genes"].apply(lambda x: set(x.split(" & ")))
    df["subunits"] = df["genes"].apply(
        lambda x: [genes_to_proteins.get(gid, gid) for gid in sorted(x)]
    )
    df["subunits"] = df["subunits"].apply(
        lambda x: list(set(combo) for combo in itertools.product(*x))
    )
    df = df.explode("subunits")
    if cofactor_genes is None:
        df["subunits"] = df["subunits"].apply(build_string)
    else:
        if isinstance(cofactor_genes, str):
            cofactor_genes = split_string(cofactor_genes)
        cofactor_genes = set(
            [genes_to_proteins.get(gid, gid) for gid in cofactor_genes]
        )
        df["cofactors"] = df["subunits"].apply(
            lambda x: build_string(sorted(x.intersection(cofactor_genes)))
        )
        df["subunits"] = df["subunits"].apply(
            lambda x: build_string(sorted(x.difference(cofactor_genes)))
        )
    df["genes"] = df["genes"].apply(lambda x: build_string(sorted(x)))

    # Set additional columns, these should be included in all tables.
    if additional_cols:
        for key, value in additional_cols.items():
            df[key] = value

    if replace_compartments:
        df["compartment"] = df["compartment"].replace(replace_compartments)

    df["compartment"] = df["compartment"].apply(
        lambda x: build_string(sorted(split_string(x)))
    )
    df = df[
        df.apply(
            lambda row: all(
                [x.endswith(row["compartment"]) for x in split_string(row["subunits"])]
            ),
            axis=1,
        )
    ]
    df = (
        df.replace("", pd.NA)
        .groupby(["genes", "subunits", "compartment"], as_index=False)
        .agg(lambda x: build_string(x.dropna()))
    )
    df = df.sort_values("genes").reset_index(drop=True)
    return df
