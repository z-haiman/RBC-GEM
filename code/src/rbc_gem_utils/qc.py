from collections import Counter

import pandas as pd
from cobra.core.group import Group
from cobra.core.metabolite import element_re

from .util import COBRA_CONFIGURATION


def standardardize_metabolite_formulas(metabolite_formulas):
    """Standardize metabolite chemical formulas by reordering elements according to the Hill system.

    Parameters
    ----------
    metabolite_formulas : dict
        Contains identifiers as keys and strings representing chemical formulas as values.

    Returns
    -------
    updated_formulas : dict
        Contains original keys and standardized chemical formulas as values.

    Raises
    ------
    ValueError
        Occurs when metabolite(s) have missing or non-alphanumeric formulas.

    """
    invalid_formula = {
        key: formula
        for key, formula in metabolite_formulas.items()
        if not formula  # validate formula is neither `None` nor an empty string
        or not formula.isalnum()  # validate formula is alphanumeric, enforces integer coefficients
    }
    if invalid_formula:
        raise ValueError(
            f"Missing or invalid formulas for the following: {list(invalid_formula)}"
        )

    updated_formulas = {}
    for key, formula in metabolite_formulas.items():
        # Parse formula
        parsed = element_re.findall(formula)
        # Identify whether carbons are in the formula
        element_counter = Counter()
        for e, n in parsed:
            if e not in element_counter:
                element_counter[e] = 0
            element_counter[e] += int(n) if n != "" else 1

        if element_counter.get("C"):
            # Presence of carbons, organize first carbons, then hydrogens, and alphabetically sort the remaining.
            carbons = [("C", element_counter.pop("C"))]
            hydrogens = (
                [("H", element_counter.pop("H"))] if element_counter.get("H") else []
            )
            parsed = carbons + hydrogens + sorted(element_counter.items())
        else:
            # Alphabetically sort, regardless of the presence of hydrogens.
            parsed = sorted(element_counter.items())

        # Add formula to return.
        updated_formulas[key] = "".join(e if n == 1 else f"{e}{n}" for e, n in parsed)

    # Return updated formulas as a dictionary with original keys
    return updated_formulas


def compare_series(old, new, to_compare=None):
    """Compare two series."""
    name = set([old.name, new.name])
    if len(name) == 2:
        raise ValueError(f"Series have different names: {name}")
    else:
        name = name.pop()

    # Determine all possible indicies, validate indicies to compare
    if to_compare is None:
        to_compare = sorted(old.index.union(new.index))

    merged_to_compare = pd.merge(
        old,
        new,
        left_index=True,
        right_index=True,
        how="outer",
        suffixes=("_OLD", "_NEW"),
    )
    # TODO cleanup code, probably a better way to do this
    # TODO update to have one location for values and colors
    compared = pd.DataFrame([], index=merged_to_compare.index, columns=[name])
    # First, record all changes that have/have not been made.
    compared.loc[
        merged_to_compare[f"{name}_NEW"].fillna("").astype(str)
        == merged_to_compare[f"{name}_OLD"].fillna("").astype(str),
        name,
    ] = "NO CHANGE"
    compared.loc[
        merged_to_compare[f"{name}_NEW"].fillna("").astype(str)
        != merged_to_compare[f"{name}_OLD"].fillna("").astype(str),
        name,
    ] = "CHANGED"
    # Determine removed entries, present with a value in the old but gone in the new
    compared.loc[old.index.difference(new.index), name] = "REMOVED"
    # Determine added entries, present with a value in the new but not present in the old
    compared.loc[new.index.difference(old.index), name] = "ADDED"
    # Determine new entries that were made as blanks
    compared.loc[new[new.isna()].index, name] = "EMPTY"
    # Entries that were intentionally changed from having a value to being blank should be marked as removed
    compared.loc[
        old[~old.isna()].index.intersection(new[new.isna()].index), name
    ] = "REMOVED"
    # All remaining entries are EMPTY in both datasets
    missing = set(to_compare).difference(compared.index)
    if missing:
        compared = pd.concat(
            (compared, pd.DataFrame([], index=list(missing), columns=[name]))
        )
    compared = compared.fillna("EMPTY").sort_index(ascending=True)
    return compared


def compare_tables(old, new, to_compare=None):
    """Compare two tables for differences."""
    # Determine all possible columns, validate columns to compare
    all_indicies = sorted(set(old.index).union(new.index))
    all_columns = sorted(set(old.columns).union(new.columns))
    if to_compare is not None:
        to_compare = set(to_compare)
        if to_compare.difference(all_columns):
            raise ValueError(
                f"columns {to_compare.difference(all_columns)} not found in either table."
            )
    else:
        to_compare = all_columns

    compared = pd.DataFrame([], columns=sorted(to_compare))
    for label in compared.columns:
        compared[label] = compare_series(
            old.get(label, pd.Series([], name=label)),
            new.get(label, pd.Series([], name=label)),
            to_compare=all_indicies,
        )

    return compared


def reset_subsystem_groups(model):
    model.remove_groups(model.groups)
    for subsystem in sorted(set(model.reactions.list_attr("subsystem"))):
        reaction_list = model.reactions.query(lambda x: x.subsystem == subsystem)
        if subsystem not in model.groups:
            group = Group(id=subsystem, name=subsystem, members=reaction_list)
            model.add_groups([group])
        else:
            group = model.groups.get_by_id(subsystem).add_members(reaction_list)


def reset_reaction_bounds(model):
    for reaction in model.reactions:
        if reaction.bounds == COBRA_CONFIGURATION.bounds or reaction.bounds == (
            0.0,
            COBRA_CONFIGURATION.upper_bound,
        ):
            # Already at default, no need to change
            continue
        print(f"Before: {reaction}")
        if reaction.boundary:
            if reaction.id.startswith("DM_"):
                reaction.bounds = (0.0, COBRA_CONFIGURATION.upper_bound)
            elif reaction.id.startswith("EX_") or reaction.id.startswith("SK_"):
                reaction.bounds = COBRA_CONFIGURATION.bounds
            else:
                print(f"Unreccognized boundary type for {reaction}")
        elif not reaction.reversibility:
            if "<--" in reaction.reaction:
                reaction.bounds = (COBRA_CONFIGURATION.lower_bound, 0.0)
            else:
                reaction.bounds = (0.0, COBRA_CONFIGURATION.upper_bound)
        else:
            reaction.bounds = (
                COBRA_CONFIGURATION.lower_bound,
                COBRA_CONFIGURATION.upper_bound,
            )
