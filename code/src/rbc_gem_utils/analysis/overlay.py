"""TODO DOCSTRING."""

import logging
from collections import defaultdict
from warnings import warn

import numpy as np
import pandas as pd
from cobra import Metabolite, Reaction, manipulation

from rbc_gem_utils.io import read_cobra_model
from rbc_gem_utils.table import parse_gprs_to_complexes
from rbc_gem_utils.util import (
    build_string,
    ensure_iterable,
    explode_column,
    log_msg,
    split_string,
    strip_plural,
)


LOGGER = logging.getLogger(__name__)


class ComplexForm(Reaction):
    """Special reaction class representing the formation/concentration of complexes."""

    def __init__(
        self,
        id: str | None = None,
        name: str = "",
        subsystem: str = "",
        lower_bound: float = 0,
        upper_bound: float | None = None,
        **kwargs,
    ) -> None:
        super().__init__(id, name, subsystem, lower_bound, upper_bound, **kwargs)


class EnzymeForm(Reaction):
    """Special reaction class representing the formation/concentration of enzymes."""

    def __init__(
        self,
        id: str | None = None,
        name: str = "",
        subsystem: str = "",
        lower_bound: float = 0,
        upper_bound: float | None = None,
        **kwargs,
    ) -> None:
        super().__init__(id, name, subsystem, lower_bound, upper_bound, **kwargs)


class ProteinDilution(Reaction):
    """Special reaction class representing the formation/concentration of proteins."""

    def __init__(
        self,
        id: str | None = None,
        name: str = "",
        subsystem: str = "",
        lower_bound: float = 0,
        upper_bound: float | None = None,
        **kwargs,
    ) -> None:
        super().__init__(id, name, subsystem, lower_bound, upper_bound, **kwargs)


class ComplexDilution(Reaction):
    """Special reaction class representing the formation/concentration of complexes."""

    def __init__(
        self,
        id: str | None = None,
        name: str = "",
        subsystem: str = "",
        lower_bound: float = 0,
        upper_bound: float | None = None,
        **kwargs,
    ) -> None:
        super().__init__(id, name, subsystem, lower_bound, upper_bound, **kwargs)


class EnzymeDilution(Reaction):
    """Special reaction class representing the concentration of enzymes."""

    def __init__(
        self,
        id: str | None = None,
        name: str = "",
        subsystem: str = "",
        lower_bound: float = 0,
        upper_bound: float | None = None,
        **kwargs,
    ) -> None:
        super().__init__(id, name, subsystem, lower_bound, upper_bound, **kwargs)


class ProteomeBudgetDilution(Reaction):
    """Special reaction class representing the concentration of enzymes."""

    def __init__(
        self,
        id: str | None = None,
        name: str = "",
        subsystem: str = "",
        lower_bound: float = 0,
        upper_bound: float | None = None,
        **kwargs,
    ) -> None:
        super().__init__(id, name, subsystem, lower_bound, upper_bound, **kwargs)


class Protein(Metabolite):
    """TODO DOCSTRING."""

    def __init__(
        self,
        id: str | None = None,
        formula: str | None = None,
        name: str | None = "",
        charge: float | None = None,
        compartment: str | None = None,
    ) -> None:
        super().__init__(id, formula, name, charge, compartment)


class Complex(Metabolite):
    """TODO DOCSTRING."""

    def __init__(
        self,
        id: str | None = None,
        formula: str | None = None,
        name: str | None = "",
        charge: float | None = None,
        compartment: str | None = None,
    ) -> None:
        super().__init__(id, formula, name, charge, compartment)


class Enzyme(Metabolite):
    """TODO DOCSTRING."""

    def __init__(
        self,
        id: str | None = None,
        formula: str | None = None,
        name: str | None = "",
        charge: float | None = None,
        compartment: str | None = None,
    ) -> None:
        super().__init__(id, formula, name, charge, compartment)


class ProteomeBudget(Metabolite):
    """TODO DOCSTRING."""

    def __init__(
        self,
        id: str | None = None,
        formula: str | None = None,
        name: str | None = "",
        charge: float | None = None,
        compartment: str | None = None,
    ) -> None:
        super().__init__(id, formula, name, charge, compartment)


ATTRIBUTE_MAPPING = {
    "metabolites": ["id", "name", "formula", "charge", "compartment"],
    "reactions": ["id", "name", "subsystem", "lower_bound", "upper_bound"],
}

TABLE_COLUMNS = {
    "proteins": {
        "required": ["genes", "protein", "compartment", "molar_mass"],
        "optional": ["proteins", "reactions", "sequence.id", "sequence"],
    },
    "complexes": {
        "required": ["complex", "compartment", "stoichiometry", "reactions", "genes"],
        "optional": ["complexes", "subunits", "coefficients", "cofactors"],
    },
    "enzymes": {
        "required": ["enzyme", "compartment", "reactions"],  # "complexes_keff",
        "optional": ["enzymes", "complexes", "direction"],  # "complex_keff",
    },
}

DEFAULT_MAX_WEIGHT_FRACTION = 150  # 150 mg / gDW
DEFAULT_CONCENTRATION_BOUND = 1e6  # either nmol / gDW or nmol / L
DEFAULT_KEFF = 65 * 3600  # 1 / s  # s / hr
DEFAULT_PROTEIN_MOLAR_MASS = 30000  # g / mol --> mg / mmol
DEFAULT_PREFIX_SUFFIX_VALUES = {
    "proteins": {
        "prefix.metabolite": "protein_",
        "prefix.dilution": "PROTDL_",
        "prefix.formation": "PROTFM_",
    },
    "complexes": {
        "prefix.metabolite": "cplx_",
        "prefix.dilution": "CPLXDL_",
        "prefix.formation": "CPLXFM_",
    },
    "enzymes": {
        "prefix.metabolite": "enzyme_",
        "prefix.dilution": "ENZDL_",
        "prefix.formation": "ENZFM_",
        "suffix.forward": "_fwd",
        "suffix.reverse": "_rev",
        "suffix.total": "_total",
    },
    "constraints": {
        "prefix.constraint": "CONS_",
        "prefix.isoform": "ISOCONS_",
        "prefix.compartent": f"COMPCONS_",
    },
}

ATTR_SUBCLASS_DICT = {
    "metabolites": {
        "proteins": Protein,
        "complexes": Complex,
        "enzymes": Enzyme,
        "proteome budget": ProteomeBudget,
    },
    "reactions": {
        "protein dilutions": ProteinDilution,
        "complex formation reactions": ComplexForm,
        "complex dilutions": ComplexDilution,
        "enzyme formation reactions": EnzymeForm,
        "enzyme dilutions": EnzymeDilution,
        "proteome budget demand": ProteomeBudgetDilution,
    },
}

AA_MOLECULAR_FORMULAS = {
    "A": "C3H7NO2",  # ala__L_c
    "R": "C6H15N4O2",  # arg__L_c
    "N": "C4H8N2O3",  # asn__L_c
    "D": "C4H6NO4",  # asp__L_c
    "C": "C3H7NO2S",  # cys__L_c
    "E": "C5H8NO4",  # glu__L_c
    "Q": "C5H10N2O3",  # gln__L_c
    "G": "C2H5NO2",  # gly_c
    "H": "C6H9N3O2",  # his__L_c
    "I": "C6H13NO2",  # ile__L_c
    "L": "C6H13NO2",  # leu__L_c
    "K": "C6H15N2O2",  # lys__L_c
    "M": "C5H11NO2S",  # met__L_c
    "F": "C9H11NO2",  # phe__L_c
    "P": "C5H9NO2",  # pro__L_c
    "S": "C3H7NO3",  # ser__L_c
    "T": "C4H9NO3",  # thr__L_c
    "W": "C11H12N2O2",  # trp__L_c
    "Y": "C9H11NO3",  # tyr__L_c
    "V": "C5H11NO2",  # val__L_c
}
# Selenocysteine
# AA_MOLECULAR_CHARGES["X"] =    # for unknown AAs
AA_MOLECULAR_FORMULAS["B"] = (
    "C4H6NO4",
)  # Asparagine/aspartate (Asx) indistinguishable, deaminated residue
AA_MOLECULAR_FORMULAS["Z"] = (
    "C5H8NO4"  # Glutamine/glutamate (Glx) indistinguishable, deaminated residue
)
AA_MOLECULAR_FORMULAS["U"] = "C3H7NO2Se"  # Selenocysteine


AA_MOLECULAR_CHARGES = {
    "A": 0,  # ala__L_c
    "R": 1,  # arg__L_c
    "N": 0,  # asn__L_c
    "D": -1,  # asp__L_c
    "C": 0,  # cys__L_c
    "E": -1,  # glu__L_c
    "Q": 0,  # gln__L_c
    "G": 0,  # gly_c
    "H": 0,  # his__L_c
    "I": 0,  # ile__L_c
    "L": 0,  # leu__L_c
    "K": 1,  # lys__L_c
    "M": 0,  # met__L_c
    "F": 0,  # phe__L_c
    "P": 0,  # pro__L_c
    "S": 0,  # ser__L_c
    "T": 0,  # thr__L_c
    "W": 0,  # trp__L_c
    "Y": 0,  # tyr__L_c
    "V": 0,  # val__L_c
}
AA_MOLECULAR_CHARGES["X"] = 0  # Neutral for unknown AAs
AA_MOLECULAR_CHARGES["B"] = (
    -1,
)  # Asparagine/aspartate (Asx) indistinguishable, deaminated residue
AA_MOLECULAR_CHARGES["Z"] = (
    -1
)  # Glutamine/glutamate (Glx) indistinguishable, deaminated residue
AA_MOLECULAR_CHARGES["U"] = 0  # Selenocysteine


AA_MOLECULAR_WEIGHTS = {
    # Weight of amino acid minus one H2O
    "A": 89.093 - 18.015,  # ala__L_c
    "R": 175.209 - 18.015,  # arg__L_c
    "N": 132.118 - 18.015,  # asn__L_c
    "D": 132.095 - 18.015,  # asp__L_c
    "C": 121.158 - 18.015,  # cys__L_c
    "E": 146.121 - 18.015,  # glu__L_c
    "Q": 146.144 - 18.015,  # gln__L_c
    "G": 75.067 - 18.015,  # gly_c
    "H": 155.155 - 18.015,  # his__L_c
    "I": 131.173 - 18.015,  # ile__L_c
    "L": 131.173 - 18.015,  # leu__L_c
    "K": 147.196 - 18.015,  # lys__L_c
    "M": 149.211 - 18.015,  # met__L_c
    "F": 165.189 - 18.015,  # phe__L_c
    "P": 115.13 - 18.015,  # pro__L_c
    "S": 105.093 - 18.015,  # ser__L_c
    "T": 119.119 - 18.015,  # thr__L_c
    "W": 204.225 - 18.015,  # trp__L_c
    "Y": 181.189 - 18.015,  # tyr__L_c
    "V": 117.146 - 18.015,  # val__L_c
}
# Average for unknown AAs
AA_MOLECULAR_WEIGHTS["X"] = sum(AA_MOLECULAR_WEIGHTS.values()) / len(
    AA_MOLECULAR_WEIGHTS
)
# Asparagine/aspartate (Asx) indistinguishable
AA_MOLECULAR_WEIGHTS["B"] = (
    sum([AA_MOLECULAR_WEIGHTS["N"], AA_MOLECULAR_WEIGHTS["D"]]) / 2
)
AA_MOLECULAR_WEIGHTS["Z"] = (
    sum([AA_MOLECULAR_WEIGHTS["Q"], AA_MOLECULAR_WEIGHTS["E"]]) / 2
)
# Selenocysteine
AA_MOLECULAR_WEIGHTS["U"] = 121.16 - 18.015
AA_MOLECULAR_WEIGHTS = {k: round(v, 3) for k, v in AA_MOLECULAR_WEIGHTS.items()}
AA_MOLECULAR_WEIGHTS


def calculate_protein_molar_mass(aa_seq):
    """Compute protein mass in Daltons (g / mol) from an amino acid sequence.

    Adds the mass of one water molecule to the sequence.

    Parameters
    ----------
    aa_seq : str
        The amino acid sequence as a string. The sequence must be in one-letter symbols.

    Returns
    -------
    float :
        The molecular weight of the protein. The return unit is g/mol.

    Notes
    -----
    Based on PROSO toolbox `calcProteinMM.m`.
    https://github.com/QCSB/PROSO-Toolbox/blob/main/src/PC-FBA/parseGeneRule.m
    """

    return (
        sum(
            [
                # Default to using the average if unknown/unreadable amino acid in sequence.
                AA_MOLECULAR_WEIGHTS.get(aa, AA_MOLECULAR_WEIGHTS["X"])
                for aa in aa_seq
            ]
        )
        + 18.015
    )  # One water not removed


def create_variable_value_mapping(*args, sep=";"):
    """Parse a formatted string or set of strings to map variables identifiers to values"""
    if isinstance(args[0], list):
        args = args[0]
    if isinstance(args, str):
        args = [x for x in split_string(args)]

    if len(args) == 1:
        args = np.array(
            [
                [x.split("(")[0], x.split("(")[1].strip(")")]
                for x in split_string(args[0])
            ]
        ).T

    if len(args) == 2:
        to_zip = []
        for arg in args:
            if isinstance(arg, str):
                arg = arg.split(sep)
            if all([x.isnumeric() for x in arg]):
                arg = [int(x) if float(x) == int(x) else float(x) for x in arg]
            to_zip += [arg]
    return dict(zip(*to_zip))


def create_protein_table(
    model,
    df_protein_data,
    id_key=None,
    prefix=None,
    optional_columns=False,
    annotation_columns=None,
    replace_compartments=None,
):
    """Create a table of genes mapped to proteins.

    Parameters
    ----------
    model : cobra.Model
        The model to use.
    df_protein_data : pandas.DataFrame, None
        DataFrame of sequence data already formatted.
    id_key : str, None
        The column to use when generating the IDs.
        Must be {"genes", "sequence.id"} or a value provided in `annotation_columns`.
        If ``None`` provided, generic identifiers are used instead and prefix will be considered ``True``,
    prefix : str, optional
        Prefix to use for identifiers. Default is defined in :const:`DEFAULT_PREFIX_SUFFIX_VALUES`.
    optional_columns : list, bool
        Optional columns to use when generating table.
        Possible values are defined in :const:`TABLE_COLUMNS`.
        Passing ``True`` will ensure all optional columns are included, while  ``False`` ensures none are included.
    annotation_columns : list, bool, optional
        Additional columns to include from gene annotations when generating table (e.g. ``"uniprot"``).
        Must be one or keys in the annotation dictionaries.
    replace_compartments : dict, optional
        Used to replaced compartment identifiers.
        Ideal for compartments merged together when a complex spans a membrane.

    Returns
    -------
    table : pandas.DataFrame
        Table of proteins with desired columns returned as a DataFrame.

    Notes
    -----
    Based on PROSO toolbox `pcModel.m`.
    https://github.com/QCSB/PROSO-Toolbox/blob/main/src/PC-FBA/pcModel.m
    """
    # Any additional annotations or attributes that should be included. Can be used to add annotations and other attribute columns to objects.
    table_type = "proteins"
    annotation_columns = (
        [] if not annotation_columns else ensure_iterable(annotation_columns)
    )
    if optional_columns is True:
        optional_columns = TABLE_COLUMNS[table_type]["optional"]
    elif optional_columns and set(optional_columns).difference(
        TABLE_COLUMNS[table_type]["optional"]
    ):
        raise ValueError(
            f"`optional_columns` must be one or more of the following: {TABLE_COLUMNS[table_type]['optional']}"
        )
    elif not optional_columns:
        optional_columns = []
    else:
        optional_columns = ensure_iterable(optional_columns)

    valid = {"genes"}
    if df_protein_data is not None:
        valid = valid.union(["sequence.id", "protein.id"])
    valid = valid.union(annotation_columns)
    if id_key is not None and id_key not in valid:
        raise ValueError(f"`id_key` must be one of the following: {valid}")

    # Extract data
    data = {}
    for idx, gene in enumerate(model.genes):
        data[idx] = {"genes": gene.id}
        data[idx].update(
            {
                attr: getattr(gene, attr)
                for attr in set(optional_columns).union(["reactions"])
                if hasattr(gene, attr)
            }
        )
        data[idx].update({"proteins": gene.id})
        data[idx].update(
            {
                key: build_string(gene.annotation.get(key, []))
                for key in annotation_columns
            }
        )
    table = pd.DataFrame.from_dict(data, orient="index")
    # Explode compartments, then seperate reactions based on compartments
    table["compartment"] = table["reactions"].apply(
        lambda reaction_list: build_string(
            sorted(
                [
                    "".join(sorted((getattr(rxn, "compartments", rxn))))
                    for rxn in reaction_list
                ]
            )
        )
    )
    table["compartment"] = table["compartment"].apply(split_string)
    table = table.explode("compartment")
    table["reactions"] = table[["reactions", "compartment"]].apply(
        lambda x: build_string(
            sorted(
                [
                    getattr(rxn, "id", rxn)
                    for rxn in x["reactions"]
                    if "".join(sorted(rxn.compartments)) == x["compartment"]
                ]
            )
        ),
        axis=1,
    )
    if replace_compartments:
        table["compartment"] = table["compartment"].replace(replace_compartments)
        table = table.groupby(["genes", "compartment"], as_index=False).agg(
            lambda x: build_string(split_string(x))
        )

    # Identifiers
    if df_protein_data is not None:
        table = table.merge(df_protein_data, how="left")

    table_type_strip = strip_plural(table_type)
    if not id_key:
        prefix = (
            prefix
            if (prefix and isinstance(prefix, str))
            else DEFAULT_PREFIX_SUFFIX_VALUES[table_type]["prefix.metabolite"]
        )
        table[table_type_strip] = [
            f"{prefix}{x}".replace("-", "_") for x in table.index
        ]
    else:
        if prefix:
            prefix = (
                prefix
                if (prefix and isinstance(prefix, str))
                else DEFAULT_PREFIX_SUFFIX_VALUES[table_type]["prefix.metabolite"]
            )
        else:
            prefix = ""
        table[table_type_strip] = table[id_key].apply(
            lambda x: f"{prefix}{x}".replace("-", "_")
        )

    # Clean up for consistency
    required_columns = TABLE_COLUMNS[table_type]["required"]
    columns = required_columns + optional_columns + annotation_columns
    for col in columns:
        if table.get(col) is None:
            table[col] = pd.NA

    # No attributes to annotate
    table["compartment"] = table["compartment"].apply(
        lambda x: "".join(sorted(set(split_string(x))))
    )

    table[table_type] = table[[table_type_strip, "compartment"]].apply(
        lambda x: "_".join(x), axis=1
    )
    table = table.replace(float("nan"), pd.NA).replace("", pd.NA)
    table = (
        table.loc[:, columns]
        .drop_duplicates()
        .sort_values(required_columns)
        .reset_index(drop=True)
    )
    return table


# Gene protein mapping is necessary for correct subunit mapping, can be obtained from protein table
def create_complex_table(
    model,
    genes_to_proteins,
    cofactor_genes=None,
    id_key=None,
    prefix=None,
    optional_columns=False,
    annotation_columns=None,
    replace_compartments=None,
):
    """Create a table of complexes mapped to proteins.

    Parameters
    ----------
    model : cobra.Model
        The model to use.
    genes_to_proteins : dict
        Maps gene identifiers to protein identifiers. For multiple identifiers, values should be seperated by a semicolin
    cofactor_genes : set, dict, optional
        Gene-encoded proteins that are necessary for the protein complex to catalyze a reaction.
        Specifically, these proteins act more in a substrate-like role and may not be formally recognized as part of the complex itself.
        Examples can include Thioredoxin (TXN), Glutaredoxin (GLRX), Polyubiquitin-B (UBB), etc.
    id_key : str, None
        The column to use when generating the IDs.
        Must be a value provided in `annotation_columns`.
        If ``None`` provided, generic identifiers are used instead and prefix will be considered ``True``,
    prefix : str, optional
        Prefix to use for identifiers. Default is defined in :const:`DEFAULT_PREFIX_SUFFIX_VALUES`.
    optional_columns : list, bool
        Optional columns to use when generating table.
        Possible values are defined in :const:`TABLE_COLUMNS`.
        Passing ``True`` will ensure all optional columns are included, while  ``False`` ensures none are included.
    annotation_columns : list, bool, optional
        Additional columns to include from gene annotations when generating table (e.g. ``"uniprot"``).
        Must be one or keys in the annotation dictionaries.
    replace_compartments : dict, optional
        Used to replaced compartment identifiers.
        Ideal for compartments merged together when a complex spans a membrane.
        If used in combination with `genes_to_proteins`, ensure compartments match.

    Returns
    -------
    table : pandas.DataFrame
        Table of complexes with desired columns returned as a DataFrame.

    Notes
    -----
    Based on PROSO toolbox `pcModel.m`.
    https://github.com/QCSB/PROSO-Toolbox/blob/main/src/PC-FBA/pcModel.m
    """
    # Any additional annotations or attributes that should be included. Can be used to add annotations and other attribute columns to objects.
    table_type = "complexes"
    annotation_columns = (
        [] if not annotation_columns else ensure_iterable(annotation_columns)
    )
    if optional_columns is True:
        optional_columns = TABLE_COLUMNS[table_type]["optional"]
    elif optional_columns and set(optional_columns).difference(
        TABLE_COLUMNS[table_type]["optional"]
    ):
        raise ValueError(
            f"`optional_columns` must be one or more of the following: {TABLE_COLUMNS[table_type]['optional']}"
        )
    elif not optional_columns:
        optional_columns = []
    else:
        optional_columns = ensure_iterable(optional_columns)

    valid = set(annotation_columns)
    if id_key is not None and id_key not in valid:
        raise ValueError(f"`id_key` must be one of the following: {valid}")
    valid = {"genes"}.union(annotation_columns + [id_key] if id_key is not None else [])

    table = parse_gprs_to_complexes(
        model,
        genes_to_proteins=genes_to_proteins,
        cofactor_genes=cofactor_genes,
        replace_compartments=replace_compartments,
    )

    table["coefficients"] = table["subunits"].apply(
        lambda x: ";".join(["1" for _ in range(len(split_string(x)))])
    )
    table["stoichiometry"] = table.apply(
        lambda x: build_string(
            [
                f"{k}({v})"
                for k, v in create_variable_value_mapping(
                    *x[["subunits", "coefficients"]].values
                ).items()
            ]
        ),
        axis=1,
    )
    table["cofactors"] = table["cofactors"].apply(
        lambda x: build_string(
            [genes_to_proteins.get(gid, gid) for gid in split_string(x) if gid]
        )
    )
    # Identifiers
    table_type_strip = strip_plural(table_type)
    if not id_key:
        prefix = (
            prefix
            if (prefix and isinstance(prefix, str))
            else DEFAULT_PREFIX_SUFFIX_VALUES[table_type]["prefix.metabolite"]
        )
        table[table_type_strip] = [
            f"{prefix}{x}".replace("-", "_") for x in table.index
        ]
    else:
        if prefix:
            prefix = (
                prefix
                if (prefix and isinstance(prefix, str))
                else DEFAULT_PREFIX_SUFFIX_VALUES[table_type]["prefix.metabolite"]
            )
        else:
            prefix = ""
        table[table_type_strip] = table[id_key].apply(
            lambda x: f"{prefix}{x}".replace("-", "_")
        )

    # Clean up for consistency and add missing annotations
    required_columns = TABLE_COLUMNS[table_type]["required"]
    columns = required_columns + optional_columns + annotation_columns
    for col in columns:
        if table.get(col) is None:
            table[col] = pd.NA

    attr = "genes"
    for col in annotation_columns:
        table[col] = table[attr].apply(
            lambda x: [
                getattr(model, attr).get_by_id(sid).annotation.get(col)
                for sid in split_string(x)
                if getattr(model, attr).has_id(sid)
                and getattr(model, attr).get_by_id(sid).annotation.get(col)
            ]
        )
        table[col] = table[col].apply(
            lambda x: build_string(
                sorted(
                    set([item for sublist in x for item in ensure_iterable(sublist)])
                )
            )
        )

    if replace_compartments:
        table["compartment"] = table["compartment"].replace(replace_compartments)
    table["compartment"] = table["compartment"].apply(
        lambda x: "".join(sorted(set(split_string(x))))
    )
    table[table_type] = table[[table_type_strip, "compartment"]].apply(
        lambda x: "_".join(x), axis=1
    )
    table = table.replace(float("nan"), pd.NA).replace("", pd.NA)
    table = (
        table.loc[:, columns]
        .drop_duplicates()
        .sort_values(required_columns)
        .reset_index(drop=True)
    )
    return table


def create_enzyme_table(
    model,
    complexes_to_reactions,
    enzyme_keff_base=None,
    enzyme_forward_suffix=None,
    enzyme_reverse_suffix=None,
    id_key=None,
    prefix=None,
    optional_columns=False,
    annotation_columns=None,
    replace_compartments=None,
):
    """Create a table of enzymes mapped to complexes.

    Parameters
    ----------
    model : cobra.Model
        The model to use.
    complexes_to_reactions : dict
        Maps complex identifiers to a list of associated reactions.
    enzyme_keff_base : float
        Baseline value for the effective rate constant of the representative "enzyme".
        Default is `234000` 1/hr (equivalently, 65 1/s).
    enzyme_reverse_suffix : str, optional
        Suffix to use for identifiers of enzymes in the reverse direction. Default is ``"_rev"``.
    id_key : str, None
        The column to use when generating the protein IDs.
        Must be {"reactions"} or a value provided in `annotation_columns`.
        If ``None`` provided, generic identifiers are used instead and prefix will be considered ``True``,
    prefix : str, optional
        Prefix to use for identifiers. Default is defined in :const:`DEFAULT_PREFIX_SUFFIX_VALUES`.
    optional_columns : list, bool
        Optional columns to use when generating table.
        Possible values are defined in :const:`TABLE_COLUMNS`.
        Passing ``True`` will ensure all optional columns are included, while  ``False`` ensures none are included.
    annotation_columns : list, bool, optional
        Additional columns to include from gene annotations when generating table (e.g. ``"uniprot"``).
        Must be one or keys in the annotation dictionaries.
    replace_compartments : dict, optional
        Used to replaced compartment identifiers.
        Ideal for compartments merged together when a complex spans a membrane.

    Returns
    -------
    table : pandas.DataFrame
        Table of enzymes with desired columns returned as a DataFrame.

    Notes
    -----
    Based on PROSO toolbox `pcModel.m`.
    https://github.com/QCSB/PROSO-Toolbox/blob/main/src/PC-FBA/pcModel.m
    """

    # Any additional annotations or attributes that should be included. Can be used to add annotations and other attribute columns to objects.
    table_type = "enzymes"
    annotation_columns = (
        [] if not annotation_columns else ensure_iterable(annotation_columns)
    )
    if optional_columns is True:
        optional_columns = TABLE_COLUMNS[table_type]["optional"]
    elif optional_columns and set(optional_columns).difference(
        TABLE_COLUMNS[table_type]["optional"]
    ):
        raise ValueError(
            f"`optional_columns` must be one or more of the following: {TABLE_COLUMNS[table_type]['optional']}"
        )
    elif not optional_columns:
        optional_columns = []
    else:
        optional_columns = ensure_iterable(optional_columns)

    valid = {"reactions"}.union(annotation_columns)
    if id_key is not None and id_key not in valid:
        raise ValueError(f"`id_key` must be one of the following: {valid}")

    # Complex to reaction is necessary for correct enzyme mapping, can be obtained from complex table
    table = pd.DataFrame.from_dict(
        {
            idx: {"complexes": k, "reactions": build_string(v)}
            for idx, (k, v) in enumerate(complexes_to_reactions.items())
        },
        orient="index",
    )
    table = explode_column(table, "reactions")
    table["compartment"] = table["reactions"].apply(
        lambda x: build_string(sorted(model.reactions.get_by_id(x).compartments))
    )
    # table["complex_keff"] = table["complexes"].apply(
    #     lambda x: ";".join([f"{enzyme_keff_base}" for _ in range(len(split_string(x)))])
    # )

    # Identifiers
    table_type_strip = strip_plural(table_type)
    if not id_key:
        prefix = (
            prefix
            if (prefix and isinstance(prefix, str))
            else DEFAULT_PREFIX_SUFFIX_VALUES[table_type]["prefix.metabolite"]
        )
        id_key = {rid: i for i, rid in enumerate(sorted(table["reactions"].unique()))}
        table[table_type_strip] = [
            f"{prefix}{id_key[rid]}".replace("-", "_")
            for rid in table["reactions"].values
        ]
    else:
        if prefix:
            prefix = (
                prefix
                if (prefix and isinstance(prefix, str))
                else DEFAULT_PREFIX_SUFFIX_VALUES[table_type]["prefix.metabolite"]
            )
        else:
            prefix = ""
        table[table_type_strip] = table[id_key].apply(
            lambda x: f"{prefix}{x}".replace("-", "_")
        )
    if not enzyme_forward_suffix:
        enzyme_forward_suffix = DEFAULT_PREFIX_SUFFIX_VALUES["enzymes"]["suffix.forward"]
    table[f"{table_type_strip}{enzyme_forward_suffix}"] = table[table_type_strip].apply(
        lambda x: f"{x}{enzyme_forward_suffix}"
    )

    if not enzyme_reverse_suffix:
        enzyme_reverse_suffix = DEFAULT_PREFIX_SUFFIX_VALUES["enzymes"]["suffix.reverse"]
    table[f"{table_type_strip}{enzyme_reverse_suffix}"] = table[table_type_strip].apply(
        lambda x: f"{x}{enzyme_reverse_suffix}"
    )
    table = (
        table.melt(
            # Exclude "direction" from optional
            id_vars=[
                "compartment",
                "reactions",
                "complexes",  # "complex_keff"
            ],
            value_vars=[
                f"{table_type_strip}{enzyme_forward_suffix}",
                f"{table_type_strip}{enzyme_reverse_suffix}",
            ],
        )
        .rename({"value": table_type_strip}, axis=1)
        .drop("variable", axis=1)
    )
    # Now add direction
    table["direction"] = [
        "reverse" if x.endswith(enzyme_reverse_suffix) else "forward"
        for x in table[table_type_strip].values
    ]
    # # Needs to be a string `'0'` for aggregation step to work
    # table["complex_keff"] = table[["reactions", "complex_keff", "direction"]].apply(
    #     lambda x: "0"
    #     if (
    #         not model.reactions.get_by_id(x["reactions"]).reversibility
    #         and x["direction"] == "reverse"
    #     )
    #     else x["complex_keff"],
    #     axis=1,
    # )
    # # Aggregate data
    # table["complexes_keff"] = table.apply(
    #     lambda x: build_string(
    #         [
    #             f"{k}({v})"
    #             for k, v in create_variable_value_mapping(
    #                 *x[["complexes", "complex_keff"]].values
    #             ).items()
    #         ]
    #     ),
    #     axis=1,
    # )
    agg_mapping = {
        x: lambda x: build_string(x.dropna())
        for x in TABLE_COLUMNS[table_type]["required"][1:]
        + TABLE_COLUMNS[table_type]["optional"][1:]
    }
    # agg_mapping.update({"complex_keff": lambda x: ";".join(x)})
    table = table.groupby([table_type_strip]).agg(agg_mapping).reset_index(drop=False)

    # Clean up for consistency and add missing annotations
    required_columns = TABLE_COLUMNS[table_type]["required"]
    columns = required_columns + optional_columns + annotation_columns
    for col in columns:
        if table.get(col) is None:
            table[col] = pd.NA

    attr = "reactions"
    for col in annotation_columns:
        table[col] = table[attr].apply(
            lambda x: [
                getattr(model, attr).get_by_id(sid).annotation.get(col)
                for sid in split_string(x)
                if getattr(model, attr).has_id(sid)
                and getattr(model, attr).get_by_id(sid).annotation.get(col)
            ]
        )
        table[col] = table[col].apply(
            lambda x: build_string(
                sorted(
                    set([item for sublist in x for item in ensure_iterable(sublist)])
                )
            )
        )

    if replace_compartments:
        table["compartment"] = table["compartment"].replace(replace_compartments)
    table["compartment"] = table["compartment"].apply(
        lambda x: "".join(sorted(set(split_string(x))))
    )
    table[table_type] = table[[table_type_strip, "compartment"]].apply(
        lambda x: "_".join(x), axis=1
    )
    table = table.replace(float("nan"), pd.NA).replace("", pd.NA)
    table = (
        table.loc[:, columns]
        .drop_duplicates()
        .sort_values(required_columns)
        .reset_index(drop=True)
    )
    return table


def create_sequence_table(
    df_sequences, df_copy_numbers=None, mapping_key=None, isoform_transform=None
):
    """Format the protein data to use in advanced creation of the protein table.

    Parameters
    ----------
    df_sequences : pandas.DataFrame
        DataFrame of proteins mapped to sequences. Expected columns include the following:
            * ``sequence.id``: Identifier for the protein sequence. If None provided, will use the column given by the `mapping_key`.
            * ``sequence``: The amino acid sequence of the protein. If None provided, will use the default molar mass.
            * ``molar_mass`` : Molar mass of the protein sequence. Provided values remain while remaining values are calculated based on the sequences if provided or set as the default
            * The additional column used for merging DataFrames as defined by the ``mapping_key``.
    df_copy_numbers : pandas.DataFrame
        DataFrame of proteins mapped to copy numbers. Expected columns include the following:
            * ``copy_number`` : The copy number of the protein. If None provided, fill with the value of 0.
            * The additional column used for merging DataFrames as defined by the ``mapping_key``.
    mapping_key : str
        The main column to use for merging DataFrames.
    isoform_transform : dict, bool
        If True, the mean value for is used to transform the molar mass and copy number columns for isoforms with identical values on the `mapping_key`.
        If False, no transform is applied to isoforms. Use this option to preserved isoforms to apply additional constraints later.
        Alternatively, a dictionary mapping colunns the `molar_mass` and/or `copy_number` columns to functions for aggregating data.
    """
    # Cannot add proteins that do not have molar mass included. Furthermore, cannot usually distingusih copy number of isoforms.
    # Therefore, sequences take priority in mapping, and any isoforms need to be addressed with additional constraints.
    df_protein_data = df_sequences.copy()
    if df_protein_data.get("sequence.id") is None:
        df_protein_data["sequence.id"] = df_protein_data[mapping_key]

    # If no sequence or molar mass defined, use the default value.
    if (
        df_protein_data.get("molar_mass") is None
        and df_protein_data.get("sequence") is None
    ):
        # Set as NA to fill with default molar mass
        df_protein_data["molar_mass"] = pd.NA
        df_protein_data["molar_mass_from"] = pd.NA
    # If sequence is defined and molar mass is not defined, calculate molar mass
    elif (
        df_protein_data.get("molar_mass") is None
        and df_protein_data.get("sequence") is not None
    ):
        df_protein_data["molar_mass"] = (
            df_protein_data["sequence"]
            .fillna("")
            .apply(lambda x: calculate_protein_molar_mass(x) if x else pd.NA)
        )
        df_protein_data["molar_mass_from"] = (
            df_protein_data["molar_mass"]
            .notna()
            .apply(lambda x: "SEQUENCE" if x else pd.NA)
        )
    else:
        # If molar mass were already defined, assume molar mass was pre-calculated.
        df_protein_data["molar_mass_from"] = (
            df_protein_data["molar_mass"]
            .notna()
            .apply(lambda x: "PROVIDED" if x else pd.NA)
        )

    df_protein_data["molar_mass"] = df_protein_data["molar_mass"].fillna(
        DEFAULT_PROTEIN_MOLAR_MASS
    )
    df_protein_data["molar_mass_from"] = df_protein_data["molar_mass_from"].fillna(
        "DEFAULT"
    )

    # Add copy numbers
    if df_copy_numbers is not None:
        df_protein_data = df_protein_data.merge(
            df_copy_numbers, left_on=mapping_key, right_on=mapping_key, how="left"
        )
        # # First, set copy number of 0 for any proteins that have sequences included but not a copy number:
        df_protein_data["copy_number"] = df_protein_data["copy_number"].fillna(0)
        df_protein_data["copy_number_from"] = df_protein_data["copy_number"].apply(
            lambda x: "MEASURED" if x else pd.NA
        )
        df_protein_data["copy_number_from"] = df_protein_data[
            "copy_number_from"
        ].fillna("NONE")

    df_protein_data = df_protein_data.convert_dtypes().infer_objects()
    df_isoforms = df_protein_data[df_protein_data[mapping_key].duplicated(False)].copy()
    df_protein_data = df_protein_data[~df_protein_data[mapping_key].duplicated(False)]
    # Transform the isforms to averages unless there are plans for special contraints on isoforms
    if isoform_transform:
        if not isinstance(isoform_transform, dict):
            isoform_transform = {"molar_mass": "mean"}
            if df_copy_numbers is not None:
                isoform_transform.update({"copy_number": "mean"})

        agg_mapping = {
            col: lambda x: build_string(x)
            for col in df_isoforms.columns
            if col not in ["molar_mass", "copy_number"]
        }
        agg_mapping.update(
            {
                "molar_mass": isoform_transform.get("molar_mass", "mean"),
                "molar_mass_from": lambda x: "TRANSFORMED",
            }
        )
        if df_copy_numbers is not None:
            agg_mapping.update(
                {
                    "copy_number": isoform_transform.get("copy_number", "mean"),
                    "copy_number_from": lambda x: "TRANSFORMED",
                }
            )
        df_isoforms = df_isoforms.groupby(mapping_key, as_index=False).agg(agg_mapping)
    df_protein_data = (
        pd.concat((df_protein_data, df_isoforms), axis=0)
        .sort_values(mapping_key)
        .reset_index(drop=True)
    )
    return df_protein_data


def add_dilution_reaction(
    pcmodel, protein, protein_type=None, bounds=None, gene_reaction_rule="", verbose=0
):
    """Add a pseudoreaction representing protein dilution/accumulation in the model.

    Parameters
    ---------

    Returns
    -------

    """
    if protein_type is None:
        protein_type = protein.__class__.__name__.lower()
    valid = {"protein", "complex", "enzyme"}
    if protein_type not in valid:
        raise ValueError(f"`protein_type` must be one of the following: {valid}.")

    cls, prefix, default_bounds = {
        "protein": (
            ProteinDilution,
            DEFAULT_PREFIX_SUFFIX_VALUES["proteins"]["prefix.dilution"],
            (0, DEFAULT_CONCENTRATION_BOUND)
        ),
        "complex": (
            ComplexDilution,
            DEFAULT_PREFIX_SUFFIX_VALUES["complexes"]["prefix.dilution"],
            (0, DEFAULT_CONCENTRATION_BOUND),
        ),
        "enzyme": (
            EnzymeDilution,
            DEFAULT_PREFIX_SUFFIX_VALUES["enzymes"]["prefix.dilution"],
            (0, DEFAULT_CONCENTRATION_BOUND),
        ),
    }[protein_type]

    if bounds is not None:
        bounds = ensure_iterable(bounds)
        bounds = tuple(
            [
                (
                    float(concentration)
                    if isinstance(concentration, str)
                    else concentration
                )
                for concentration in bounds
            ]
        )
        if any([bound < 0 for bound in bounds]):
            raise ValueError("`bounds` must be a tuple of non-negative numbers.")
    else:
        bounds = default_bounds

    protein = pcmodel.metabolites.get_by_id(str(protein))
    rid = f"{prefix}{protein.id}"
    if not pcmodel.reactions.has_id(rid):
        dilution_rxn = cls(
            id=rid,
            name=f"{protein.name if protein.name else protein.id} concentration",
            subsystem=f"Pseudoreactions, {protein_type.capitalize()} concentrations",
        )
        pcmodel.add_reactions([dilution_rxn])
    else:
        log_msg(
            LOGGER,
            logging.WARNING,
            "%s %s already exists, updating with table values.",
            cls.__name__,
            rid,
            print_lvl=np.floor(10 * verbose),
        )
    dilution_rxn = pcmodel.reactions.get_by_id(rid)
    # Set to default concentrations. COBRA default bounds are typically 1000 mmol / hr / gDW
    # However for protein and enzyme pseudoreactions, the "flux" is not real, but actually represents the concentration in nmol / gDW / cell
    # All proteins/enzymes are initialized at 100000 nmol / gDW / cell (or equivalently, 1 mmol / gDW / cell)
    # If initialized as "EX_protein_i: protein_i <-- ", then concentration given as: protein_i = -v_{EX_protein_i}
    # If initialized as "EX_protein_i: --> protein_i"., then concentration given as: protein_i = v_{EX_protein_i}
    # For simplicity, initialize all proteins to have positive values. Complexes and enzymes are already formatted to be positive.
    sign = 1 if protein_type == "protein" else -1
    dilution_rxn.add_metabolites({protein: sign * 1}, combine=False)
    dilution_rxn.bounds = bounds
    if gene_reaction_rule:
        dilution_rxn.gene_reaction_rule = gene_reaction_rule
    # Relaxation on constraint for protein metabolite
    # protein.constraint.lb = 0
    # protein.constraint.ub = float(concentration)

    return dilution_rxn


def add_complex_formation_reaction(
    pcmodel,
    cplx,
    complex_type,
    coeff_map=None,
    gene_reaction_rule="",
    verbose=0,
    bounds=None,
):
    """Add a reaction representing the formation of a protein complex to model.

    Parameters
    ----------

    Returns
    -------

    """
    valid = {"complex", "enzyme"}
    if complex_type not in valid:
        raise ValueError(f"`complex_type` must be one of the following: {valid}.")

    cls, prefix, default_bounds = {
        "complex": (
            ComplexForm,
            DEFAULT_PREFIX_SUFFIX_VALUES["complexes"]["prefix.formation"],
            (0, DEFAULT_CONCENTRATION_BOUND),
        ),
        "enzyme": (
            EnzymeForm,
            DEFAULT_PREFIX_SUFFIX_VALUES["enzymes"]["prefix.formation"],
            (0, DEFAULT_CONCENTRATION_BOUND),
        ),
    }[complex_type]
    coeff_mapping_attr = (
        "complexes_keff" if complex_type == "enzyme" else "stoichiometry"
    )

    # Add complex formation reactions
    cplx = pcmodel.metabolites.get_by_id(str(cplx))

    # 2.3 Set stoichiometry
    base_val = 1
    if coeff_map:
        if "(" and ")" in coeff_map:
            coeff_map = create_variable_value_mapping(coeff_map)
        else:
            log_msg(
                LOGGER,
                logging.INFO,
                "No coefficients provided for '%s', assigning base value of '%s'",
                coeff_mapping_attr,
                base_val,
                print_lvl=verbose,
            )
            coeff_map = {variable: base_val for variable in split_string(coeff_map)}
    else:
        warn(
            f"No {coeff_mapping_attr} provided for {complex_type.capitalize()} `{formation_rxn.id}`, unable to set {cls.__name__} `{formation_rxn.id}` {coeff_mapping_attr}."
        )
        coeff_map = {}

    rid = f"{prefix}{cplx.id}"
    if complex_type == "enzyme":
        rid = f"{rid.replace(f'_{cplx.compartment}', '')}_{list(coeff_map)[0]}"
    if not pcmodel.reactions.has_id(rid):
        formation_rxn = cls(
            id=rid,
            name=f"{cplx.name if cplx.name else cplx.id} concentration",
            subsystem=f"Pseudoreactions, {complex_type.capitalize()} concentrations",
        )
        pcmodel.add_reactions([formation_rxn])
    else:
        log_msg(
            LOGGER,
            logging.WARNING,
            "%s %s already exists, updating with table values.",
            cls.__name__,
            rid,
            print_lvl=np.floor(10 * verbose),
        )
    formation_rxn = pcmodel.reactions.get_by_id(rid)

    # Only add metabolites if proteins are asociated, otherwise creates a sink which may create errors
    coeff_map = {
        pcmodel.metabolites.get_by_id(sid): -abs(float(coeff))
        for sid, coeff in coeff_map.items()
    }
    coeff_map.update({cplx: 1})
    formation_rxn.add_metabolites(coeff_map, combine=False)
    # Set to default concentrations. COBRA default bounds are typically 1000 mmol / hr / gDW
    # However for complex formation reactions, the "flux" is not real, but actually represents the concentration in nmol / gDW / cell
    # All complexes have an upper limit initialized at 100000 nmol / gDW / cell (or equivalently, 1 mmol / gDW / cell)
    # If initialized as "CPLXFORM_cplx_j: protein_i1 + protein_i2 <=> cplx_j", then concentration given as: cplx_j = v_{CPLXFORM_cplx_j}
    # If initialized as "CPLXFORM_cplx_j: cplx_j <=> protein_i1 + protein_i2", then concentration given as: cplx_j = -v_{CPLXFORM_cplx_j}

    #  print({
    #         met: (abs(formation_rxn.metabolites[met]), max(list(map(abs, reaction.bounds))))
    #         for met in formation_rxn.reactants
    #         for reaction in list(met.reactions)
    #         if reaction.id.endswith(met.id)
    # })
    if bounds is not None:
        bounds = ensure_iterable(bounds)
        bounds = tuple(
            [
                (
                    float(concentration)
                    if isinstance(concentration, str)
                    else concentration
                )
                for concentration in bounds
            ]
        )
        if any([bound < 0 for bound in bounds]):
            raise ValueError("`bounds` must be a tuple of non-negative numbers.")
    else:
        bounds = default_bounds
    if gene_reaction_rule:
        formation_rxn.gene_reaction_rule = gene_reaction_rule
    return formation_rxn


def construct_pcmodel_from_tables(
    model,
    protein_table,
    complex_table,
    enzyme_table,
    max_weight_fraction=DEFAULT_MAX_WEIGHT_FRACTION,
    enzyme_keff_base=DEFAULT_KEFF,
    enzyme_forward_suffix=None,
    enzyme_reverse_suffix=None,
    enzyme_total_suffix=None,
    include_complex_dilutions=False,
    irrev_rxn_complex_keff=None,
    verbose=0,
):
    """Construct a protein-constrained model using the OVERLAY method.

    Parameters
    ----------

    Returns
    -------

    """
    # 0. Parse input
    if max_weight_fraction > 500 or max_weight_fraction < 50:
        if max_weight_fraction <= 0:
            raise ValueError("`max_weight_fraction` must be a positive non-zero value.")
        msg = "over 50" if max_weight_fraction > 500 else "under 5"
        warn(f"Maximum proteome weight fraction is {msg}% (={max_weight_fraction}).")

    pcmodel = model.copy()
    pcmodel.id += "_PC"
    tables = {}
    add_table_cols = defaultdict(dict)
    direction_dict = {
        "forward": (DEFAULT_PREFIX_SUFFIX_VALUES["enzymes"]["suffix.forward"] if not enzyme_forward_suffix else enzyme_forward_suffix, -1),
        "reverse": (DEFAULT_PREFIX_SUFFIX_VALUES["enzymes"]["suffix.reverse"] if not enzyme_reverse_suffix else enzyme_reverse_suffix, 1),
    }
    for table_type, table, additional in zip(
        ["proteins", "complexes", "enzymes"],
        [protein_table, complex_table, enzyme_table],
        [
            None,
            ("stoichiometry", ["subunits", "coefficients"]),
            None,
            # ("complexes_keff", ["complexes", "complex_keff"]),
        ],
    ):
        table = _format_table_input(table, table_type=table_type, additional=additional)
        # Special for enzyme direction
        if table_type == "enzymes":
            if table.get("reactions") is None:
                raise KeyError(
                    "`enzyme_table` must have a column containing the associated reactions."
                )

            if table.get("direction") is None:
                log_msg(
                    LOGGER,
                    logging.INFO,
                    "Attempting to determine enzyme direction from forward suffix `%s` and revere suffix  `%s`.",
                    enzyme_forward_suffix,
                    enzyme_reverse_suffix,
                    print_lvl=verbose,
                )

                table["direction"] = table["enzyme"].apply(
                    lambda x: [
                        k for k, v in direction_dict.items() if v[0] == x.split("_")[-1]
                    ].pop()
                )
        tables[table_type] = table

    # 1. Initialize proteome budget constraint and conversion factor
    pcmodel.add_metabolites(
        [
            ProteomeBudget(
                id="proteome_budget",
                name="Proteome Budget Constraint",
                compartment="c",  # In RBCs, the cytosol is the only compartment
            )
        ]
    )
    proteome_budget = pcmodel.metabolites.get_by_id("proteome_budget")
    pbid = f"PBDL_{proteome_budget.id}"
    if not pcmodel.reactions.has_id(pbid):
        pcmodel.add_reactions(
            [ProteomeBudgetDilution(id=pbid, name="Proteome budget demand")]
        )
    budget_demand = pcmodel.reactions.get_by_id(pbid)
    budget_demand.add_metabolites({proteome_budget: -1}, combine=False)
    budget_demand.bounds = (0, max_weight_fraction)

    # Add as a boundary reaction instead of setting constraint bounds, can be updated later.

    cf = 1 / 1e6  # Conversion factor from nmol to mmol

    # 2. Add proteins to model
    table_type = "proteins"
    cls = ATTR_SUBCLASS_DICT["metabolites"][table_type]
    log_msg(
        LOGGER,
        logging.INFO,
        "Adding %s to model.",
        table_type,
        print_lvl=verbose,
    )
    try:
        for sid, row in tables[table_type].set_index(table_type).iterrows():
            # 2.1 Create objects
            # ID and compartment are required by nearly all modeling standards.
            if not pcmodel.metabolites.has_id(sid):
                item = cls(sid, compartment=row.get("compartment"))
                pcmodel.add_metabolites([item])
            else:
                log_msg(
                    LOGGER,
                    logging.WARNING,
                    "%s %s already exists, updating with table values.",
                    cls.__name__,
                    sid,
                    print_lvl=verbose,
                )
            item = pcmodel.metabolites.get_by_id(sid)
            # Initialize entry for item
            add_table_cols[table_type][item.id] = {}
            # 2.2 Set inherited attributes (if any)
            for attr in ATTRIBUTE_MAPPING["metabolites"]:
                value = row.get(attr)
                if value is None or value == "":
                    continue
                setattr(item, attr, value)
                add_table_cols[table_type][item.id].update({attr: value})

            # 2.3 Get molar mass from sequence or provided.
            # 2.3.1 Protein sequences
            sequences = row.get("sequence", {})
            if sequences:
                if row.get("sequence.id", ""):
                    iter_seq = zip(
                        split_string(row["sequence.id"]), split_string(sequences)
                    )
                else:
                    iter_seq = enumerate(split_string(sequences))
                sequences = {seq_id: seq for seq_id, seq in iter_seq}
            else:
                sequences = {}
            # 2.3.2 Molar mass
            # Always use the molar mass column if provided.
            molar_mass = row.get("molar_mass")
            if not molar_mass:
                # If no mass column exists, try computing from sequence(s).
                # Multiple sequences are seperated by semicolons and the average molar mass is used.
                if sequences:
                    log_msg(
                        LOGGER,
                        logging.DEBUG,
                        "Calculating molar mass from amino acid sequence for %s `%s`.",
                        cls.__name__,
                        sid,
                        print_lvl=verbose,
                    )
                    molar_mass = np.mean(
                        [
                            calculate_protein_molar_mass(seq)
                            for seq in sequences.values()
                        ]
                    )
                else:
                    log_msg(
                        LOGGER,
                        logging.WARNING,
                        "No molar mass or amino acid sequence provided for %s `%s`, using default instead.",
                        cls.__name__,
                        item.id,
                        print_lvl=verbose,
                    )
                    molar_mass = DEFAULT_PROTEIN_MOLAR_MASS
            molar_mass = float(molar_mass)
            # 2.4 Add pseudoreactions representing concentration/dilution (if any)
            # Pseudoreactions use the "flux bounds" to be representative of concentrations in nmol / gDW
            # Concetrations are always positive, therefore pseudoreactions are irreversible in the direction that facilitates a positive concentration.
            bounds = (
                (
                    float(row.get("lower_bound"))
                    if row.get("lower_bound")
                    else 0
                ),
                (
                    float(row.get("upper_bound"))
                    if row.get("upper_bound")
                    else DEFAULT_CONCENTRATION_BOUND
                ),
            )
            gene_reaction_rule = (
                " and ".join(row["genes"].split(";")) if row.get("genes") else ""
            )
            dilution_rxn = add_dilution_reaction(
                pcmodel,
                item,
                strip_plural(table_type),
                bounds=bounds,
                gene_reaction_rule=gene_reaction_rule,
            )
            # 2.5 Add proteome budget constraint using molar mass of proteins
            # mg / mmol (=g / mol) * (1 mmol / 1e6 nmol) --> mg / nmol
            dilution_rxn.add_metabolites(
                {proteome_budget: molar_mass * cf}, combine=True
            )

            # Concentration should always be positive
            add_table_cols[table_type][item.id].update(
                {
                    "sequence.id": build_string(list(sequences.keys())),
                    "sequence": build_string(list(sequences.values())),
                    "molar_mass": molar_mass,
                    "lower_bound": dilution_rxn.lower_bound,
                    "upper_bound": dilution_rxn.upper_bound,
                }
            )
    except Exception as e:
        print(f"Problem with row {sid} in {table_type} table")
        raise e

    # 3. Add complexes to model
    table_type = "complexes"
    cls = ATTR_SUBCLASS_DICT["metabolites"][table_type]
    log_msg(
        LOGGER,
        logging.INFO,
        "Adding %s to model.",
        table_type,
        print_lvl=verbose,
    )
    try:
        for sid, row in tables[table_type].set_index(table_type).iterrows():
            # 3.1 Create objects
            # ID and compartment are required by nearly all modeling standards.
            if not pcmodel.metabolites.has_id(sid):
                item = cls(sid, compartment=row.get("compartment"))
                pcmodel.add_metabolites([item])
            else:
                log_msg(
                    LOGGER,
                    logging.WARNING,
                    "%s %s already exists, updating with table values.",
                    cls.__name__,
                    sid,
                    print_lvl=verbose,
                )
            item = pcmodel.metabolites.get_by_id(sid)
            # Initialize entry for item
            add_table_cols[table_type][item.id] = {}
            # 3.2 Set inherited attributes (if any)
            for attr in ATTRIBUTE_MAPPING["metabolites"]:
                value = row.get(attr)
                if value is None or value == "":
                    continue
                setattr(item, attr, value)
                add_table_cols[table_type][item.id].update({attr: value})

            # 3.3 Add formation reaction (if any)
            gene_reaction_rule = (
                " and ".join(row["genes"].split(";")) if row.get("genes") else ""
            )
            formation_rxn = add_complex_formation_reaction(
                pcmodel,
                item,
                strip_plural(table_type),
                coeff_map=row["stoichiometry"],
                gene_reaction_rule=gene_reaction_rule,
            )
            bounds = (
                (
                    float(row.get("lower_bound"))
                    if row.get("lower_bound")
                    else 0
                ),
                (
                    float(row.get("upper_bound"))
                    if row.get("upper_bound")
                    else DEFAULT_CONCENTRATION_BOUND
                ),
            )
            formation_rxn.bounds = bounds
            if row.get("molar_mass") is None:
                molar_mass = sum(
                    [
                        abs(coeff) * add_table_cols["proteins"][protein]["molar_mass"]
                        for protein, coeff in create_variable_value_mapping(
                            row["stoichiometry"]
                        ).items()
                    ]
                )
                add_table_cols[table_type][item.id].update(
                    {
                        "molar_mass": molar_mass,
                    }
                )
            # print({
            #     met: (abs(formation_rxn.metabolites[met]), "protein_mass")
            #     for met in formation_rxn.reactants
            #     for reaction in list(met.reactions)
            #     if reaction.id.endswith(met.id)
            # })
            # 3.3 Add pseudoreactions representing concentration/dilution (if any)
            # Pseudoreactions use the "flux bounds" to be representative of concentrations in nmol / gDW / cell.
            # Concetrations are always positive, therefore pseudoreactions are irreversible in the direction that facilitates a positive concentration.
            # In the case of complexes, a dilution reaction provides flexibility in the constraints
            if include_complex_dilutions:
                dilution_rxn = add_dilution_reaction(
                    pcmodel,
                    item,
                    strip_plural(table_type),
                    bounds=formation_rxn.bounds,
                    gene_reaction_rule=gene_reaction_rule,
                )
                # add_table_cols[table_type][item.id].update({
                #     "lower_bound": dilution_rxn.lower_bound,
                #     "upper_bound": dilution_rxn.upper_bound,
                # })
    except Exception as e:
        print(f"Problem with row {sid} in {table_type} table")
        raise e

    df = pd.DataFrame.from_dict(add_table_cols["complexes"], orient="index")
    df = enzyme_keff_base * (df["molar_mass"] / df["molar_mass"].mean()) ** 0.75
    keff_from_molar_mass = df.to_dict()
    # Rescale based on new average
    mean_val = np.mean(list(keff_from_molar_mass.values()))
    keff_from_molar_mass = {
        cplx: cplx_keff * enzyme_keff_base / mean_val
        for cplx, cplx_keff in keff_from_molar_mass.items()
    }

    # 4. Add enzymes to model
    table_type = "enzymes"
    cls = ATTR_SUBCLASS_DICT["metabolites"][table_type]
    log_msg(
        LOGGER,
        logging.INFO,
        "Adding %s to model.",
        table_type,
        print_lvl=verbose,
    )
    cplx_prefix = DEFAULT_PREFIX_SUFFIX_VALUES["complexes"]["prefix.formation"]
    try:
        for sid, row in tables[table_type].set_index(table_type).iterrows():
            # 4.1 Create objects
            # ID and compartment are required by nearly all modeling standards.
            if not pcmodel.metabolites.has_id(sid):
                item = cls(sid, compartment=row.get("compartment"))
                pcmodel.add_metabolites([item])
            else:
                log_msg(
                    LOGGER,
                    logging.WARNING,
                    "%s %s already exists, updating with table values.",
                    cls.__name__,
                    sid,
                    print_lvl=verbose,
                )
            item = pcmodel.metabolites.get_by_id(sid)
            # Initialize entry for item
            add_table_cols[table_type][item.id] = {}
            # 4.2 Set inherited attributes (if any)
            for attr in ATTRIBUTE_MAPPING["metabolites"]:
                value = row.get(attr)
                if value is None or value == "":
                    continue
                setattr(item, attr, value)
                add_table_cols[table_type][item.id].update({attr: value})
            reaction, direction = row[["reactions", "direction"]].values
            reaction = pcmodel.reactions.get_by_id(reaction)

            # 4.3 Add formation reactions (if any)
            # 4.3.1 Parse complexes and keffs to determine number of formation reactions to add.
            coeff_map = row.get("complexes_keff", row["complexes"])
            if "(" and ")" in coeff_map:
                # Mapping with coefficients exists
                coeff_map = create_variable_value_mapping(coeff_map)
            else:
                # No coefficients, only variables provided
                # Set up dictionary with variables and default keffs
                coeff_map = {
                    variable: keff_from_molar_mass[variable]
                    for variable in split_string(coeff_map)
                }
                if (
                    not reaction.reversibility and direction == "reverse"
                ) and irrev_rxn_complex_keff is not None:
                    coeff_map.update(
                        {
                            variable: irrev_rxn_complex_keff
                            for variable in coeff_map.keys()
                        }
                    )

            # If an enzyme keff is provided, then coefficients are scaled against it.
            # A coefficient of 0 means that the complex should not be included with a enzyme formation reaction.
            enzyme_keff = float(row.get("enzyme_keff", enzyme_keff_base))
            for cplx, keff in coeff_map.items():
                if keff == 0 or enzyme_keff == 0:
                    continue

                keff = float(keff) / float(enzyme_keff)
                cplx_formation_rxn = pcmodel.reactions.get_by_id(f"{cplx_prefix}{cplx}")
                formation_rxn = add_complex_formation_reaction(
                    pcmodel,
                    item,
                    complex_type=strip_plural(table_type),
                    coeff_map=f"{cplx}({keff})",
                    bounds=cplx_formation_rxn.bounds,
                    gene_reaction_rule=cplx_formation_rxn.gene_reaction_rule,
                )
            # 4.4 Add pseudoreactions representing concentration/dilution (if any)
            # Pseudoreactions use the "flux bounds" to be representative of concentrations in nmol / gDW / cell.
            # Concetrations are always positive, therefore pseudoreactions are irreversible in the direction that facilitates a positive concentration.
            dilution_rxn = add_dilution_reaction(
                pcmodel,
                item,
                strip_plural(table_type),
                gene_reaction_rule=reaction.gene_reaction_rule,
            )

            # Summation variable
            if not enzyme_total_suffix:
                enzyme_total_suffix = DEFAULT_PREFIX_SUFFIX_VALUES["enzymes"]["suffix.total"]
            sum_var = item.id.replace(
                f"{direction_dict[direction][0]}_", f"{enzyme_total_suffix}_"
            )
            if not pcmodel.metabolites.has_id(sum_var):
                sum_var = cls(sum_var, compartment=item.compartment)
                pcmodel.add_metabolites([sum_var])
                add_dilution_reaction(
                    pcmodel,
                    sum_var,
                    strip_plural(table_type),
                    bounds=formation_rxn.bounds,
                    gene_reaction_rule=reaction.gene_reaction_rule,
                )
            sum_var = pcmodel.metabolites.get_by_id(str(sum_var))
            dilution_rxn.add_metabolites({sum_var: 1})
            # 4.5. Add proteomic constraints to model reactions
            # Constrain each enzyme by the base rate constant.
            # 1 / (1 / hr) / (mmol / nmol) --> (hr / (mmol/nmol)) --> hr / mmol
            sign = direction_dict[direction][1]
            reaction.add_metabolites(
                {item: sign * (1 / enzyme_keff / cf)}, combine=True
            )
            # Systemic effective rate constant
            add_table_cols[table_type][item.id].update(
                {
                    "complexes_keff": build_string(
                        [f"{cplx}({keff})" for cplx, keff in coeff_map.items()]
                    ),
                    "complex_keff": ";".join(
                        ([str(keff) for keff in coeff_map.values()])
                    ),
                    "enzyme_keff": enzyme_keff,
                    # "lower_bound": dilution_rxn.lower_bound,
                    # "upper_bound": dilution_rxn.upper_bound,
                }
            )
    except Exception as e:
        print(f"Problem with row {sid} in {table_type} table")
        raise e

    # Add additional columns to table, replacing previous ones
    for table_type, table in tables.items():
        table = table.merge(
            pd.DataFrame.from_dict(add_table_cols[table_type], orient="index"),
            left_on=table_type,
            right_index=True,
            how="left",
            suffixes=("_original", ""),
        )
        # Remove original
        table = table.drop(
            [col for col in table.columns if col.endswith("_original")], axis=1
        )
        # Reset index to make table_type initial column
        table = table.set_index(table_type).reset_index(drop=False)
        table = table.replace("", pd.NA).replace(float("nan"), pd.NA)
        # Replace with updated table
        tables[table_type] = table

    pcmodel, _ = manipulation.prune_unused_metabolites(pcmodel)
    pcmodel, _ = manipulation.prune_unused_reactions(pcmodel)
    pcmodel.repair()
    return pcmodel, tables


def _format_table_input(table, table_type, additional=None):
    """Format the table input."""
    # Make a copy
    formatted_table = table.copy().fillna("")
    table_type_strip = strip_plural(table_type)
    # Ensure compartment column exists, always essential
    if formatted_table.get("compartment") is None:
        raise KeyError(
            f"`{table_type}_table` must have a column `compartment` containing compartment identifiers."
        )

    # If no ID column found, generate one using the "protein" and "compartment" columns
    if formatted_table.get(table_type) is None:
        try:
            formatted_table[table_type] = (
                formatted_table[[table_type_strip, "compartment"]]
                .apply(lambda x: "_".join(x), axis=1)
                .values
            )
        except KeyError as e:
            if formatted_table.get("id") is None:
                raise KeyError(
                    f"`{table_type}_table` must have a column of {table_type_strip} identifiers. "
                    f"Either a column of {table_type_strip} identifiers named `{table_type}` or `id` must exist, or"
                    f"the `{table_type_strip}` and `compartment` columns must exist in order to generate an ID."
                )
            formatted_table[table_type] = formatted_table["id"].values

    if additional:
        key, columns = additional[0], additional[1]
        if formatted_table.get(key) is None:
            try:
                formatted_table[key] = formatted_table[columns].apply(
                    lambda x: ";".join(
                        [
                            f"{k}({v})"
                            for k, v in dict(zip(*x.str.split(";").values)).items()
                        ]
                    ),
                    axis=1,
                )
            except KeyError as e:
                raise KeyError(
                    f"`{table_type}_table` must have a column for `{key}`. "
                    f"Either a column named `{key}` must exist, or"
                    f"the `{columns[0]}` and `{columns[1]}` columns must exist in order to generate the `{key}`."
                )

    return formatted_table


def compute_complex_molar_mass(protein_table, complex_table, verbose=0):
    """Compute the molar mass of complexes"""
    protein_table = _format_table_input(
        protein_table, table_type="proteins", additional=None
    )
    complex_table = _format_table_input(
        complex_table,
        table_type="complexes",
        additional=("stoichiometry", ["subunits", "coefficients"]),
    )

    if protein_table.get("molar_mass") is not None:
        log_msg(
            LOGGER,
            logging.INFO,
            "Using given molar mass values for proteins.",
            print_lvl=verbose,
        )
    elif (
        protein_table.get("molar_mass") is None
        and protein_table.get("sequence") is not None
    ):
        log_msg(
            LOGGER,
            logging.DEBUG,
            "Calculating molar mass from amino acid sequences.",
            print_lvl=verbose,
        )
        protein_table["molar_mass"] = protein_table["sequence"].apply(
            calculate_protein_molar_mass
        )
    else:
        raise KeyError(
            "Must provide a column of molar masses or protein amino acid sequences."
        )

    protein_to_mass_dict = protein_table.set_index("proteins")["molar_mass"].to_dict()
    complex_table["molar_mass"] = complex_table["stoichiometry"].apply(
        lambda x: sum(
            [
                coeff * protein_to_mass_dict[protein]
                for protein, coeff in create_variable_value_mapping(x).items()
            ]
        )
    )
    series = complex_table.set_index("complexes")["molar_mass"]
    return series


def estimate_keff_from_molar_mass(
    protein_table, complex_table, enzyme_keff_base=DEFAULT_KEFF, verbose=0
):
    """Estimate the effective rate constants for the given complexes."""
    series = compute_complex_molar_mass(protein_table, complex_table, verbose=verbose)
    series = enzyme_keff_base * (series / series.mean()) ** 0.75
    series = series.apply(lambda x: x * enzyme_keff_base / series.mean())
    return series


def load_overlay_model(
    filename,
    protein_prefixes=None,
    complex_prefixes=None,
    enzyme_prefixes=None,
    budget_id="budget",
    proteome_budget_prefix="PBDL_",
    relaxation_prefix="RELAX_",
):
    model = read_cobra_model(filename)
    if protein_prefixes is None:
        protein_met_prefix = DEFAULT_PREFIX_SUFFIX_VALUES["proteins"]["prefix.metabolite"]
        protein_dil_prefix = DEFAULT_PREFIX_SUFFIX_VALUES["proteins"]["prefix.dilution"]
    if complex_prefixes is None:
        complex_met_prefix = DEFAULT_PREFIX_SUFFIX_VALUES["complexes"]["prefix.metabolite"]
        complex_form_prefix = DEFAULT_PREFIX_SUFFIX_VALUES["complexes"]["prefix.formation"]
        complex_dil_prefix = DEFAULT_PREFIX_SUFFIX_VALUES["complexes"]["prefix.dilution"]
    if enzyme_prefixes is None:
        enzyme_met_prefix = DEFAULT_PREFIX_SUFFIX_VALUES["enzymes"]["prefix.metabolite"]
        enzyme_form_prefix = DEFAULT_PREFIX_SUFFIX_VALUES["enzymes"]["prefix.formation"]
        enzyme_dil_prefix = DEFAULT_PREFIX_SUFFIX_VALUES["enzymes"]["prefix.dilution"]

    for prefix, cls in zip(
        [protein_met_prefix, complex_met_prefix, enzyme_met_prefix],
        [Protein, Complex, Enzyme],
    ):
        for old in model.metabolites.query(lambda x: x.id.startswith(prefix)):
            new = cls(old)
            new.__dict__.update(old.__dict__)
            model.metabolites._replace_on_id(new)
            new.annotation.update(old.annotation)

    for old in model.metabolites.get_by_any(
        ensure_iterable([x for x in model.metabolites if budget_id in x.id])
    ):
        new = ProteomeBudget(old)
        new.__dict__.update(old.__dict__)
        model.metabolites._replace_on_id(new)
        new.annotation.update(old.annotation)

    for prefix, cls in zip(
        [
            protein_dil_prefix,
            complex_dil_prefix,
            enzyme_dil_prefix,
            complex_form_prefix,
            enzyme_form_prefix,
            proteome_budget_prefix,
            relaxation_prefix,
        ],
        [
            ProteinDilution,
            ComplexDilution,
            EnzymeDilution,
            ComplexForm,
            EnzymeForm,
            ProteomeBudgetDilution,
            ProteinDilution,
        ],
    ):
        if not prefix:
            continue

        for old in model.reactions.query(lambda x: x.id.startswith(prefix)):
            new = cls(old)
            new.__dict__.update(old.__dict__)
            model.reactions._replace_on_id(new)
            new._metabolites = {
                model.metabolites.get_by_id(m.id): coeff
                for m, coeff in new.metabolites.items()
            }
            new.annotation.update(old.annotation)

    model.repair()
    return model


def add_relaxation_budget(pcmodel, slack_value, verbose=False):
    # Get budget metabolites
    proteome_budget = pcmodel.metabolites.get_by_id("proteome_budget")
    hemoglobin_budget = pcmodel.metabolites.get_by_id("hemoglobin_budget")
    total_budget = pcmodel.metabolites.get_by_id("total_budget")
    # Add a relaxation budget to the model
    pcmodel.add_metabolites(
        [
            ProteomeBudget(
                id="relaxation_budget",
                name="Relaxation Budget Constraint",
                compartment=proteome_budget.compartment,
            )
        ]
    )
    relaxation_budget = pcmodel.metabolites.get_by_id("relaxation_budget")
    pbid = f"PBDL_{relaxation_budget.id}"
    pcmodel.add_reactions(
        [
            ProteomeBudgetDilution(
                id=pbid,
                name="Relaxation budget demand",
                lower_bound=0,
            )
        ]
    )
    relaxation_demand = pcmodel.reactions.get_by_id(pbid)
    relaxation_demand.add_metabolites(
        {relaxation_budget: -1, total_budget: 1}, combine=False
    )
    relax_to_protmw_dict = {}
    total_slack_mg_prot_per_gDW = 0
    hb_slack_mg_prot_per_gDW = 0
    for protdl in pcmodel.reactions.query(lambda x: isinstance(x, ProteinDilution)):
        protein = [m for m in protdl.metabolites if "budget" not in m.id][0]
        protmw = (
            calculate_protein_molar_mass(protdl.annotation.get("uniprot.aa_sequence"))
            / 1e6
        )

        # Upper bound represents protein concentration measurement, slack is introduced via lower bound
        prot_bound = protdl.upper_bound
        protdl.bounds = (prot_bound * (1 - slack_value), prot_bound)
        relax_extra = prot_bound * slack_value * protmw
        total_slack_mg_prot_per_gDW += relax_extra

        # Add relaxation reaction
        rid = protdl.id.replace("PROTDL", "RELAX")
        relaxdl = ProteinDilution(
            id=rid,
            name=protdl.name,
            subsystem=f"Pseudoreactions, Relaxation concentrations",
            lower_bound=0,
            upper_bound=DEFAULT_CONCENTRATION_BOUND,
        )
        # Store MW for relaxation bound
        relax_to_protmw_dict[rid] = protmw
        if hemoglobin_budget in protdl.metabolites:
            coeff = protdl.get_coefficient(hemoglobin_budget)  # mg/nmol
            hb_slack_mg_prot_per_gDW += relax_extra
        else:
            coeff = protdl.get_coefficient(proteome_budget)  # mg/nmol

        relaxdl.add_metabolites(
            {
                relaxation_budget: coeff,  # nmol/mg conversion factor
                protein: 1,
            },
            combine=False,
        )
        pcmodel.add_reactions([relaxdl])

    relaxation_demand.bounds = (0, total_slack_mg_prot_per_gDW)
    for relaxdl, protmw in relax_to_protmw_dict.items():
        relaxdl = pcmodel.reactions.get_by_id(relaxdl)
        relaxdl.bounds = (0, (total_slack_mg_prot_per_gDW / protmw) * slack_value)

    if verbose:
        print(
            f"Relaxation budget added to {pcmodel}, extra {total_slack_mg_prot_per_gDW:.4f} mg/gDW ({hb_slack_mg_prot_per_gDW:.4f} mg HB/gDW) from {100 * slack_value:.4f}% slack"
        )


def update_slack_value(pcmodel, slack_value, verbose):
    relax_to_protmw_dict = {}
    total_slack_mg_prot_per_gDW = 0
    hb_slack_mg_prot_per_gDW = 0
    for protdl in pcmodel.reactions.query(lambda x: x.id.startswith("PROTDL")):
        protmw = (
            calculate_protein_molar_mass(protdl.annotation.get("uniprot.aa_sequence"))
            / 1e6
        )
        # Store MW for relaxation bound
        relax_to_protmw_dict[protdl.id.replace("PROTDL", "RELAX")] = protmw

        # Upper bound represents protein concentration measurement, slack is introduced via lower bound
        prot_bound = protdl.upper_bound
        protdl.bounds = (prot_bound * (1 - slack_value), prot_bound)
        relax_extra = prot_bound * slack_value * protmw
        total_slack_mg_prot_per_gDW += relax_extra
        hemoglobin_budget = pcmodel.metabolites.get_by_id("hemoglobin_budget")
        if hemoglobin_budget in protdl.metabolites:
            hb_slack_mg_prot_per_gDW += relax_extra

    relaxation_demand = pcmodel.reactions.get_by_id(f"PBDL_relaxation_budget")
    relaxation_demand.bounds = (0, total_slack_mg_prot_per_gDW)
    for relaxdl, protmw in relax_to_protmw_dict.items():
        relaxdl = pcmodel.reactions.get_by_id(relaxdl)
        relaxdl.bounds = (
            0,
            (total_slack_mg_prot_per_gDW / protmw) * slack_value,
        )

    if verbose:
        print(
            f"Relaxation budget updated for {pcmodel}, extra {total_slack_mg_prot_per_gDW:.4f} mg/gDW ({hb_slack_mg_prot_per_gDW:.4f} mg HB/gDW) from {100 * slack_value:.4f}% slack"
        )
