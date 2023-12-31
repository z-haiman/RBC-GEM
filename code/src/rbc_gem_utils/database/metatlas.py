"""Functions to extract relevant annotation information from the HumanGEM reconstruction.

Notes
-----
Code based on Human-GEM (1.18.0)

"""
import logging
import requests
import pathlib
import pandas as pd
from rbc_gem_utils.util import RAW_GH_URL, DATABASE_PATH, ROOT_PATH, check_if_valid
from warnings import warn

LOGGER = logging.getLogger(__name__)
HUMANGEM_VERSION_EXPECTED = "1.18.0"
HUMANGEM_PATH = "/Human-GEM"
HUMANGEM_URL = f"{RAW_GH_URL}/SysBioChalmers/Human-GEM"
HUMANGEM_MODEL_FILETYPES = {'mat', 'txt', 'xlsx', 'xml', 'yml'}

# https://github.com/SysBioChalmers/Human-GEM/blob/v1.18.0/code/annotateGEM.m#L75C2-L75C2
HUMANGEM_MIRIAM = {
    'reactions': {
        'rxns': 'metatlas',
        'rxnKEGGID': 'kegg.reaction',
        'rxnBiGGID': 'bigg.reaction',
        'rxnREACTOMEID': 'reactome',
        'rxnRecon3DID': 'vmhreaction',
        'rxnMetaNetXID': 'metanetx.reaction',
        'rxnTCDBID': 'tcdb',
        'rxnRheaID': 'rhea',
        'rxnRheaMasterID': 'rhea',
    },
    'metabolites': {
        'mets': 'metatlas',
        'metBiGGID': 'bigg.metabolite',
        'metKEGGID': 'kegg.compound',
        'metHMDBID': 'hmdb',
        'metChEBIID': 'chebi',
        'metPubChemID': 'pubchem.compound',
        'metLipidMapsID': 'lipidmaps',
        'metRecon3DID': 'vmhmetabolite',
        'metMetaNetXID': 'metanetx.chemical',
    },
    'genes': {
        'genes': 'ensembl',
        'geneENSTID': 'ensembl',
        'geneENSPID': 'ensembl',
        'geneUniProtID': 'uniprot',
        'geneSymbols': 'hgnc.symbol',
        'geneEntrezID': 'ncbigene',
    },
}

def get_version_HumanGEM(branch="main"):
    """Return the version of the Human-GEM model at the specified branch.
    
    Parameters
    ----------
    branch : str
        The branch to use. Default is ``"main"`` to return the latest version of the Human-GEM.
        
    Returns
    -------
    str :
        The retrieved version as a string.
    
    """
    # Make sure this points to main branch
    response = requests.get(f"{HUMANGEM_URL}/{branch}/version.txt")
    # Only text in file is the version
    version = response.text
    return version


def get_annotations_HumanGEM(
    annotation_type,
    database_dirpath,
    annotation_columns=None,
):
    """Return columns containing annotations from the Human-GEM annotation files.

    Parameters
    ----------
    annotation_type : {'reactions', 'metabolites', 'genes'}
        The type of annotation data to use. 
    database_dirpath : str or file-like
        Path to the directory containing annotation files, 
        which are supposed to be named as 'reactions.tsv', 'metabolites.tsv', and 'genes.tsv'
    columns : list, optional
        The columns to return. If ``None`` provided, defaults to all.
    """
    valid_types = {
        'reactions': "rxns", 
        'metabolites': "mets", 
        'genes': "genes"
    }
    check_if_valid(annotation_type, valid_types, "Must be one of the following")
    df_annotations = pd.read_table(f"{database_dirpath}/{annotation_type}.tsv",)

    # Only keep a specific set of columns after ensuring they are valid
    if annotation_columns is not None:
        check_if_valid(annotation_columns, df_annotations.columns, "Unrecognized columns")
        annotation_type = valid_types[annotation_type]
        # Use columns provided after verifying they exist.
        df_annotations= df_annotations.loc[:, annotation_columns]
    
    return df_annotations


def download_database_HumanGEM(annotation_type=None, database_dirpath=None, model_filetype=None, model_version=None):
    """Download the HumanGEM database files. Requires internet connection.

    Default values are used based on the RBC-GEM repository format.
    
    Parameters
    ----------
    annotation_type : {'reactions', 'metabolites', 'genes'}
        The type(s) of annotation data to use. 
        Default is to utilize all possibile values.
    database_dirpath : str or file-like
        Path or descriptor to the directory where files should be written.
        Default value is ``"{ROOT_PATH}{DATABASE_PATH}{HUMANGEM_PATH}"``
    model_filetype : {'mat', 'txt', 'xlsx', 'xml', 'yml'}
        The type of model file(s) to download. Default value is `xml`.
        Valid values exist in ``:const:HUMANGEM_MODEL_FILETYPES``
    model_version : str
        The version of the Human-GEM model and associated annotation values. 
        Default is the value of ``:const:HUMANGEM_VERSION_EXPECTED``


    """
    # Check inputs
    if model_version is None:
        # Use expected version if None provided.
        model_version = get_version_HumanGEM(f"v{HUMANGEM_VERSION_EXPECTED}")

    valid = {'reactions', 'metabolites', 'genes'}
    if annotation_type is None:
        annotation_type = valid
    else:
        annotation_type = check_if_valid(annotation_type, valid, "Unrecognized annotations for Human-GEM")

    if model_filetype is None:
        model_filetype = ["xml"]
    else:
        model_filetype = check_if_valid(model_filetype, HUMANGEM_MODEL_FILETYPES, "Unrecognized filetypes for Human-GEM")

    if database_dirpath is not None:
        # Ensure the path exists
        pathlib.Path(f'{database_dirpath}').mkdir(parents=False, exist_ok=True)
    else:
        database_dirpath = f"{ROOT_PATH}{DATABASE_PATH}{HUMANGEM_PATH}"


    for ann_type in annotation_type:
        # FIXME probably a better way to do this instead of erroring out.

        response = requests.get(f"{HUMANGEM_URL}/v{model_version}/model/{ann_type}.tsv")
        response.raise_for_status()

        filepath = f"{database_dirpath}/{ann_type}.tsv"
        with open(filepath, "w") as file:
            file.write(response.text)

        LOGGER.info("`%s.tsv` saved at `%s`", ann_type, database_dirpath)
    
    for ftype in model_filetype:
        filename = f"Human-GEM.{ftype}"
        # FIXME probably a better way to do this instead of erroring out.
        response = requests.get(f"{HUMANGEM_URL}/v{model_version}/model/{filename}")
        response.raise_for_status()

        # Write file
        filepath = f"{database_dirpath}/{filename}"
        # Is there a better way of checking whether binary file?
        if not response.encoding:
            with open(filepath, "wb") as file:
                file.write(response.content)
        else:
            with open(filepath, "w") as file:
                file.write(response.text)
                
        LOGGER.info("`%s` saved at `%s`", filename, database_dirpath)
