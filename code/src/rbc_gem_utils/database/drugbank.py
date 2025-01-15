"""Functions to extract relevant annotation information from DrugBank

Notes
-----
Code based on DrugBank (5.1.11)

"""

import pathlib
import zipfile
from collections import defaultdict
from xml.etree import ElementTree

import pandas as pd
import requests
from bs4 import BeautifulSoup

from rbc_gem_utils.util import DATABASE_PATH, ROOT_PATH, build_string, strip_plural


DRUGBANK_VERSION_EXPECTED = "5.1.13"
DRUGBANK_URL = "https://go.drugbank.com"
DRUGBANK_DB_TAG = "DrugBank"
DRUGBANK_PATH = "/DrugBank"
DRUGBANK_NS = "{http://www.drugbank.ca}"

# Fields for the DrugBank XML schema are found [here](https://docs.drugbank.com/xml/#introduction).
DRUGBANK_GENERAL_ELEMENTS = [
    "drugbank-id",
    "name",
    "description",
    # 'simple-description',
    # 'clinical-description',
    # 'therapeutically-significant',
    "cas-number",
    "unii",
    "average-mass",
    "monoisotopic-mass",
    "state",
    "groups",
    "categories",
    "affected-organisms",
    "ahfs-codes",
    "pdb-entries",
    "fda-label",
    "msds",
    "food-interactions",
    "general-references",
    "synthesis-reference",
]
DRUGBANK_PHARMACOLOGY_ELEMENTS = [
    "indication",
    "pharmacodynamics",
    "mechanism-of-action",
    "toxicity",
    "metabolism",
    "absorption",
    "half-life",
    "protein-binding",
    "route-of-elimination",
    "volume-of-distribution",
    "clearance",
]

DRUGBANK_CLASSIFICATION_ELEMENTS = [
    "kingdom",
    "superclass",
    "direct-parent",
    "subclass",
    "substituent",
    "description",
    "alternative-parent",
]

DRUGBANK_REFERENCE_ELEMENTS = [
    "articles",
    "textbooks",
    "links",
    "attachments",
]

DRUGBANK_MIXTURES_ELEMENTS = [
    "name",
    "ingredients",
    "supplemental-ingredients",
]

DRUGBANK_SALTS_ELEMENTS = [
    "drugbank-id",
    "name",
    "unii",
    "cas-number",
    "inchikey",
    "average-mass",
    "monoisotopic-mass",
    "smiles",
    # "inchi",
    "formula",
]

DRUGBANK_PRICE_ELEMENTS = ["description", "cost", "unit"]

DRUGBANK_DOSAGE_ELEMENTS = ["form", "route", "strength"]

DRUGBANK_PATHWAY_ELEMENTS = [
    "smpdb-id",
    "name",
    "category",
    "drugs",
    "enzymes",
]
DRUGBANK_PATENT_ELEMENTS = [
    "number",
    "country",
    "approved",
    "expires",
    "pediatric-extension",
]
ATC_CODE_LEVELS = {
    1: "anatomical",
    2: "therapeutic",
    3: "pharmacological",
    4: "chemical",
    5: "substance",
}


def get_version_DrugBank():
    # Get the table of releases
    response = requests.get(f"{DRUGBANK_URL}/releases")
    response.raise_for_status()

    # Parse the table to create the soup
    soup = BeautifulSoup(response.text, "html.parser")
    header = soup.table.find_next("tr")
    keys = [x.contents[0] for x in header.children if x.contents]

    # First row has latest release
    first_row = header.find_next("tr")
    values = [x.contents[0] for x in list(first_row.children)[: len(keys)]]
    # Zip the information together and return
    release_entry = dict(zip(keys, values))

    return release_entry["Version"]


def download_database_DrugBank(username, password, database_dirpath=None, version=None):
    """Download and extract the DrugBank database for the given version.

    Parameters
    ----------
    username : str
    password : str
    database_dirpath : str or path-like, optional

    """
    if database_dirpath is None:
        # Ensure the path exists
        database_dirpath = f"{ROOT_PATH}{DATABASE_PATH}{DRUGBANK_PATH}"
    else:
        pathlib.Path(f"{database_dirpath}").mkdir(parents=False, exist_ok=True)

    if version is None:
        version = DRUGBANK_VERSION_EXPECTED
    # Download and write the zip file
    version = version.replace(".", "-")
    response = requests.get(
        f"{DRUGBANK_URL}/releases/{version}/downloads/all-full-database",
        auth=(username, password),
    )
    filename = "drugbank_all_full_database.xml"
    filepath = f"{database_dirpath}/{filename}"
    response.raise_for_status()
    with open(f"{filepath}.zip", "wb") as file:
        file.write(response.content)

    # Rename
    with zipfile.ZipFile(f"{filepath}.zip", "r") as file:
        file.getinfo("full database.xml").filename = filename
        file.extract("full database.xml", path=database_dirpath)

    return version


def strip_ns_DrugBank(value):
    """Strip the DrugBank Namespace prefix from the text."""
    return value.replace(DRUGBANK_NS, "")
