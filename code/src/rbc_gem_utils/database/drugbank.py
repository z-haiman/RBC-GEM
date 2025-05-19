"""Functions to extract relevant annotation information from DrugBank

Notes
-----
* Main site: https://go.drugbank.com/
* Code written written or updated on DrugBank (5.1.13), released 02-Jan-25
* Code last updated: May 2025
"""

import re
import zipfile
from pathlib import Path

import requests
from bs4 import BeautifulSoup

from rbc_gem_utils.util import DATABASE_PATH, ROOT_PATH, build_string, strip_plural


DRUGBANK_RELEASE_EXPECTED = "5.1.13"
DRUGBANK_URL = "https://go.drugbank.com"
DRUGBANK_DB_TAG = "DrugBank"
DRUGBANK_PATH = Path(DRUGBANK_DB_TAG)
DRUGBANK_NS = "{http://www.drugbank.ca}"
DRUGBANK_RELEASE_RE = re.compile(r"DrugBank Release Version (\d+.\d+.\d+)")

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


def get_release_DrugBank():
    # Get the table of releases
    response = requests.get(f"{DRUGBANK_URL}/releases/latest")
    response.raise_for_status()

    # Parse the table to create the soup
    soup = BeautifulSoup(response.text, "html.parser")
    header = soup.find("title")
    release = DRUGBANK_RELEASE_RE.search(header.text).group(1)

    return release


def download_database_DrugBank(username, password, database_dirpath=None, release=None):
    """Download and extract the DrugBank database for the given release.

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
        Path(f"{database_dirpath}").mkdir(parents=False, exist_ok=True)

    if release is None:
        release = DRUGBANK_RELEASE_EXPECTED
    # Download and write the zip file
    release = release.replace(".", "-")
    response = requests.get(
        f"{DRUGBANK_URL}/releases/{release}/downloads/all-full-database",
        auth=(username, password),
    )
    response.raise_for_status()
    filename = "drugbank_all_full_database.xml"
    filepath = f"{database_dirpath}/{filename}"
    with open(f"{filepath}.zip", "wb") as file:
        file.write(response.content)

    # Rename
    with zipfile.ZipFile(f"{filepath}.zip", "r") as file:
        file.getinfo("full database.xml").filename = filename
        file.extract("full database.xml", path=database_dirpath)

    return release


def strip_ns_DrugBank(value):
    """Strip the DrugBank Namespace prefix from the text."""
    return value.replace(DRUGBANK_NS, "")
