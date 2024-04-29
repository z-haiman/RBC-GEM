"""Functions to extract relevant information from Complex Portal

Notes
-----
* Main site: https://www.ebi.ac.uk/complexportal/home
* Code written based on release of 27-Mar-24 (6759 active entries)
* Code last updated: May 2024

"""

import requests
from bs4 import BeautifulSoup
from datetime import datetime
import pandas as pd

COMPLEXPORTAL_URL = "https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/"
COMPLEXPORTAL_VERSION_EXPECTED = "2024-02-14"
COMPLEXPORTAL_DB_TAG = "ComplexPortal"
COMPLEXPORTAL_PATH = "/ComplexPortal"

def get_version_ComplexPortal(taxomony_int=9606):
    """Return the current version of ComplexPortal database."""
    # Get the table of releases
    response = requests.get(f"{COMPLEXPORTAL_URL}")
    response.raise_for_status()

    # Parse the table to create the soup
    soup = BeautifulSoup(response.text, "html.parser")
    header = soup.table
    taxomony_int = str(taxomony_int)
    for tr in header.find_all("tr"):
        ref = [ref.parent.parent for ref in tr.find_all("a") if ref.text.split(".")[0] == f"{taxomony_int}"]
        if not ref:
            continue
        # Third item has date and time stamp
        date_time = ref[0].find_all("td")[2].contents[0].strip()
        break

    release = datetime.strptime(date_time, "%Y-%m-%d %H:%M").strftime("%Y-%m-%d")
    return release

def parse_complex_participants(series):
    """Parse complex participants and stoichiometry. A stoichiometric cefficient of 0 indicates unclear stoichiometry."""
    series = series.str.split("|").explode()
    series.name = "participants_stoichiometry"
    col1 = series.apply(lambda x: x.split("(")[0])
    col1.name = "participants"

    col2 = series.apply(lambda x: x.split("(")[-1].rstrip(")"))
    col2.name = "stoichiometry"

    df = pd.concat((
        series, 
        col1,
        col2, 
    ), axis=1)
    
    return df