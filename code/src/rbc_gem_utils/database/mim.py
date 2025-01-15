"""Functions to extract relevant annotation information from OMIM.

Notes
-----
Main site: https://omim.org/

"""

import re
from datetime import datetime
from pathlib import Path
from warnings import warn

import pandas as pd

from rbc_gem_utils.util import DATABASE_PATH, RAW_GH_URL, ROOT_PATH, check_if_valid


MIM_PATH = "/MIM"
MIM_FILENAMES = ["mim2gene", "genemap2", "morbidmap", "mimTitles"]
MIM_NUMBER_RE = re.compile(r"^\d{6}\Z")
MIM_DB_TAG = "MIM"


def load_data_MIM(mim_filename, mim_directory=None, print_footer_notes=False):
    """Load the MIM file, ignoring any initial lines before the header that are not relevant."""
    if mim_directory is None:
        mim_directory = f"{ROOT_PATH}{DATABASE_PATH}{MIM_PATH}"
    else:
        Path(f"{mim_directory}").mkdir(parents=False, exist_ok=True)

    check_if_valid(mim_filename, MIM_FILENAMES, "Invalid MIM filename:")
    header = {
        "mim2gene": 4,
        "genemap2": 3,
        "morbidmap": 3,
        "mimTitles": 2,
    }.get(mim_filename)
    df = pd.read_csv(
        f"{mim_directory}/{mim_filename}.txt", sep="\t", header=header, dtype=str
    )
    try:
        if print_footer_notes:
            df_footer_comments = df[df["MIM Number"].isna()]
            df_footer_comments = df_footer_comments[
                df_footer_comments.columns[0]
            ].replace("#", "")
            df_footer_comments = df_footer_comments[4:].reset_index(drop=True)
            for value in df_footer_comments.values:
                print(value)

        df = df[~df["MIM Number"].isna()]
    except KeyError:
        pass
    return df.replace("", float("nan"))


def get_last_updated_dates_MIM(mim_directory=None, mim_filenames=None, verbose=True):
    """Check the MIM files to ensure that they were all generated at the same time"""
    if mim_filenames is None:
        mim_filenames = MIM_FILENAMES
    else:
        check_if_valid(mim_filenames, MIM_FILENAMES, "Invalid MIM filename:")

    if mim_directory is None:
        mim_directory = f"{ROOT_PATH}{DATABASE_PATH}{MIM_PATH}"
    else:
        Path(f"{mim_directory}").mkdir(parents=False, exist_ok=True)

    last_updated_dict = {}
    for mim_file in mim_filenames:
        with open(f"{mim_directory}/{mim_file}.txt") as file:
            for line in file:
                if line.startswith("# Generated"):
                    break
            date = line.split(" ")[-1]
            date = date.strip()
            date = datetime.strptime(date, "%Y-%m-%d").strftime("%Y-%m-%d")
            last_updated_dict[mim_file] = date

    last_updated = set(last_updated_dict.values())
    if not len(last_updated) == 1:
        warn(f"Files were not all downloaded at the same time: {last_updated_dict}")
    elif verbose:
        print(f"Files last generated on: {list(last_updated).pop()}")
    return last_updated_dict
