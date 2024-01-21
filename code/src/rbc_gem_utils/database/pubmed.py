"""Functions to extract metadata for PubMed articles."""
import re

import numpy as np
from Bio import Entrez

from rbc_gem_utils.util import build_string


PUBMED_ERYTHROCYTE_TERMS = [
    "red blood cell",
    "red cell",
    "erythrocyte",
    "erythro",
    "rbc",
]


def search_erythrocyte_terms_PubMed(df, text_columns, search_terms=None):
    df = df.loc[:, text_columns].copy()
    key = "RBC Terms"
    df[key] = ""
    for idx, row in df.iterrows():
        hits = set()
        for value in row.values:
            hits.update(re.findall("|".join(search_terms), value))
        df.loc[idx, key] = build_string(hits)
    return df


def fetch_results_PubMed(email, pubmed_ids):
    """Query PubMed and return results. Email address required to query pubmed.

    If query is large, use the batch querying which utilizes pagination to obtain results.
    """
    Entrez.email = email
    search = Entrez.efetch(db="pubmed", sort="relevance", retmode="xml", id=pubmed_ids)
    results = Entrez.read(search)
    search.close()
    return results


def fetch_batch_results_PubMed(email, pubmed_ids, batch_size=500, return_failed=True):
    """Query PubMed in batches of a given side and return all results."""
    all_results = []
    batch_size = 500
    for idx, batch in enumerate(np.arange(0, len(pubmed_ids), batch_size), start=1):
        query_ids = pubmed_ids[batch : batch + batch_size]
        print(
            f"Fetching results for batch {idx}  ({batch + len(query_ids)}/{len(pubmed_ids)})"
        )
        results = fetch_results_PubMed(email=email, pubmed_ids=query_ids)
        all_results += results["PubmedArticle"]

    results_pmids = [f"{article['MedlineCitation']['PMID']}" for article in all_results]
    failed_ids = set(pubmed_ids).difference(results_pmids)
    if failed_ids:
        print(f"Failed {len(failed_ids)} IDs")

    if return_failed:
        return all_results, failed_ids

    return all_results


def get_mesh_terms(mesh_heading, only_major=True, use_ids=True):
    if (
        only_major
        and mesh_heading["DescriptorName"].attributes.get("MajorTopicYN") != "Y"
    ):
        return []

    if use_ids:
        descriptor = mesh_heading["DescriptorName"].attributes.get("UI")
    else:
        descriptor = mesh_heading["DescriptorName"]

    if use_ids:
        mesh_terms = [
            f"{descriptor}/{mesh.attributes.get('UI')}"
            for mesh in mesh_heading["QualifierName"]
            if (only_major and mesh.attributes.get("MajorTopicYN") == "Y")
            or not only_major
        ]
    else:
        mesh_terms = [
            f"{descriptor}/{mesh}"
            for mesh in mesh_heading["QualifierName"]
            if (only_major and mesh.attributes.get("MajorTopicYN") == "Y")
            or not only_major
        ]
    if not mesh_terms:
        mesh_terms = [str(descriptor)]

    return mesh_terms


def get_value_PubMed(item, key, subkey=None, default_value=""):
    try:
        if subkey is None:
            value = item[key]
        else:
            value = item[key][subkey]
    except KeyError as e:
        return default_value
    except IndexError as e:
        return default_value
    else:
        if isinstance(value, str):
            return str(value)
        return str(value[0])
