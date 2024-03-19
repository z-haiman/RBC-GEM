"""Functions to extract relevant annotation information from UniProt-KB.

Notes
-----
Code based on samples provided from https://www.uniprot.org/help/id_mapping
Page last modified date: Thu Oct. 13 2022

"""

import csv
import json
import logging
import re
import time
import zlib
from urllib.parse import parse_qs, urlencode, urlparse
from warnings import warn
from xml.etree import ElementTree

import pandas as pd
import requests
from requests.adapters import HTTPAdapter, Retry

from rbc_gem_utils.util import build_string, split_string


LOGGER = logging.getLogger(__name__)

POLLING_INTERVAL = 3
UNIPROT_API_URL = "https://rest.uniprot.org"
UNIPROT_RELEASE_RE = re.compile(r"UniProt Release (?P<release>\d+_\d+)")
UNIPROT_ID_RE = re.compile(
    r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"
)
UNIPROT_ISOFORM_ID_RE = re.compile(
    "|".join([s + "-[0-9]{1,2}" for s in UNIPROT_ID_RE.pattern.split("|")])
)
UNIPROT_PATH = "/UniProt"
# Fields that are completely irrelevant to human RBCs are commented out below, but they
# are left here nonetheless to indicate that they were intentionally excluded.
# Extracted from https://www.uniprot.org/help/return_fields,
# Extracted from https://www.uniprot.org/help/return_fields_databases
# Note some corrections were made due to what looked like typos in table
# TODO Add the BeautifulSoup methods for table extraction
UNIPROT_VERSION_EXPECTED = "2024_01"
UNIPROT_DB_TAG = "UniProt"
UNIPROT_QUERY_LABEL_MIRIAM = {
    # query field: label, https://identifiers.org/
    ## Names & Taxonomy
    "accession": ("Entry", "uniprot"),
    "id": ("Entry Name", ""),
    "gene_names": ("Gene Names", ""),
    "gene_primary": ("Gene Names (primary)", "hgnc.symbol"),
    "gene_synonym": ("Gene Names (synonym)", ""),
    "gene_oln": ("Gene Names (ordered locus)", ""),
    "gene_orf": ("Gene Names (ORF)", ""),
    "organism_name": ("Organism", ""),
    "organism_id": ("Organism (ID)", "taxonomy"),
    "protein_name": ("Protein names", ""),
    "xref_proteomes": ("Proteomes", ""),
    "lineage": ("Taxonomic lineage", ""),
    "lineage_ids": ("Taxonomic lineage (Ids)", ""),
    "virus_hosts": ("Virus hosts", ""),
    ## Sequences
    "cc_alternative_products": ("Alternative products (isoforms)", "uniprot.isoform"),
    "ft_var_seq": ("Alternative sequence", ""),
    "error_gmodel_pred": ("Erroneous gene model prediction", ""),
    "fragment": ("Fragment", ""),
    "organelle": ("Gene encoded by", ""),
    "length": ("Length", ""),
    "mass": ("Mass", ""),
    "cc_mass_spectrometry": ("Mass spectrometry", ""),
    "ft_variant": ("Natural variant", ""),
    "ft_non_cons": ("Non-adjacent residues", ""),
    "ft_non_std": ("Non-standard residue", ""),
    "ft_non_ter": ("Non-terminal residue", ""),
    "cc_polymorphism": ("Polymorphism", ""),
    "cc_rna_editing": ("RNA Editing", ""),
    "sequence": ("Sequence", ""),
    "cc_sequence_caution": ("Sequence caution", ""),
    "ft_conflict": ("Sequence conflict", ""),
    "ft_unsure": ("Sequence uncertainty", ""),
    "sequence_version": ("Sequence version", ""),
    ## Function
    "absorption": ("Absorption", ""),
    "ft_act_site": ("Active site", ""),
    "cc_activity_regulation": ("Activity regulation", ""),
    "ft_binding": ("Binding site", ""),
    "cc_catalytic_activity": ("Catalytic activity", ""),
    "cc_cofactor": ("Cofactor", ""),
    "ft_dna_bind": ("DNA binding", ""),
    "ec": ("EC number", "ec-code"),
    "cc_function": ("Function [CC]", ""),
    "kinetics": ("Kinetics", ""),
    "cc_pathway": ("Pathway", ""),
    "ph_dependence": ("pH dependence", ""),
    "redox_potential": ("Redox potential", ""),
    "rhea": ("Rhea ID", "rhea"),
    "ft_site": ("Site", ""),
    "temp_dependence": ("Temperature dependence", ""),
    ## Miscellaneous
    "annotation_score": ("Annotation", ""),
    "cc_caution": ("Caution", ""),
    "comment_count": ("Comments", ""),
    "feature_count": ("Features", ""),
    "keywordid": ("Keyword ID", ""),
    "keyword": ("Keywords", ""),
    "cc_miscellaneous": ("Miscellaneous [CC]", ""),
    "protein_existence": ("Protein existence", ""),
    "reviewed": ("Reviewed", ""),
    "tools": ("Tools", ""),
    "uniparc_id": ("UniParc", "uniparc"),
    ## Interaction
    "cc_interaction": ("Interacts with", ""),
    "cc_subunit": ("Subunit structure", ""),
    ## Expression
    "cc_developmental_stage": ("Developmental stage", ""),
    "cc_induction": ("Induction", ""),
    "cc_tissue_specificity": ("Tissue specificity", ""),
    ## Gene Ontology (GO)
    "go_p": ("Gene Ontology (biological process)", ""),
    "go_c": ("Gene Ontology (cellular component)", ""),
    "go": ("Gene Ontology (GO)", ""),
    "go_f": ("Gene Ontology (molecular function)", ""),
    "go_id": ("Gene Ontology IDs", "go"),
    ## Pathology & Biotech
    "cc_allergen": ("Allergenic Properties", ""),
    "cc_biotechnology": ("Biotechnological use", ""),
    "cc_disruption_phenotype": ("Disruption phenotype", ""),
    "cc_disease": ("Involvement in disease", ""),
    "ft_mutagen": ("Mutagenesis", ""),
    "cc_pharmaceutical": ("Pharmaceutical use", ""),
    "cc_toxic_dose": ("Toxic dose", ""),
    ## Subcellular location
    "ft_intramem": ("Intramembrane", ""),
    "cc_subcellular_location": ("Subcellular location [CC]", ""),
    "ft_topo_dom": ("Topological domain", ""),
    "ft_transmem": ("Transmembrane", ""),
    ## PTM / Processsing
    "ft_chain": ("Chain", "uniprot.chain"),
    "ft_crosslnk": ("Cross-link", ""),
    "ft_disulfid": ("Disulfide bond", ""),
    "ft_carbohyd": ("Glycosylation", ""),
    "ft_init_met": ("Initiator methionine", ""),
    "ft_lipid": ("Lipidation", ""),
    "ft_mod_res": ("Modified residue", ""),
    "ft_peptide": ("Peptide", ""),
    "cc_ptm": ("Post-translational modification", ""),
    "ft_propep": ("Propeptide", ""),
    "ft_signal": ("Signal peptide", ""),
    "ft_transit": ("Transit peptide", ""),
    ## Structure
    "structure_3d": ("3D", ""),
    "ft_strand": ("Beta strand", ""),
    "ft_helix": ("Helix", ""),
    "ft_turn": ("Turn", ""),
    ## Publications
    "lit_pubmed_id": ("PubMed ID", "pubmed"),
    ## Date of
    "date_created": ("Date of creation", ""),
    "date_modified": ("Date of last modification", ""),
    "date_sequence_modified": ("Date of last sequence modification", ""),
    "version": ("Entry version", ""),
    ## Family & Domains
    "ft_coiled": ("Coiled coil", ""),
    "ft_compbias": ("Compositional bias", ""),
    "cc_domain": ("Domain [CC]", ""),
    "ft_domain": ("Domain [FT]", ""),
    "ft_motif": ("Motif", ""),
    "protein_families": ("Protein families", ""),
    "ft_region": ("Region", ""),
    "ft_repeat": ("Repeat", ""),
    "ft_zn_fing": ("Zinc finger", ""),
    ## Cross-references
    ## Sequences databases
    "xref_ccds": ("CCDS", "ccds"),
    "xref_embl": ("EMBL", "ena.embl"),
    "xref_refseq": ("RefSeq", "refseq"),
    "xref_pir": ("PIR", ""),
    ## 3D structure databases
    "xref_bmrb": ("BMRB", "bmrb"),
    "xref_pcddb": ("PCDDB", ""),
    "xref_pdb": ("PDB", "pdb"),
    "xref_pdbsum": ("PDBsum", ""),
    "xref_sasbdb": ("SASBDB", "sasbdb"),
    "xref_smr": ("SMR", "smr"),
    ## Protein-protein interaction databases
    "xref_biogrid": ("BioGRID", "biogrid"),
    "xref_corum": ("CORUM", "corum"),
    "xref_complexportal": ("ComplexPortal", "complexportal"),
    "xref_dip": ("DIP", "dip"),
    "xref_elm": ("ELM", "elm"),
    "xref_intact": ("IntAct", "intact"),
    "xref_mint": ("MINT", "mint"),
    "xref_string": ("STRING", "string"),
    ## Chemistry databases
    "xref_bindingdb": ("BindingDB", "bindingdb"),
    "xref_chembl": ("ChEMBL", "chembl.target"),
    "xref_drugbank": ("DrugBank", "drugbank"),
    "xref_drugcentral": ("DrugCentral", "drugcentral"),
    "xref_guidetopharmacology": ("GuidetoPHARMACOLOGY", "iuphar.receptor"),
    "xref_swisslipids": ("SwissLipids", "slm"),
    ## Protein family/group databases
    "xref_allergome": ("Allergome", "allergome"),
    "xref_cazy": ("CAZy", "cazy"),
    "xref_clae": ("CLAE", ""),
    "xref_esther": ("ESTHER", ""),
    "xref_ideal": ("IDEAL", "ideal"),
    "xref_imgt_gene-db": ("IMGT_GENE-DB", ""),
    "xref_moondb": ("MoonDB", ""),
    "xref_moonprot": ("MoonProt", ""),
    "xref_merops": ("MEROPS", "merops"),
    "xref_peroxibase": ("PeroxiBase", "peroxibase"),
    "xref_rebase": ("REBASE", "rebase"),
    "xref_tcdb": ("TCDB", "tcdb"),
    "xref_unilectin": ("UniLectin", ""),
    ## PTM databases
    "xref_carbonyldb": ("CarbonylDB", ""),
    "xref_depod": ("DEPOD", "depod"),
    "xref_glyconnect": ("GlyConnect", ""),
    "xref_glygen": ("GlyGen", ""),
    "xref_iptmnet": ("iPTMnet", ""),
    "xref_metosite": ("MetOSite", ""),
    "xref_phosphositeplus": ("PhosphoSitePlus", "phosphosite.protein"),
    "xref_swisspalm": ("SwissPalm", ""),
    "xref_unicarbkb": ("UniCarbKB", ""),
    ## Genetic variation/Polymorphism and mutation databases
    "xref_biomuta": ("BioMuta", ""),
    "xref_dbsnp": ("dbSNP", "dbsnp"),
    "xref_dmdm": ("DMDM", ""),
    ## 2D gel
    "xref_compluyeast-2dpage": ("COMPLUYEAST-2DPAGE", "compluyeast"),
    "xref_dosac-cobs-2dpage": ("DOSAC-COBS-2DPAGE", ""),
    "xref_ogp": ("OGP", ""),
    "xref_reproduction-2dpage": ("REPRODUCTION-2DPAGE", ""),
    "xref_swiss-2dpage": ("SWISS-2DPAGE", ""),
    "xref_ucd-2dpage": ("UCD-2DPAGE", ""),
    "xref_world-2dpage": ("World-2DPAGE", ""),
    ## Proteomic databases
    "xref_cptac": ("CPTAC", ""),
    "xref_epd": ("EPD", "epd"),
    "xref_jpost": ("jPOST", ""),
    "xref_massive": ("MassIVE", "massive"),
    "xref_maxqb": ("MaxQB", "maxqb"),
    "xref_pride": ("PRIDE", "pride"),
    "xref_paxdb": ("PaxDb", "paxdb.protein"),
    "xref_peptideatlas": ("PeptideAtlas", "peptideatlas"),
    "xref_proteomicsdb": ("ProteomicsDB", "proteomicsdb.protein"),
    "xref_promex": ("ProMEX", ""),
    "xref_topdownproteomics": ("TopDownProteomics", ""),
    ## Protocols and materials databases
    "xref_abcd": ("ABCD", ""),
    "xref_antibodypedia": ("Antibodypedia", ""),
    "xref_cptc": ("CPTC", ""),
    "xref_dnasu": ("DNASU", ""),
    ## Genome annotation databases
    "xref_ensembl": ("Ensembl", "ensembl"),
    "xref_ensemblbacteria": ("EnsemblBacteria", "ensembl.bacteria"),
    "xref_ensemblfungi": ("EnsemblFungi", "ensembl.fungi"),
    "xref_ensemblmetazoa": ("EnsemblMetazoa", "ensembl.metazoa"),
    "xref_ensemblplants": ("EnsemblPlants", "ensembl.plant"),
    "xref_ensemblprotists": ("EnsemblProtists", "ensembl.protist"),
    # "xref_genedb": ("GeneDB", "genedb"), # Broken?
    "xref_geneid": ("GeneID", "ncbigene"),
    "xref_gramene": ("Gramene", "gramene.genes"),
    "xref_kegg": ("KEGG", "kegg.genes"),
    "xref_patric": ("PATRIC", ""),
    "xref_ucsc": ("UCSC", ""),
    "xref_vectorbase": ("VectorBase", "vectorbase"),
    "xref_wbparasite": ("WBParaSite", "wbparasite"),
    ## Organism-specific
    "xref_arachnoserver": ("ArachnoServer", "arachnoserver"),
    "xref_araport": ("Araport", ""),
    "xref_cgd": ("CGD", "cgd"),
    "xref_conoserver": ("ConoServer", "conoserver"),
    "xref_ctd": ("CTD", "ctd"),
    "xref_dictybase": ("dictyBase", "dictybase.gene"),
    "xref_disgenet": ("DisGeNET", ""),
    "xref_echobase": ("EchoBASE", "echobase"),
    "xref_euhcvdb": ("euHCVdb", ""),
    "xref_flybase": ("FlyBase", ""),
    "xref_genecards": ("GeneCards", "genecards"),
    "xref_genereviews": ("GeneReviews", ""),
    "xref_hgnc": ("HGNC", "hgnc"),
    "xref_hpa": ("HPA", "hpa"),
    "xref_legiolist": ("LegioList", ""),
    "xref_leproma": ("Leproma", ""),
    "xref_maizegdb": ("MaizeGDB", ""),
    "xref_malacards": ("MalaCards", ""),
    "xref_mgi": ("MGI", "MGI"),
    "xref_mim": ("MIM", "mim"),
    "xref_nextprot": ("neXtProt", "nextprot"),
    "xref_niagads": ("NIAGADS", ""),
    "xref_opentargets": ("OpenTargets", ""),
    "xref_orphanet": ("Orphanet", "orphanet"),
    "xref_pharmgkb": ("PharmGKB", "pharmgkb.gene"),
    "xref_phi-base": ("PHI-base", ""),
    "xref_pombase": ("PomBase", "pombase"),
    "xref_pseudocap": ("PseudoCAP", ""),
    "xref_rgd": ("RGD", "rgd"),
    "xref_sgd": ("SGD", "sgd"),
    "xref_tair": ("TAIR", "tair.gene"),
    "xref_tuberculist": ("TubercuList", "myco.tuber"),
    "xref_veupathdb": ("VEuPathDB", ""),
    "xref_vgnc": ("VGNC", "vgnc"),
    "xref_wormbase": ("WormBase", "wb"),
    "xref_xenbase": ("Xenbase", "xenbase"),
    "xref_zfin": ("ZFIN", "zfin"),
    ## Phylogenomic databases
    "xref_eggnog": ("eggNOG", "eggnog"),
    "xref_genetree": ("GeneTree", "genetree"),
    "xref_hogenom": ("HOGENOM", "hogenom"),
    "xref_inparanoid": ("InParanoid", ""),
    "xref_ko": ("KO", ""),
    "xref_oma": ("OMA", "oma.grp"),
    "xref_orthodb": ("OrthoDB", "orthodb"),
    "xref_phylomedb": ("PhylomeDB", "phylomedb"),  # Not a direct mapping
    "xref_treefam": ("TreeFam", "treefam"),
    ## Enzyme and pathway databases
    "xref_biocyc": ("BioCyc", "biocyc"),
    "xref_brenda": ("BRENDA", "brenda"),
    "xref_pathwaycommons": ("PathwayCommons", "pathwaycommons"),
    ## "xref_plantreactome": ("PlantReactome", ""),
    "xref_reactome": ("Reactome", "reactome"),
    "xref_sabio-rk": ("SABIO-RK", ""),
    "xref_signalink": ("SignaLink", ""),
    "xref_signor": ("SIGNOR", "signor"),
    "xref_unipathway": ("UniPathway", ""),
    ## Miscellaneous databases
    "xref_biogrid-orcs": ("BioGRID-ORCS", ""),
    "xref_chitars": ("ChiTaRS", ""),
    "xref_evolutionarytrace": ("EvolutionaryTrace", ""),
    "xref_genewiki": ("GeneWiki", "genewiki"),
    "xref_genomernai": ("GenomeRNAi", ""),
    "xref_pharos": ("Pharos", ""),
    "xref_pro": ("PRO", ""),
    "xref_rnact": ("RNAct", ""),
    ## Gene expression databases
    "xref_bgee": ("Bgee", "bgee.gene"),
    "xref_collectf": ("CollecTF", ""),
    "xref_expressionatlas": ("ExpressionAtlas", "gxa"),
    "xref_genevisible": ("Genevisible", ""),
    "xref_cleanex": ("CleanEx", ""),
    ## Family and domain databases
    "xref_cdd": ("CDD", "cdd"),
    "xref_disprot": ("DisProt", "disprot"),
    "xref_gene3d": ("Gene3D", ""),
    "xref_hamap": ("HAMAP", "hamap"),
    "xref_interpro": ("InterPro", "interpro"),
    "xref_ncbifam": ("NCBIfam", ""),
    "xref_panther": ("PANTHER", "panther.family"),
    "xref_pfam": ("Pfam", "pfam"),
    "xref_pirsf": ("PIRSF", "pirsf"),
    "xref_prints": ("PRINTS", "prints"),
    # "xref_prodom": ("ProDom", "prodom"),
    "xref_prosite": ("PROSITE", "prosite"),
    "xref_sfld": ("SFLD", ""),
    "xref_smart": ("SMART", "smart"),
    "xref_supfam": ("SUPFAM", "supfam"),
}


def get_isoform_value_from_entry_UniProt(entry, isoform_id=None, fillna=False):
    """Extract the values associated with a given UniProt isoform"""
    if str(entry).lower() == "nan":
        if not fillna:
            return entry
        entry = entry.fillna("")

    values = split_string(entry.rstrip(";"))
    isoform_values = [x for x in values if UNIPROT_ISOFORM_ID_RE.search(x)]
    if isoform_id:
        values = [
            x
            for x in isoform_values
            if UNIPROT_ISOFORM_ID_RE.search(x).group() == isoform_id
        ]

    if isoform_values:
        values = [x.split(" ")[0] for x in values]
    entry = build_string([s for value in values for s in value.split(" ")])
    return entry


def get_label_miriam_mapping_UniProt(fields):
    """Return UniProt labels mapped to MIRIAM.

    Parameters
    ----------
    fields : list
        Query fields to get mapping for.

    Returns
    -------
    dict
        Contains UniProt lables mapped to MIRIAM values.
        Fields without MIRIAM are mapped to empty strings.

    """
    return dict(UNIPROT_QUERY_LABEL_MIRIAM[k] for k in fields)


def get_query_fields_UniProt(miriam_only=True):
    """Return UniProt query fields.

    Parameters
    ----------
    miriam_only : bool
        Whether to return only the query fields that have MIRIAM mappings associated.
        Default is ``True``

    Returns
    -------
    list
        Contains query fields for UniProt.
    """
    if miriam_only:
        return [k for k, v in UNIPROT_QUERY_LABEL_MIRIAM.items() if v[1]]
    else:
        return list(UNIPROT_QUERY_LABEL_MIRIAM)


def query_UniProt(
    query_ids,
    query_parameters,
    from_db="UniProtKB",
    to_db="UniProtKB",
    return_failed=True,
):
    """Query the UniProt database and download the results. Requires internet connection.

    Default values are used based on the RBC-GEM repository format.

    Parameters
    ----------
    query_ids : iterable
        The IDs to search in the query.
    from_db : str
        The database where the original IDs are from.
        See https://www.uniprot.org/help/return_fields_databases more
    query_parameters : dict
        Contains parameters for the query.
        See https://www.uniprot.org/help/api for more information.
    return_failed : bool
        Whether the failed IDs should be returned. Default is ``True``.

    Returns
    -------
    DataFrame
        Contains the results from the query.

    """
    try:
        assert len(set(query_ids)) == len(
            query_ids
        ), "Duplicate IDs in list to query, will correct before query"
    except AssertionError as e:
        warn(str(e))
        query_ids = set(query_ids)

    # Create session
    session = create_session()
    # Map to UniProt, IDs that fail typically are secondary accessions, considered obsolete, or are unreviewed
    job_id = submit_id_mapping(
        from_db=from_db, to_db=to_db, ids=query_ids  # To UniProtKB IDs
    )
    if check_id_mapping_results_ready(session, job_id):
        link = get_id_mapping_results_link(session, job_id)
        link = add_query_parameters(link, **query_parameters)
        results = get_id_mapping_results_search(session, link)

        df_results = get_data_frame_from_tsv_results(results)
        failed_ids = sorted(set(query_ids).difference(df_results["From"]))
        if failed_ids:
            LOGGER.warning(f"Number of failed query IDs : {len(failed_ids)}")
            request = session.get(f"{UNIPROT_API_URL}/idmapping/status/{job_id}")
            check_response(request)
            result = request.json()
            unmapped_ids = set(failed_ids).difference(result.get("failedIds", {}))
            failed_ids = set(result.get("failedIds", {}))
            if len(failed_ids) != 0:
                LOGGER.warning(f"Number of failed IDs : {len(failed_ids)}")
            obsolete_count = result.get("obsoleteCount", 0)
            if obsolete_count != 0:
                LOGGER.warning(f"Number of obsolete IDs : {obsolete_count}")

            uniparc = uniparc = {
                suggestion["from"]: suggestion["to"]
                for suggestion in result.get("suggestedIds", {})
                if suggestion
            }
        else:
            unmapped_ids, failed_ids, uniparc = set(), set(), dict()

        if return_failed:
            return df_results, uniparc, failed_ids, unmapped_ids

        return df_results


def get_version_UniProt():
    """Return the current version of UniProt."""
    response = requests.get(
        "https://ftp.uniprot.org/pub/databases/uniprot/relnotes.txt"
    )
    lines = response.text.split("\n")
    # First line should have the version information
    match = UNIPROT_RELEASE_RE.match(lines[0])
    if match is None:
        warn("Could not retrieve version")
        return None

    return match.group("release")


def get_annotation_to_from_db_UniProt(miriam_only=True):
    """Return possible databases for UniProtID mapping.

    Parameters
    ----------
    miriam_only : bool
        Return only the databases fields that have MIRIAM mappings associated.
        Default is ``True``.

    Returns
    -------
    dict
        Contains UniProt lables mapped to MIRIAM values.
        Fields without MIRIAM are mapped to empty strings

    """
    # Request for all cross-referenced databases, all pulled in at once (good for <500)
    response = requests.get(
        f"{UNIPROT_API_URL}/database/stream?compressed=False&format=json&query=%28*%29"
    )
    response.raise_for_status()
    results_json = json.loads(response.text)

    # Ensure UniProt understands that it can map to itself :)
    annotation_to_fromdb = {
        "uniprot": "UniProtKB",
        "uniprot.isoform": "UniProtKB-Swiss-Prot",
    }
    label_miriam_mapping = get_label_miriam_mapping_UniProt(
        get_query_fields_UniProt(miriam_only)
    )
    for result in results_json["results"]:
        abbrev = result["abbrev"]
        if label_miriam_mapping.get(abbrev):
            annotation_to_fromdb[label_miriam_mapping.get(abbrev)] = abbrev

    return annotation_to_fromdb


def parse_chains_UniProt(df_chains):
    chain_re = re.compile(r"id=\W(?P<chain_id>PRO_\d+)")
    for idx, (uniprot_id, uniprot_chain_value) in df_chains.loc[
        :, ["uniprot", "uniprot.chain"]
    ].iterrows():
        df_chains.loc[idx, "uniprot.chain"] = build_string(
            [
                chain_re.search(chain).group("chain_id")
                for chain in uniprot_chain_value.split("; /")
                if chain_re.search(chain)
            ]
        )
    return df_chains


def parse_isoforms_UniProt(
    df_isoforms, add_canonical=False, fill_missing_isoform=False
):
    """TODO DOCSTRING"""
    total_count = (
        df_isoforms["uniprot.isoform"]
        .replace("ALTERNATIVE PRODUCTS:  ", float("nan"))
        .count()
    )
    LOGGER.info("Number of entries with defined isoforms: {total_count}", total_count)
    if add_canonical:
        df_isoforms["uniprot.canonical"] = None
    for idx, (uniprot_id, uniprot_isoform_value) in df_isoforms.loc[
        :, ["uniprot", "uniprot.isoform"]
    ].iterrows():
        num_isoforms = re.search(r"Named isoforms=(\d+);", uniprot_isoform_value)
        isoforms = []
        canonical = uniprot_id
        if num_isoforms:
            expected_num = int(num_isoforms.group(1))
            for isoform in re.finditer(
                f"IsoId=(?P<isoid>{UNIPROT_ISOFORM_ID_RE.pattern})",
                uniprot_isoform_value,
            ):
                isoform_id = isoform.group("isoid")
                if not isoform_id.startswith(uniprot_id):
                    LOGGER.info(
                        "Isoform is an alternate entry, decreasing expected number.",
                        len(isoforms),
                        expected_num,
                    )
                    expected_num -= 1
                    continue
                s = re.search(isoform_id, uniprot_isoform_value).start()
                match = re.search("Sequence=(.+?(?=;))", uniprot_isoform_value[s:])
                if match.group(1) == "Displayed":
                    canonical = isoform_id
                isoforms += [isoform_id]
                if len(isoforms) != expected_num:
                    LOGGER.info(
                        "Number of parsed isoforms  (%s) does not match expected number (%s).",
                        len(isoforms),
                        expected_num,
                    )
        else:
            # Add the `-1` suffix for the isoform, but do not use it for the canonical form
            if fill_missing_isoform:
                isoforms += [f"{uniprot_id}-1"]
        df_isoforms.loc[idx, "uniprot.isoform"] = build_string(isoforms)
        if add_canonical:
            df_isoforms.loc[idx, "uniprot.canonical"] = canonical

    return df_isoforms


# Code taken from https://www.uniprot.org/help/id_mapping and adapted.
def create_session():
    retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
    session = requests.Session()
    session.mount("https://", HTTPAdapter(max_retries=retries))
    return session


def check_response(response):
    try:
        response.raise_for_status()
    except requests.HTTPError:
        print(response.text)
        raise


def submit_id_mapping(from_db, to_db, ids):
    request = requests.post(
        f"{UNIPROT_API_URL}/idmapping/run",
        data={"from": from_db, "to": to_db, "ids": ",".join(ids)},
    )
    check_response(request)
    return request.json()["jobId"]


def check_id_mapping_results_ready(session, job_id):
    while True:
        request = session.get(f"{UNIPROT_API_URL}/idmapping/status/{job_id}")
        check_response(request)
        j = request.json()
        if "jobStatus" in j:
            if j["jobStatus"] == "RUNNING":
                print(f"Retrying in {POLLING_INTERVAL}s")
                time.sleep(POLLING_INTERVAL)
            else:
                raise Exception(j["jobStatus"])
        else:
            return bool(j["results"] or j["failedIds"])


def get_id_mapping_results_link(session, job_id):
    url = f"{UNIPROT_API_URL}/idmapping/details/{job_id}"
    request = session.get(url)
    check_response(request)
    return request.json()["redirectURL"]


def get_query_id_categories_UniProt(session, job_id):
    request = session.get(f"{UNIPROT_API_URL}/idmapping/status/{job_id}")
    check_response(request)
    result = request.json()
    print(request.keys())
    try:
        obselete_count = result["obseleteCount"]
        failed_ids = result["failedIds"]
        uniparc = {
            suggestion["from"]: suggestion["to"]
            for suggestion in result["suggestedIds"]
        }
    except KeyError:
        return dict(), [], 0

    else:
        return uniparc, failed_ids, obselete_count


def get_next_link(headers):
    re_next_link = re.compile(r'<(.+)>; rel="next"')
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)


def get_batch(session, batch_response, file_format, compressed):
    batch_url = get_next_link(batch_response.headers)
    while batch_url:
        batch_response = session.get(batch_url)
        batch_response.raise_for_status()
        yield decode_results(batch_response, file_format, compressed)
        batch_url = get_next_link(batch_response.headers)


def combine_batches(all_results, batch_results, file_format):
    if file_format == "json":
        for key in ("results", "failedIds"):
            if key in batch_results and batch_results[key]:
                all_results[key] += batch_results[key]
    elif file_format == "tsv":
        return all_results + batch_results[1:]
    else:
        return all_results + batch_results
    return all_results


def decode_results(response, file_format, compressed):
    if compressed:
        decompressed = zlib.decompress(response.content, 16 + zlib.MAX_WBITS)
        if file_format == "json":
            j = json.loads(decompressed.decode("utf-8"))
            return j
        elif file_format == "tsv":
            return [line for line in decompressed.decode("utf-8").split("\n") if line]
        elif file_format == "xlsx":
            return [decompressed]
        elif file_format == "xml":
            return [decompressed.decode("utf-8")]
        else:
            return decompressed.decode("utf-8")
    elif file_format == "json":
        return response.json()
    elif file_format == "tsv":
        return [line for line in response.text.split("\n") if line]
    elif file_format == "xlsx":
        return [response.content]
    elif file_format == "xml":
        return [response.text]
    return response.text


def get_xml_namespace(element):
    m = re.match(r"\{(.*)\}", element.tag)
    return m.groups()[0] if m else ""


def merge_xml_results(xml_results):
    merged_root = ElementTree.fromstring(xml_results[0])
    for result in xml_results[1:]:
        root = ElementTree.fromstring(result)
        for child in root.findall("{http://uniprot.org/uniprot}entry"):
            merged_root.insert(-1, child)
    ElementTree.register_namespace("", get_xml_namespace(merged_root[0]))
    return ElementTree.tostring(merged_root, encoding="utf-8", xml_declaration=True)


def print_progress_batches(batch_index, size, total):
    n_fetched = min((batch_index + 1) * size, total)
    print(f"Fetched: {n_fetched} / {total}")


def get_id_mapping_results_search(session, url):
    parsed = urlparse(url)
    query = parse_qs(parsed.query)
    file_format = query["format"][0] if "format" in query else "json"
    if "size" in query:
        size = int(query["size"][0])
    else:
        size = 500
        query["size"] = size
    compressed = (
        query["compressed"][0].lower() == "true" if "compressed" in query else False
    )
    parsed = parsed._replace(query=urlencode(query, doseq=True))
    url = parsed.geturl()
    request = session.get(url)
    check_response(request)
    results = decode_results(request, file_format, compressed)
    total = int(request.headers["x-total-results"])
    print_progress_batches(0, size, total)
    for i, batch in enumerate(get_batch(session, request, file_format, compressed), 1):
        results = combine_batches(results, batch, file_format)
        print_progress_batches(i, size, total)
    return results


def get_data_frame_from_tsv_results(tsv_results):
    reader = csv.DictReader(tsv_results, delimiter="\t", quotechar='"')
    return pd.DataFrame(list(reader))


def add_query_parameters(link, **query_parameters):
    parsed = urlparse(link)
    query = parse_qs(parsed.query)

    for key, value in query_parameters.items():
        if key in query:
            if key not in query_parameters:
                continue
            print(f"Replacing '{key}' in query with provided query parameter.")
        query[key] = str(value).lower()

    parsed = parsed._replace(query=urlencode(query, doseq=True))
    return parsed.geturl()
