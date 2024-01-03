# RBC-GEM external/databases/Human-GEM

Contains the downloaded Human-GEM model and annotation files.

## Relevant notebooks
* [Extract annotations from Metabolic Atlas via Human-GEM](../../../../code/notebooks/annotation/ExtractAnnotations_MetAtlas.ipynb)

## File Descriptions
### HumanGEM: 1.18.0
URL: https://github.com/SysBioChalmers/Human-GEM/tree/v1.18.0

* `Human-GEM.xml`: Human-GEM model file (SBML L3V1, fbc-v2, groups-v1)
* `reactions.tsv`: Contains additional annotation information and external identifiers for Human-GEM reactions.

|# |fieldname      |annotation               |Prefixes (https://identifiers.org/)|
|--|---------------|-------------------------|-----------------------------------|
|1 |rxns           |identical to `model.rxns`|metatlas                           |
|2 |rxnKEGGID      |KEGG reaction ID         |kegg.reaction                      |
|3 |rxnBiGGID      |BiGG reaction ID         |bigg.reaction                      |
|4 |rxnEHMNID      |EHMN reaction ID         |                                   |
|5 |rxnHepatoNET1ID|HepatoNET1 reaction ID   |                                   |
|6 |rxnREACTOMEID  |REACTOME ID              |reactome                           |
|7 |rxnRecon3DID   |Recon3D reaction ID      |vmhreaction                        |
|8 |rxnMetaNetXID  |MetaNetX reaction ID     |metanetx.reaction                  |
|9 |rxnHMR2ID      |HMR2 reaction ID         |                                   |
|10|rxnRatconID    |Ratcon reaction ID       |                                   |
|11|rxnTCDBID      |TCDB ID                  |tcdb                               |
|12|spontaneous    |Spontaneous status       |                                   |
|13|rxnRheaID      |Rhea ID                  |rhea                               |
|14|rxnRheaMasterID|Master Rhea ID           |rhea                               |
|15|rxnRetired     |Retired reaction IDs     |                                   |


* `metabolites.tsv`: Contains additional annotation information and external identifiers for Human-GEM metabolites.

|# |fieldname      |annotation                             |Prefixes (https://identifiers.org/)|
|--|---------------|---------------------------------------|-----------------------------------|
|1 |mets           |identical to `model.mets`              |metatlas                           |
|2 |metsNoComp     |`model.mets` without compartment suffix|                                   |
|3 |metBiGGID      |BiGG metabolite ID                     |bigg.metabolite                    |
|4 |metKEGGID      |KEGG metabolite ID                     |kegg.compound                      |
|5 |metHMDBID      |HMDB ID                                |hmdb                               |
|6 |metChEBIID     |ChEBI ID                               |chebi                              |
|7 |metPubChemID   |PubChem ID                             |pubchem.compound                   |
|8 |metLipidMapsID |LipidMaps ID                           |lipidmaps                          |
|9 |metEHMNID      |EHMN metabolite ID                     |                                   |
|10|metHepatoNET1ID|HepatoNET1 metabolite ID               |                                   |
|11|metRecon3DID   |Recon3D metabolite ID                  |vmhmetabolite                      |
|12|metMetaNetXID  |MetaNetX metabolite ID                 |metanetx.chemical                  |
|13|metHMR2ID      |HMR2 metabolite ID                     |                                   |
|14|metRetired     |Retired metabolite IDs                 |                                   |

* `genes.tsv`: Contains additional annotation information and external identifiers for Human-GEM genes.

|# |fieldname     |  annotation          |Prefixes (https://identifiers.org/)|
|--|--------------|----------------------|-----------------------------------|
|1 |genes         |Ensembl gene ID       |ensembl                            |
|2 |geneENSTID    |Ensembl transcript ID |ensembl                            |
|3 |geneENSPID    |Ensembl protein ID    |ensembl                            |
|4 |geneUniProtID |UniProt ID            |uniprot                            |
|5 |geneSymbols   |Gene Symbol           |hgnc.symbol                        |
|6 |geneEntrezID  |NCBI Entrez ID        |ncbigene                           |
|7 |geneNames     |Gene Name             |                                   |
|8 |geneAliases   |Alias Names           |                                   |
|9 |compartments  |Subcellular location  |                                   |
|10|compDataSource|Source for compartment|                                   |
