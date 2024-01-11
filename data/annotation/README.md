# RBC-GEM data/annotation


Contains model content mapped to different databases. Files stored here are intended to be used repeatedly and updated with any model or database changes. Databases and brief descriptions of each file are found below:


## All Annotations
### Reactions
|annotation                                    |Prefixes (https://identifiers.org/)|cobrapy          |RAVEN Toolbox  |COBRA Toolbox    |
|----------------------------------------------|-----------------------------------|-----------------|---------------|-----------------|
|identical to `model.reactions`                |                                   |reactions        |rxns           |rxns             |
|MetAtlas reaction ID                          |metatlas                           |metatlas         |rxns           |rxns             |
|KEGG reaction ID                              |kegg.reaction                      |kegg.reaction    |rxnKEGGID      |rxnKEGGID        |
|BiGG reaction ID                              |bigg.reaction                      |bigg.reaction    |rxns           |rxns             |
|REACTOME reaction ID                          |reactome                           |reactome         |rxnREACTOMEID  |rxnReactomeID    |
|Recon3D reaction ID                           |vmhreaction                        |vmhreaction      |rxnRecon3DID   |rxnRecon3DID     |
|MetaNetX reaction ID                          |metanetx.reaction                  |metanetx.reaction|rxnMetaNetXID  |rxnMetaNetXID    |
|TCDB reaction ID                              |tcdb                               |tcdb             |rxnTCDBID      |rxnTCDBID        |
|Master Rhea ID                                |rhea                               |rhea             |rxnRheaMasterID|rxnRheaID        |
|spontaneous                                   |                                   |spontaneous      |spontaneous    |                 |
|Enzyme Classification code for each reaction  |ec-code                            |ec-code          |eccodes        |rxnECNumbers     |
|KEGG ID of manually drawn KEGG pathway map    |kegg.pathway                       |kegg.pathway     |               |rxnKEGGPathways  |
|BRENDA ID                                     |brenda                             |brenda           |               |rxnBRENDAID      |
|BioCyc ID                                     |biocyc                             |biocyc           |               |rxnBioCycID      |
|SABIO-RK reaction ID                          |sabiork.reaction                   |sabiork.reaction |               |rxnSABIORKID     |
|SEED reaction ID                              |seed.reaction                      |seed.reaction    |               |rxnSEEDID        |
|Systems Biology Ontology                      |sbo                                |sbo              |               |rxnSBOTerms      |
|Gibbs Free Energy of reaction                 |                                   |                 |rxnDeltaG      |                 |
|                                              |                                   |                 |               |                 |

### Metabolites

|annotation                                    |Prefixes (https://identifiers.org/)|cobrapy          |RAVEN Toolbox  |COBRA Toolbox    |
|----------------------------------------------|-----------------------------------|-----------------|---------------|-----------------|
|identical to `model.metabolites`              |                                   |metabolites      |mets           |mets             |
|`model.metabolites` without compartment suffix|                                   |                 |metsNoComp     |metsNoComp       |
|MetAtlas metabolite ID                        |metatlas                           |metatlas         |mets           |mets             |
|BiGG metabolite ID                            |bigg.metabolite                    |bigg.metabolite  |metBiGGID      |metBiGGID        |
|KEGG compound ID                              |kegg.compound                      |kegg.compound    |metKEGGID      |metKEGGID        |
|HMDB ID                                       |hmdb                               |hmdb             |metHMDBID      |metHMDBID        |
|ChEBI ID                                      |chebi                              |chebi            |metChEBIID     |metChEBIID       |
|PubChem ID                                    |pubchem.compound                   |pubchem.compound |metPubChemID   |metPubChemID     |
|LipidMaps ID                                  |lipidmaps                          |lipidmaps        |metLipidMapsID |metLIPIDMAPSID   |
|Recon3D metabolite ID                         |vmhmetabolite                      |vmhmetabolite    |metRecon3DID   |metRecon3DID     |
|MetaNetX metabolite ID                        |metanetx.chemical                  |metanetx.chemical|metMetaNetXID  |metMetaNetXID    |
|Metabolite formula in InCHI string format     |inchi                              |inchi            |inchis         |metInChIString   |
|Metabolite formula in InCHI key format        |inchikey                           |inchikey         |               |                 |
|SEED metabolite ID                            |seed.compound                      |seed.compound    |               |metSEEDID        |
|BioCyc ID                                     |biocyc                             |biocyc           |               |metBioCycID      |
|SABIO-RK metabolite ID                        |sabiork.compound                   |sabiork.compound |               |metSABIORKID     |
|SwissLipids ID                                |slm                                |slm              |               |metSLMID         |
|Systems Biology Ontology                      |sbo                                |sbo              |               |metSBOTerms      |
|Metabolite formula in SMILES format           |                                   |                 |metSmiles      |metSmiles        |
|Gibbs Free Energy of formation                |                                   |                 |metDeltaG      |                 |
|                                              |                                   |                 |               |                 |
### Genes

|annotation                                    |Prefixes (https://identifiers.org/)|cobrapy          |RAVEN Toolbox  |COBRA Toolbox    |
|----------------------------------------------|-----------------------------------|-----------------|---------------|-----------------|
|identical to `model.genes`                    |                                   |genes            |genes          |genes            |
|Ensembl gene ID                               |ensembl                            |ensembl          |geneEnsemblID  |geneEnsemblID    |
|Ensembl transcript ID                         |ensembl                            |                 |geneENSTID     |                 |
|Ensembl protein ID                            |ensembl                            |                 |geneENSPID     |                 |
|UniProt ID                                    |uniprot                            |uniprot          |geneUniProtID  |geneUniprotID    |
|HGNC Gene Symbol                              |hgnc.symbol                        |hgnc.symbol      |geneSymbols    |geneSymbols      |
|NCBI Entrez ID                                |ncbigene                           |ncbigene         |geneEntrezID   |geneEntrezID     |
|RefSeq ID                                     |refseq                             |refseq           |               |geneRefSeqID     |
|KEGG gene ID                                  |kegg.genes                         |kegg.genes       |               |geneKEGGID       |
|Human Protein Reference Database ID           |hprd                               |hprd             |               |geneHPRDID       |
|Concensus CDS gene ID                         |ccds                               |ccds             |               |geneCCDSID       |
|NCBI protein ID                               |ncbiprotein                        |ncbiprotein      |               |geneNCBIProteinID|
|Systems Biology Ontology                      |sbo                                |sbo              |               |geneSBOTerms     |
|                                              |                                   |                 |               |                 |

## Descriptions
* `reactions_MetAtlas.tsv`: Contains additional annotation information and external identifiers extracted from Human-GEM.
* `metabolites_MetAtlas.tsv`: Contains additional annotation information and external identifiers for Human-GEM metabolites.
* `genes_MetAtlas.tsv`: Contains additional annotation information and external identifiers for Human-GEM genes.


##### Notes
* Often times, a database contains cross-references to other databases, making it possible to obtain additional annotation information for those databases. Data contained in files here may also have additional cross-reference information. 
* Consequently, there are occassions where a cross-reference within one database does not match the cross-reference in another database. This can occur for various reasons, such as differences in database update schedules, whether any databases have undergone backwards incompatible changes, etc.
* All information is coalesced into a single main file using a [notebook](../../code/notebooks/annotation/ReconcileAnnotations_RBC-GEM.ipynb).
* Annotations extracted directly from their primary database should always be prioritized over those extracted as a secondary cross-reference.

##### References
* Fields for COBRApy are dictionary keys for the `Object.annotation` attribute
* Fields for the [RAVEN Toolbox (2.8.6)](https://github.com/SysBioChalmers/RAVEN/blob/v2.8.6/readme/RAVEN_structure_fields.csv)
* Fields for the [COBRA Toolbox (3.33)](https://github.com/opencobra/cobratoolbox/blob/v3.33/docs/source/notes/COBRAModelFields.md)
