# History
### RBC-GEM 1.2.1
* Updates to repository README files
* Filepaths are corrected to be independent of platform 
* Pre-commit hooks and GitHub workflows updated for newer versions


### RBC-GEM 1.2.0
* RBC-GEM has been officially published!

### RBC-GEM 1.1.1
* Patch to fix issues recognizing the ISSUE_TEMPLATE folder due to capitalization sensitivity.
* Update the version in version.txt to match

### RBC-GEM 1.1.0
* Develop by @z-haiman in https://github.com/z-haiman/RBC-GEM/pull/6
* Notebooks updated to generate preprint figures
* Added reactions, genes, and metabolites not included in the 1.0.0 update (Specifically, SLC7A5 [PMID: 37976448])
* Others additions indicated by proteomic data and some metabolomic tracing, filling in genes for previously defined reactions
* Added additional drug metabolic pathways found in KeGG

### RBC-GEM 1.0.0
First major release of RBC-GEM, changes will be outlined in upcoming manuscript!

**Full Changelog**: https://github.com/z-haiman/RBC-GEM/compare/0.3.0...1.0.0

### RBC-GEM 0.3.0
* Not an official release, but changes are significant enough to warrant an increase in minor version.
* Remove reactions that are "duplicated" other than having different directionality
* Remove pseudoreactions that enable flux consistency to identify dead-ends. Leave exchanges.
* Remove distinction for transcripts, ensuring only unique genes in model.
* Change all gene identifiers to HGNC symbols.
* Chemical formulas and charges are updated for some metabolites.
* Metabolite formulas are standardized
* Stoichiometric corrections for reactions
* Lipids reactions are pooled.
* As the model has stoichiometrically altered from the iAB-RBC-283 model, the ID of the model has been officially changed to RBC-GEM.


### RBC-GEM 0.2.0
Updates to iAB-RBC-283 reconstruction for the following purposes:
1. To serve as a base/draft model for the initial expansion of iAB-RBC-283 to RBC-GEM.
2. To serve as a model for comparison to the new reconstruction.
3. To link to the HumanGEM reconstruction annotations where possible in the model.

No stoichiometric changes are made in this minor update, however BIGG identifiers are altered in order to enable better comparision to the RBC-GEM. Furthermore, associations to Recon3D and HumanGEM are made here, as subsequent changes past 0.2.0 will be in an effort to modernize to current COBRA standards before iterative expansion.

RBC-GEM 0.2.0 will remain identified as iAB-RBC-283 in order to help future comparisions.

Bordbar, A., Jamshidi, N. & Palsson, B.O. iAB-RBC-283: A proteomically derived knowledge-base of erythrocyte metabolism that can be used to simulate its physiological and patho-physiological states. BMC Syst Biol 5, 110 (2011). https://doi.org/10.1186/1752-0509-5-110


### RBC-GEM 0.1.1
Repairs made to downloaded iAB-RBC-283 model, restoring it to a state similar to publication. Includes updating missing metabolite formulas, charges, and fixing identifiers

### RBC-GEM 0.1.0
First release of the RBC model in GitHub.
* Memote setup with `iAB-RBC-283.xml` downloaded from BiGG (last updated 10/31/2019).
* Repository initialized from [standard-GEM](https://github.com/MetabolicAtlas/standard-GEM)
* Repository setup to follow format of [Yeast-GEM](https://github.com/SysBioChalmers/Yeast-GEM) and [Human-GEM](https://github.com/SysBioChalmers/Human-GEM) GitHub repositories.
