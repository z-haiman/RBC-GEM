# RBC-GEM: Genome-scale metabolic model for the erythrocyte of _Homo sapiens_

[![GitHub version](https://badge.fury.io/gh/z-haiman%2Frbc-gem.svg)](https://badge.fury.io/gh/z-haiman%2Frbc-gem)
[![Zenodo](https://zenodo.org/badge/733772184.svg)](https://zenodo.org/doi/10.5281/zenodo.10836860)
[![Gitter chat](https://badges.gitter.im/z-haiman/RBC-GEM.svg)](https://gitter.im/z-haiman/RBC-GEM)
[![memote tested](https://img.shields.io/badge/memote-tested-blue.svg?style=plastic)](https://z-haiman.github.io/RBC-GEM)


#### Citation
If you use RBC-GEM in your research, please cite the following:

> Haiman ZB, Key A, D'Alessandro A, Palsson BO. RBC-GEM: A genome-scale metabolic model for systems biology of the human red blood cell. PLoS Comput Biol. 2025 Mar 12;21(3):e1012109. doi: 10.1371/journal.pcbi.1012109. PMID: 40072998; PMCID: PMC11925312.

All the releases are also archived in [Zenodo](https://doi.org/10.5281/zenodo.10836860) from which specific version can be cited if used.

#### Keywords

**Utilisation:** experimental data reconstruction; multi-omics integrative analysis; model template
**Field:** metabolic-network reconstruction
**Type of model:** reconstruction; curated
**Model source:** [iAB-RBC-283](https://doi.org/10.1186/1752-0509-5-110),   [Human-GEM](https://doi.org/10.5281/zenodo.10303455)
**Omic source:** genomics; proteomics; metabolomics; lipidomics
**Taxonomic name:** _Homo sapiens_
**Taxonomy ID:** [taxonomy:9606](https://identifiers.org/taxonomy:9606)
**Metabolic system:** general metabolism
**Tissue:** Blood
**Cell type:**  Red Blood Cell (erythrocyte)
**Condition:** generic metabolism

### Installation
#### With Python

The recommended method is to install **rbc_gem_utils** is to use ``pip`` to
install the software. It is recommended to do this inside a [virtual environment](http://docs.python-guide.org/en/latest/dev/virtualenvs/). If you have [conda](https://docs.conda.io/en/latest/), follow the instructions for [managing environments](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)

1.  Clone the [main branch](https://github.com/z-haiman/RBC-GEM/tree/main) of this repository, or [download the latest release](https://github.com/z-haiman/RBC-GEM/releases/latest).
2.  Navigate to the `code` directory containing the `pyproject.toml` file and install the package:

        cd "/my/path/RBC-GEM/code"
        pip install "."

3. Test your install:

        python -c "from rbc_gem_utils import show_versions; show_versions()"



### Contributing

Contributions are always welcome! Please read the [contributing guideline](.github/CONTRIBUTING.md) to get started.


### Contributors

Code contributors are reported automatically by GitHub under [Contributors](https://github.com/z-haiman/RBC-GEM/graphs/contributors), while other contributions come in as [Issues](https://github.com/z-haiman/RBC-GEM/issues).
