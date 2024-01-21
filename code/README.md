# RBC-GEM Code

This directory contains functions and scripts that facilitate work with the RBC-GEM model and other content in the repository.

This directory is organized as follows:

```
code
├── notebooks   # For Python 3
└── src         # For Python 3
```

##### Note

There are several tools that enable working with genome scale metabolic reconstructions (e.g., [cobrapy](https://github.com/opencobra/cobrapy), [COBRA Toolbox](https://github.com/opencobra/cobratoolbox/), [RAVEN Toolbox](https://github.com/SysBioChalmers/RAVEN)).

Therefore, to help maintain code clarity across multiple programming languages:

* Python code should be written using `snake_case` nomenclature and follow the [PEP 8 – Style Guide for Python Code](https://peps.python.org/pep-0008/).
* MATLAB code should be written using `mixedCase` nomenclature and follow the [MATLAB Style Guidelines](https://www.mathworks.com/matlabcentral/fileexchange/46056-matlab-style-guidelines-2-0).

## Installation
### With Python

The recommended method is to install **rbc_gem_utils** is to use ``pip`` to
install the software. It is recommended to do this inside a [virtual environment](http://docs.python-guide.org/en/latest/dev/virtualenvs/). If you have [conda](https://docs.conda.io/en/latest/), follow the instructions for [managing environments](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)

1.  Clone the [main branch](https://github.com/z-haiman/RBC-GEM/tree/main) of this repository, or [download the latest release](https://github.com/z-haiman/RBC-GEM/releases/latest).
2.  Navigate to the `code` directory containing the `pyproject.toml` file and install the package:

        cd "/my/path/RBC-GEM/code"
        pip install "."

3. Test your install:

       python -c "from rbc_gem_utils import show_versions; show_versions()"

## Description of directory contents

In an effort to keep this repository compatibile with

Brief description of directory contents are provided below.
### notebooks (Python)

**TODO**

Contains iPython (`.ipynb`) notebooks

### src (Python)

Contains the source code for  `rbc_gem_utils`, a helper package designed to simplify working with the RBC-GEM repository content. Can be installed by following the provided [installation instructons](#installation).
