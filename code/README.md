# RBC-GEM Code

This directory contains functions and scripts that facilitate work with the RBC-GEM model and other content in the repository. 

This directory is organized as follows:

```
code
├── notebooks   # For Python 3
└── src         # For Python 3
```

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

Brief description of directory contents are provided below.
### Notebooks (Python)

**TODO** 

Contains iPython (`.ipynb`) notebooks

### src (Python)

Contains the source code for  `rbc_gem_utils`, a helper package designed to simplify working with the RBC-GEM repository content. Can be installed by following the provided [installation instructons](#installation).