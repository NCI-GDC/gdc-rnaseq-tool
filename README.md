# gdc-rnaseq-tool

Utility scripts for GDC RNA-seq workflows. The docker file also installs Trimmomatic and fqvendorfail.

## Installation

```sh
pip install .
```

## Development

* Clone this repository
* Requirements:
  * Python >= 3.8
  * Tox
* `make venv` to create a virtualenv
* `source .venv/bin/activate` to activate new virtualenv
* `make init` to install dependencies and pre-commit hooks
