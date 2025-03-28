# mgnify-pipelines-toolkit

This Python package contains a collection of scripts and tools for including in MGnify pipelines. Scripts stored here are mainly for:

- One-off production scripts that perform specific tasks in pipelines
- Scripts that have few dependencies
- Scripts that don't have existing containers built to run them
- Scripts for which building an entire container would be too bulky of a solution to deploy in pipelines

This package is built and uploaded to PyPi and bioconda. The package bundles scripts and makes them executable from the command-line when this package is installed.

## How to install

This package is available both on [PyPi](https://pypi.org/project/mgnify-pipelines-toolkit/) and bioconda.

To install from PyPi with pip:

`pip install mgnify-pipelines-toolkit`

To install from bioconda with conda/mamba:

`conda install -c bioconda mgnify-pipelines-toolkit`

You should then be able to run the packages from the command-line. For example to run the `get_subunits.py` script:

`get_subunits -i ${easel_coords} -n ${meta.id}`

## Adding a new script to the package

### Local development requirements
Before starting any development, you should do these few steps:
- Clone the repo if you haven't already and create a feature branch from the `dev` branch (NOT `main`).
- Create a virtual environment with the tool of your choice (i.e. `conda create --name my_new_env`)
- Activate you new environment (i.e. `conda activate my_new_env`)
- Install dev dependencies `pip install -e '.[tests,dev]'`
- Install pre-commit hooks `pre-commit install`
- Run unit tests `pytest`

When doing these steps above, you ensure that the code you add will be linted and formatted properly.

### New script requirements

There are a few requirements for your script:

- It needs to have a named main function of some kind. See `mgnify_pipelines_toolkit/analysis/shared/get_subunits.py` and the `main()` function for an example
- Because this package is meant to be run from the command-line, make sure your script can easily pass arguments using tools like `argparse` or `click`
- A small amount of dependencies. This requirement is subjective, but for example if your script only requires a handful of basic packages like `Biopython`, `numpy`, `pandas`, etc., then it's fine. However if the script has a more extensive list of dependencies, a container is probably a better fit.

### How to add a new script

To add a new Python script, first copy it over to the `mgnify_pipelines_toolkit` directory in this repository, specifically to the subdirectory that makes the most sense. If none of the subdirectories make sense for your script, create a new one. If your script doesn't have a `main()` type function yet, write one.

Then, open `pyproject.toml` as you will need to add some bits. First, add any missing dependencies (include the version) to the `dependencies` field.

Then, if you created a new subdirectory to add your script in, go to the `packages` line under `[tool.setuptools]` and add the new subdirectory following the same syntax.

Then, scroll down to the `[project.scripts]` line. Here, you will create an alias command for running your script from the command-line. In the example line:

`get_subunits = "mgnify_pipelines_toolkit.analysis.shared.get_subunits:main"`

- `get_subunits` is the alias
- `mgnify_pipelines_toolkit.analysis.shared.get_subunits` will link the alias to the script with the path `mgnify_pipelines_toolkit/analysis/shared/get_subunits.py`
- `:main` will specifically call the function named `main()` when the alias is run.

When you have setup this command, executing `get_subunits` on the command-line will be the equivalent of doing:

`from mgnify_pipelines_toolkit.analysis.shared.get_subunits import main; main()`

You should then write at least one unit test for your addition. This package uses `pytest` at the moment for this purpose. A GitHub Action workflow will run all of the unit tests whenever a commit is pushed to any branch.

Finally, you will need to bump up the version in the `version` line.

At the moment, these should be the only steps required to setup your script in this package (which is subject to change).

### Building and uploading to PyPi

The building and pushing of the package is automated by GitHub Actions, which will activate only on a new release. Bioconda should then automatically pick up the new PyPi release and push it to their recipes, though it's worth keeping an eye on their automated pull requests just in case [here](https://github.com/bioconda/bioconda-recipes/pulls).
