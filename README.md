# Verkko-fillet : post-Verkko graph and assembly cleaning in Python

Please post an issue before running the consensus again. I need to check whether your folder structure was generated correctly before version 1.0. Thanks!

verkko-fillet is an easy-to-use, Python-based toolkit for cleaning graph paths generated from [Verkko](https://github.com/marbl/verkko) assembler. It is designed to be run within a Jupyter notebook to enable interactive use. Verkko-fillet includes tools for performing assembly quality checks, identifying and resolving gaps, assigning chromosomes, and generating a corrected path (in a GAF-like format), as is required to generate a Verkko consensus run.

This Python-based implementation streamlines the entire process, starting right after the Verkko assembly is completed and preparing for the CNS run.


### Installation

#### Other tools
The gap-filling steps are easier when done interactively, such as in Jupyter Notebook or JupyterLab. Please install one and add the environment you generate to the Jupyter kernel to enable use.
During gap filling, we highly recommend viewing the graph, coverage, and trio markers if available using [BandageNG](https://github.com/asl/BandageNG). Please install this as well.

Dependencies (will not be installed using the command below):
* [gfacpp](https://github.com/snurk/gfacpp)

#### Install Verkko-fillet and other dependencies
Using `pip` is recommended. [link](https://pypi.org/project/verkkofillet/)

The default name of the Mamba or Conda environment is `verkko-fillet`. If you want to use a different name, please update the name field in the `environment.yaml` file before proceeding. All required packages and external tools are installed using the `vf_environment.Jun252025.yml` file. This file specifies the exact tools and their versions, ensuring a reproducible environment setup.

```
# Generate a mamba or conda environment and install dependencies.
mamba create -n verkko-fillet -f vf_environment.Jun252025.yml -vvv --dry-run --channel-priority flexible

# Once you’ve checked the generated environment, please re-run this without the --dry-run option. It may take some time.
mamba create -n verkko-fillet -f vf_environment.Jun252025.yml -vvv --channel-priority flexible

# Acivate the environment.
mamba activate verkko-fillet # or the name you desired

# Add python jupyter kernel.
python -m ipykernel install --user --name verkko-fillet --display-name verkko-fillet
pip install verkkofillet
```

