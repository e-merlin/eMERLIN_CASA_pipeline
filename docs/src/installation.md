# Installation

The e-MERLIN CASA Pipeline (eMCP) has the following dependencies:

- Python 3.8
- CASA v6.5+
- aoflagger v2.9+ (needed for L-band data)
- wsclean (optional, for improved imaging)

There are three ways to install eMCP, with conda being the recommended approach.

## Conda Installation (Recommended)

This is the recommended approach as it will install all dependencies, including non-Python ones like aoflagger and wsclean.

### Step 1: Install Miniconda (if not already installed)

```bash
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash ./Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc  # or restart your terminal
```

### Step 2: Clone the repository and create the environment

```bash
git clone https://github.com/e-merlin/eMERLIN_CASA_pipeline.git
cd eMERLIN_CASA_pipeline
conda env create -f environment.yml
conda activate emcp
```

This will create a new conda environment called `emcp` with all necessary dependencies installed.

## Pip Installation

If you already have CASA, aoflagger and wsclean installed on your system, you can install eMCP using pip:

```bash
pip install git+https://github.com/e-merlin/eMERLIN_CASA_pipeline.git
```

Or to install from a local copy:

```bash
git clone https://github.com/e-merlin/eMERLIN_CASA_pipeline.git
cd eMERLIN_CASA_pipeline
pip install .
```

## Docker Installation

For those who prefer containerized applications:

```bash
docker pull emerlin/emcp:latest
docker run -it --rm -v $(pwd):/data emerlin/emcp:latest
```

## Verifying Installation

After installation, you can verify that eMCP is installed correctly by running:

```bash
emcp -v
```

This should display the version of eMCP that's installed.
