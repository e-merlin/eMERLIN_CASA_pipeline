# Installation

You can install de dependencies with `conda`, and then install the pipeline with `pip`

If you don't have conda installed:

```bash
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash ./Miniconda3-latest-Linux-x86_64.sh
rm ./Miniconda3-latest-Linux-x86_64.sh
conda activate eMCP
# We recommend to install mamba:
conda install -c conda-forge mamba
```

Clone the pipeline repository:

```bash
git clone https://github.com/e-merlin/eMERLIN_CASA_pipeline.git
cd emerlin_casa_pipeline
pip install .
```

