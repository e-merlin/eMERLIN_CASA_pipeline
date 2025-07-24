FROM continuumio/miniconda3:latest

# Set environment variables to avoid prompts during install
ENV DEBIAN_FRONTEND=noninteractive

# Configure conda and install mamba
RUN conda config --add channels conda-forge && \
    conda config --set channel_priority strict && \
    conda install -y mamba && \
    conda clean --all -f -y

# Create and activate environment with exact dependencies
RUN mamba create -y -n emcp \
    astropy=6.1.7 \
    python=3.10.17 \
    casacore=3.7.1 \
    python-casacore=3.7.1 \
    cmasher=1.9.2 \
    ipython=8.36.0 \
    libboost-python=1.86.0 \
    matplotlib=3.10.3 \
    mpi4py=4.0.3 \
    numpy=2.2.6 \
    openmpi=5.0.7 \
    scipy=1.15.2 \
    pip=25.1.1 \
    setuptools=80.8.0 \
    setuptools-scm \
    wheel && \
    conda clean --all -f -y

# Activate environment and install pip-only dependencies
SHELL ["conda", "run", "-n", "emcp", "/bin/bash", "-c"]

RUN pip install \
    casaconfig==1.0.2 \
    casatools==6.7.0.31 \
    casatasks==6.7.0.31 \
    casaplotms==2.6.2 \
    casaviewer==2.3.2 \
    casashell==6.7.0.31 \
    casaplotserver==1.9.2 \
    casatestutils==6.7.0.31 \
    casatablebrowser==0.0.37 \
    casalogger==1.0.21 \
    casafeather==0.0.24 \
    casampi==0.5.6

# Set working directory and copy files
WORKDIR /app
COPY . /app/

# Install your package
RUN pip install -e .

# Entrypoint and default CMD
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "emcp", "emcp"]
CMD ["--help"]

# Metadata
LABEL org.opencontainers.image.source="https://github.com/e-merlin/eMERLIN_CASA_pipeline"
LABEL org.opencontainers.image.description="Container image for eMERLIN CASA pipeline"
LABEL org.opencontainers.image.licenses="GPL-3.0-or-later"

