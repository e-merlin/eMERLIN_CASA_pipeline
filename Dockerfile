FROM continuumio/miniconda3:latest

# Configure conda channels
RUN conda config --add channels conda-forge \
    && conda config --add channels pkgw-forge \
    && conda config --add channels i4ds \
    && conda config --set channel_priority strict

# Install base conda dependencies
RUN conda install -y python=3.8 mamba

# Install CASA and radio astronomy tools
RUN mamba install -y \
    casacore \
    openmpi \
    mpi4py \
    aoflagger \
    wsclean \
    numpy \
    pandas \
    matplotlib \
    scipy \
    astropy \
    python-casacore

# Install Python dependencies
RUN mamba install -y \
    setuptools>=62.6.0 \
    setuptools-scm \
    pip

# Clean up conda cache
RUN conda clean --all -f -y

# Copy the package files
WORKDIR /app
COPY . /app/

# Install the package
RUN pip install -e .

# Create entry points
ENTRYPOINT ["emcp"]

# Default command
CMD ["--help"]

RUN python3 --version
RUN pip3 --version
RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install mpi4py --no-cache-dir

LABEL org.opencontainers.image.source="https://github.com/e-merlin/eMERLIN_CASA_pipeline"
LABEL org.opencontainers.image.description="Container image for eMERLIN CASA pipeline"
LABEL org.opencontainers.image.licenses=GPL3
