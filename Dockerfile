# Base image
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies (union of all requirements)
RUN apt-get update && apt-get install -y \
    git \
    cmake \
    build-essential \
    g++ \
    pkg-config \
    casacore-data casacore-dev \
    libblas-dev liblapack-dev \
    liblua5.3-dev \
    libpython3-dev \
    libboost-date-time-dev \
    libboost-system-dev \
    libboost-test-dev \
    libboost-program-options-dev \
    libboost-filesystem-dev \
    libgtkmm-3.0-dev \
    libcfitsio-dev \
    libfftw3-dev \
    libgsl-dev \
    libhdf5-serial-dev \
    libhdf5-dev \
    libopenmpi-dev \
    libpng-dev \
    python3 \
    python3-dev \
    python3-pip \
    python3-numpy \
    python3-pytest \
    python3-sphinx \
    wget \
 && rm -rf /var/lib/apt/lists/*

# --- Install Miniconda and mamba ---
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh
ENV PATH="/opt/conda/bin:$PATH"

# Install mamba
RUN conda install -c conda-forge mamba -y

# --- Create conda environment for emcp ---
RUN mamba create -n emcp python=3.10 pip -c conda-forge -y

# Activate conda for following RUN commands
SHELL ["/bin/bash", "-c"]
ENV CONDA_DEFAULT_ENV=emcp
ENV PATH="/opt/conda/envs/emcp/bin:${PATH}:/opt/conda/bin:${PATH}"

# --- Install emcp from GitHub ---
RUN source /opt/conda/bin/activate emcp && \
    pip install git+https://github.com/e-merlin/eMERLIN_CASA_pipeline.git@casa6

# --- Build and install aoflagger ---
WORKDIR /opt
RUN git clone https://git.code.sf.net/p/aoflagger/code aoflagger-src
WORKDIR /opt/aoflagger-src
RUN mkdir build && cd build && cmake .. && make -j$(nproc) && make install

# Test aoflagger Python bindings
RUN cd /opt/aoflagger-src/python && echo "import aoflagger" | python3

# --- Build and install IDG (required for wsclean) ---
WORKDIR /opt
RUN git clone https://git.astron.nl/RD/idg.git
RUN mkdir -p /opt/idg/build && cd /opt/idg/build && cmake .. && make install -j$(nproc)

# --- Build and install wsclean ---
RUN git clone https://gitlab.com/aroffringa/wsclean.git
WORKDIR /opt/wsclean
RUN mkdir build && cd build && cmake .. && make -j$(nproc) && make install

# Test wsclean installation
RUN wsclean --version

# Set a working directory for data
WORKDIR /data

CMD ["/bin/bash"]

