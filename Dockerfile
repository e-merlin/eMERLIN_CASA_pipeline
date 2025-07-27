FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# Install all system dependencies in one layer (union of all requirements)
RUN apt-get update && apt-get install -y \
    # Build tools
    git \
    cmake \
    build-essential \
    g++ \
    pkg-config \
    wget \
    # Python and pip
    python3 \
    python3-dev \
    python3-pip \
    python3-numpy \
    python3-pytest \
    python3-sphinx \
    # Core libraries for astronomical software
    casacore-data \
    casacore-dev \
    # Math libraries
    libblas-dev \
    liblapack-dev \
    libfftw3-dev \
    libgsl-dev \
    # I/O libraries
    libcfitsio-dev \
    libhdf5-serial-dev \
    libhdf5-dev \
    libpng-dev \
    # Boost libraries
    libboost-date-time-dev \
    libboost-system-dev \
    libboost-test-dev \
    libboost-program-options-dev \
    libboost-filesystem-dev \
    # GUI and other libraries
    libgtkmm-3.0-dev \
    liblua5.3-dev \
    libopenmpi-dev \
 && rm -rf /var/lib/apt/lists/*

# Install emcp directly with pip (simpler than conda environment)
RUN pip3 install git+https://github.com/e-merlin/eMERLIN_CASA_pipeline.git@casa6

# Build and install IDG (dependency for wsclean)
WORKDIR /external
RUN git clone https://git.astron.nl/RD/idg.git && \
    mkdir /external/idg/build && \
    cd /external/idg/build && \
    cmake ../ && \
    make install -j$(nproc)

# Build and install aoflagger (following official Dockerfile)
WORKDIR /external
RUN git clone https://git.code.sf.net/p/aoflagger/code aoflagger
WORKDIR /external/aoflagger
RUN mkdir /build-aoflagger && \
    cd /build-aoflagger && \
    cmake ../aoflagger && \
    make -j$(nproc) && \
    make install && \
    cd /build-aoflagger/python && \
    echo "import aoflagger" | python3

# Build and install wsclean (following official Dockerfile)
WORKDIR /external
RUN git clone https://gitlab.com/aroffringa/wsclean.git
WORKDIR /external/wsclean
RUN mkdir /build-wsclean && \
    cd /build-wsclean && \
    cmake ../wsclean && \
    make -j$(nproc) && \
    make install && \
    wsclean --version

# Set working directory for user
WORKDIR /data

# Default command
CMD ["/bin/bash"]
