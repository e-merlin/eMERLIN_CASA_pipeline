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
    libxml2-dev \  
 && rm -rf /var/lib/apt/lists/*

# Install Python 3.10 and emcp
RUN apt-get update && apt-get install -y software-properties-common && \
    add-apt-repository ppa:deadsnakes/ppa && \
    apt-get update && apt-get install -y python3.10 python3.10-dev python3.10-distutils && \
    rm -rf /var/lib/apt/lists/*

# Install pip for Python 3.10 and then emcp
RUN wget https://bootstrap.pypa.io/get-pip.py && \
    python3.10 get-pip.py && \
    rm get-pip.py && \
    python3.10 -m pip install git+https://github.com/e-merlin/eMERLIN_CASA_pipeline.git@casa6

# Build and install IDG (dependency for wsclean)
WORKDIR /external
RUN git clone https://git.astron.nl/RD/idg.git && \
    mkdir /external/idg/build && \
    cd /external/idg/build && \
    cmake ../ && \
    make install -j$(nproc)

# Build and install aoflagger
RUN apt-get update && apt-get install -y \
    libboost-python1.74-dev \
    libboost-numpy1.74-dev \
    python3.10-numpy && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /external
RUN git clone https://git.code.sf.net/p/aoflagger/code aoflagger-src
RUN mkdir /external/build && \
    cd /external/build && \
    cmake ../aoflagger-src -DPYTHON_EXECUTABLE=/usr/bin/python3.10 && \
    make -j$(nproc) && \
    make install && \
    cd /external/build/python && \
    echo "import aoflagger" | python3.10

# Build and install wsclean
WORKDIR /external
RUN git clone https://gitlab.com/aroffringa/wsclean.git wsclean-src
RUN mkdir /external/build-wsclean && \
    cd /external/build-wsclean && \
    cmake ../wsclean-src && \
    make -j$(nproc) && \
    make install && \
    wsclean --version

# Set working directory for user
WORKDIR /data
