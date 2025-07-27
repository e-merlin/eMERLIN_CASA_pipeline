FROM ubuntu:22.04

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:/usr/local/lib"
ENV PYTHONPATH="${PYTHONPATH}:/usr/local/lib/python3.10/site-packages"

# Install system dependencies for all three packages
RUN apt-get update && apt-get install -y \
    # Build tools
    git \
    cmake \
    build-essential \
    g++ \
    clang \
    gfortran \
    pkg-config \
    bison \
    flex \
    # Python and development
    python3.10 \
    python3.10-dev \
    python3.10-distutils \
    python3-pip \
    libpython3.10-dev \
    python3-numpy \
    python3-pytest \
    python3-sphinx \
    # Libraries for casacore and radio astronomy
    casacore-data \
    casacore-dev \
    wcslib-dev \
    # Math and science libraries
    libblas-dev \
    liblapack-dev \
    libfftw3-dev \
    libgsl-dev \
    libhdf5-dev \
    # Boost libraries
    libboost-date-time-dev \
    libboost-test-dev \
    libboost-program-options-dev \
    libboost-system-dev \
    libboost-filesystem-dev \
    # I/O libraries
    libcfitsio-dev \
    libpng-dev \
    # Parallel processing
    libopenmpi-dev \
    # C++ standard library for aoflagger
    libc++-dev \
    libc++abi-dev \
    # Lua for aoflagger
    liblua5.3-dev \
    # Utilities
    wget \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Create software directory
RUN mkdir -p /software /external

# Install conda/mamba for Python environment management
RUN curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -o miniforge.sh && \
    bash miniforge.sh -b -p /opt/conda && \
    rm miniforge.sh
ENV PATH="/opt/conda/bin:${PATH}"

# Install specific Boost version for aoflagger compatibility
WORKDIR /software
RUN wget -nv https://sourceforge.net/projects/boost/files/boost/1.78.0/boost_1_78_0.tar.bz2/download -O boost_1_78_0.tar.bz2 && \
    tar xjf boost_1_78_0.tar.bz2 && \
    cd boost_1_78_0/ && \
    ./bootstrap.sh --with-toolset=clang && \
    ./b2 toolset=clang cxxflags="-stdlib=libc++" linkflags="-stdlib=libc++" install && \
    cd .. && rm -rf boost_1_78_0 boost_1_78_0.tar.bz2

# Install HDF5 with clang for aoflagger compatibility
RUN wget -nv -O - https://github.com/HDFGroup/hdf5/archive/refs/tags/hdf5-1_12_1.tar.gz | tar xz && \
    cd hdf5-hdf5-1_12_1 && \
    CC=/usr/bin/clang CXX=/usr/bin/clang++ CXXFLAGS="-stdlib=libc++" LDFLAGS="-stdlib=libc++" \
    ./configure --prefix /usr/local --enable-cxx && \
    make -j$(nproc) install

# Install casacore from source for aoflagger compatibility
RUN mkdir -p casacore/build && \
    cd casacore && \
    wget -nv -O - https://github.com/casacore/casacore/archive/refs/tags/v3.4.0.tar.gz | tar xz && \
    cd build && \
    cmake \
      -DBUILD_TESTING=OFF \
      -DCMAKE_C_COMPILER=/usr/bin/clang \
      -DCMAKE_CXX_COMPILER=/usr/bin/clang++ \
      -DCMAKE_CXX_FLAGS="-stdlib=libc++" \
      -DBUILD_PYTHON=OFF \
      -DBUILD_PYTHON3=OFF \
      ../casacore-3.4.0 && \
    make install -j$(nproc)

# Install EveryBeam for wsclean
WORKDIR /external
RUN git clone https://git.astron.nl/RD/EveryBeam.git && \
    mkdir EveryBeam/build && \
    cd EveryBeam/build && \
    cmake ../ && \
    make install -j$(nproc)

# Install Dysco for wsclean
RUN git clone https://github.com/aroffringa/dysco.git && \
    mkdir dysco/build && \
    cd dysco/build && \
    cmake ../ && \
    make install -j$(nproc)

# Install IDG for wsclean
RUN git clone https://git.astron.nl/RD/idg.git && \
    mkdir idg/build && \
    cd idg/build && \
    cmake ../ && \
    make install -j$(nproc)

# Install wsclean
WORKDIR /software
RUN git clone https://gitlab.com/aroffringa/wsclean.git && \
    mkdir wsclean-build && \
    cd wsclean-build && \
    cmake ../wsclean && \
    make -j$(nproc) && \
    make install

# Install aoflagger
RUN git clone https://gitlab.com/aroffringa/aoflagger.git && \
    mkdir aoflagger-build && \
    cd aoflagger-build && \
    cmake ../aoflagger \
      -DCMAKE_C_COMPILER=/usr/bin/clang \
      -DCMAKE_CXX_COMPILER=/usr/bin/clang++ \
      -DCMAKE_CXX_FLAGS="-stdlib=libc++" && \
    make -j$(nproc) && \
    make install

# Test aoflagger Python bindings
RUN cd aoflagger-build/python && echo "import aoflagger" | python3

# Create Python 3.10 environment for eMERLIN CASA pipeline
RUN mamba create -n emcp python=3.10 pip -c conda-forge -y

# Activate the environment for subsequent commands
SHELL ["/bin/bash", "-c"]
RUN echo "source /opt/conda/bin/activate emcp" >> ~/.bashrc
ENV CONDA_DEFAULT_ENV=emcp
ENV PATH="/opt/conda/envs/emcp/bin:${PATH}"

# Install eMERLIN CASA pipeline from the casa6 branch
RUN source /opt/conda/bin/activate emcp && \
    pip install git+https://github.com/e-merlin/eMERLIN_CASA_pipeline.git@casa6

# Verify installations
RUN wsclean --version
RUN aoflagger --version
RUN source /opt/conda/bin/activate emcp && emcp -h

# Set working directory
WORKDIR /data

# Create entrypoint script to activate conda environment
RUN echo '#!/bin/bash' > /entrypoint.sh && \
    echo 'source /opt/conda/bin/activate emcp' >> /entrypoint.sh && \
    echo 'exec "$@"' >> /entrypoint.sh && \
    chmod +x /entrypoint.sh

ENTRYPOINT ["/entrypoint.sh"]

# Default command runs the eMERLIN CASA pipeline help
CMD ["emcp", "-h"]

# Metadata
LABEL org.opencontainers.image.source="https://github.com/e-merlin/eMERLIN_CASA_pipeline"
LABEL org.opencontainers.image.description="Container image for eMERLIN CASA pipeline"
LABEL org.opencontainers.image.licenses="GPL-3.0-or-later"

