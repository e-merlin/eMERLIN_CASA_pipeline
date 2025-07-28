FROM ubuntu:22.04

# Set non-interactive frontend and install all dependencies in one layer
RUN export DEBIAN_FRONTEND=noninteractive && apt-get update && \
    apt-get install -y \
    git \
    cmake \
    build-essential \
    g++ \
    pkg-config \
    casacore-data casacore-dev \
    libblas-dev liblapack-dev \
    python3 \
    python3-pip \
    python3-dev \
    python3-numpy \
    python3-pytest \
    python3-sphinx \
    libpython3-dev \
    libboost-date-time-dev libboost-test-dev \
    libboost-program-options-dev libboost-system-dev libboost-filesystem-dev \
    libcfitsio-dev \
    libfftw3-dev \
    libgsl-dev \
    libhdf5-dev \
    libhdf5-serial-dev \
    libopenmpi-dev \
    libpng-dev \
    liblua5.3-dev \
    libgtkmm-3.0-dev \
    wget

# Install IDG (dependency for WSClean)
RUN mkdir /external && cd /external && git clone https://git.astron.nl/RD/idg.git && \
    mkdir /external/idg/build && cd /external/idg/build && cmake ../ && make install -j`nproc`

# Clone AOFlagger source
RUN cd /external && git clone https://gitlab.com/aroffringa/aoflagger.git

# Build and install AOFlagger
RUN cd /external/aoflagger && mkdir /build && cd /build && cmake ../src
RUN cd /external/aoflagger/build && make -j`nproc --all` && make install
RUN cd /external/aoflagger/build/python && echo "import aoflagger" | python3

# Clone WSClean source  
RUN cd /external && git clone https://gitlab.com/aroffringa/wsclean.git

# Build and install WSClean
RUN cd /external/wsclean && mkdir /build && cd /build && cmake ../src && \
    make -j`nproc` && make install && wsclean --version

# Install emcp from pip (eMERLIN CASA pipeline)
RUN pip3 install git+https://github.com/e-merlin/eMERLIN_CASA_pipeline.git@casa6

# Set working directory
WORKDIR /workspace
