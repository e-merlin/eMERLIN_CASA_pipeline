FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

# Base system
RUN apt-get update && apt-get install -y \
    software-properties-common \
    wget \
    curl \
    lsb-release \
    ca-certificates \
    gnupg

# Build tools
RUN apt-get install -y \
    cmake \
    g++ \
    git \
    pkg-config \
    libblas-dev  \
    libboost-date-time-dev \
    libboost-filesystem-dev \
    libboost-program-options-dev \
    libboost-system-dev \
    libcfitsio-dev \
    libfftw3-dev \
    libgsl-dev \
    libhdf5-dev \
    liblapack-dev \
    libopenmpi-dev \
    libpython3-dev \
    pkg-config \
    wcslib-dev \
    bison \
    flex \
    gfortran \
    python3-dev \
    python3-numpy \
    libboost-python-dev \
    libblas-dev \
    liblapack-dev \
    libreadline-dev \
        python3-dev \
    python3-pip \
    python3-numpy \
    python3-pytest \
    python3-sphinx \
    # ecmp dependencies
    libgfortran5 \
    libqt5core5a \
    libqt5gui5 \
    libqt5widgets5 \
    xvfb \
    # wsclean & aoflagger dependencies
    casacore-data \
    casacore-dev \
    libblas-dev \
    libboost-date-time-dev \
    libboost-filesystem-dev \
    libboost-program-options-dev \
    libboost-system-dev \
    libboost-test-dev \
    libcfitsio-dev \
    libfftw3-dev \
    libgsl-dev \
    libgtkmm-3.0-dev \
    libhdf5-dev \
    liblapack-dev \
    liblua5.3-dev \
    libopenmpi-dev \
    libpng-dev \
    wget && \
  rm -rf /var/lib/apt/lists/*

# Install the casacore measures data. We purposely do not install these from
# the Ubuntu repository, but download the latest version directly from the
# ASTRON ftp site.
# Note: The file on the ftp site is updated daily. When warnings regarding
# leap seconds appear, ignore them or regenerate the docker image.
RUN mkdir -p /usr/share/casacore/data && \
    ln -s /usr/share/casacore /var/lib/casacore && \
    wget -qO - https://www.astron.nl/iers/WSRT_Measures.ztar | \
        tar -C /usr/share/casacore/data -xzf -

# The casacore version in Ubuntu is too old to support C++20, so install a more recent one.
RUN mkdir /external && \
  cd /external && \
  git clone https://github.com/casacore/casacore.git && \
  cd /external/casacore && \
  git checkout ${CASACORE_VERSION} && \
  mkdir build && \
  cd build && \
  cmake .. -DBUILD_TESTING=OFF -DDATA_DIR=/usr/share/casacore/data && \
  make -j`nproc` && \
  make install -j`nproc` && \
  cd /external && \
  rm -rf /external/casacore


# Install EveryBeam
RUN cd /external && git clone https://git.astron.nl/RD/EveryBeam.git && \
  mkdir /external/EveryBeam/build && cd /external/EveryBeam/build && \
  cmake ../ && make install -j`nproc` && rm -rf /external/EveryBeam

# Install IDG
RUN cd /external && git clone https://git.astron.nl/RD/idg.git && \
  mkdir /external/idg/build && cd /external/idg/build && cmake ../ && \
  make install -j`nproc` && rm -rf /external/idg

# WSClean
RUN mkdir /src/ && cd /src/ && git clone https://gitlab.com/aroffringa/wsclean.git
WORKDIR /src/wsclean

RUN \
  mkdir /src/build && \
  cd /src/build && \
  cmake ../wsclean && \
  make -j`nproc` && \
  make install && \
  cd / && \
  rm -rf /src/build && \
  wsclean --version

# === Install ecmp (eMERLIN CASA Pipeline) ===
RUN pip install --no-cache-dir git+https://github.com/e-merlin/eMERLIN_CASA_pipeline.git@casa6 --break-system-packages
RUN mkdir -p /root/.casa/data

# Extract plotms AppImage included with ecmp and patch the call paths
RUN PLOTMS_DIR=$(find /usr/local -name "casaplotms-x86_64.AppImage" -exec dirname {} \; | head -1) && \
    cd "$PLOTMS_DIR" && \
    ./casaplotms-x86_64.AppImage --appimage-extract && \
    rm casaplotms-x86_64.AppImage && \
    find squashfs-root -type d | xargs chmod 775 && \
    chmod +x squashfs-root/AppRun && \
    find /usr/local -name "*.py" -exec grep -l "casaplotms-x86_64.AppImage" {} \; | \
    xargs sed -i 's/casaplotms-x86_64.AppImage/squashfs-root\/AppRun/g'

# === Clone, Build, and Install AOFlagger ===
# Clones the default branch. For reproducibility, add `--branch <tag>
RUN pip install pybind11 --break-system-packages
RUN  git clone https://gitlab.com/aroffringa/aoflagger.git /src/aoflagger
WORKDIR /src/aoflagger
RUN mkdir /buildao && cd /buildao && cmake -Dpybind11_DIR=$(python3.12 -m pybind11 --cmakedir)  ../src/aoflagger
RUN cd /buildao && make -j`nproc --all` && make install
# Copy aoflagger+pybind11 library to python packages folder
RUN cp /buildao/python/aoflagger*.so /usr/local/lib/python3.12/dist-packages/

# === Final Cleanup & Verification ===
# Remove source and build directories to save space
RUN rm -rf /src
# Set a final working directory
WORKDIR /root
# Verify all installations
RUN echo "Verifying installations..." && \
    wsclean --version && \
    #emcp -l && \
    python3 -c "import aoflagger; print('AOFlagger Python bindings imported successfully.')"
