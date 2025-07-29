FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive \
    LC_ALL=C.UTF-8 \
    QT_QPA_PLATFORM=offscreen

# Install all system dependencies for the three packages in a single layer
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    # Common build tools
    build-essential \
    cmake \
    g++ \
    git \
    pkg-config \
    wget \
    # Python dependencies
    python3 \
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
    # Clean up apt cache to reduce image size
    && apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# === Install ecmp (eMERLIN CASA Pipeline) ===
RUN pip3 install --no-cache-dir --upgrade pip
RUN pip3 install --no-cache-dir git+https://github.com/e-merlin/eMERLIN_CASA_pipeline.git@casa6
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

# === Install IDG (Dependency for WSClean) ===
RUN mkdir -p /src/idg && \
    git clone https://git.astron.nl/RD/idg.git /src/idg && \
    mkdir -p /src/idg/build && cd /src/idg/build && \
    cmake .. && \
    make -j`nproc` && \
    make install

# === Clone, Build, and Install WSClean ===
# Clones the default branch. For reproducibility, add `--branch <tag>`
RUN git clone https://gitlab.com/aroffringa/wsclean.git /src/wsclean
WORKDIR /src/wsclean
RUN mkdir -p build && cd build && \
    cmake .. && \
    make -j`nproc` && \
    make install

# === Clone, Build, and Install AOFlagger ===
# Clones the default branch. For reproducibility, add `--branch <tag>`
RUN git clone https://gitlab.com/aroffringa/aoflagger.git /src/aoflagger
WORKDIR /src/aoflagger
RUN mkdir -p build && cd build && \
    cmake .. && \
    make -j`nproc` && \
    make install

# === Final Cleanup & Verification ===
# Remove source and build directories to save space
RUN rm -rf /src
# Set a final working directory
WORKDIR /root
# Verify all installations
RUN echo "Verifying installations..." && \
    wsclean --version && \
    emcp -l && \
    python3 -c "import aoflagger; print('AOFlagger Python bindings imported successfully.')"
