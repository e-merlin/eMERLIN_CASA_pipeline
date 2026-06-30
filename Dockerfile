FROM ubuntu:24.04 AS native-builder

ENV \
    AOFLAGGER_TAG=v3.5.0 \
    CASACORE_TAG=v3.8.0 \
    DEBIAN_FRONTEND=noninteractive \
    WSCLEAN_TAG=v3.7

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      bison \
      ca-certificates \
      cmake \
      flex \
      g++ \
      gfortran \
      git \
      libblas-dev \
      libboost-dev \
      libboost-filesystem-dev \
      libboost-program-options-dev \
      libboost-python-dev \
      libcfitsio-dev \
      libfftw3-dev \
      libgsl-dev \
      libhdf5-dev \
      liblapack-dev \
      liblua5.3-dev \
      libopenmpi-dev \
      libpng-dev \
      libpython3-dev \
      libreadline-dev \
      make \
      pkg-config \
      pybind11-dev \
      python3 \
      python3-dev \
      wcslib-dev \
      wget && \
    rm -rf /var/lib/apt/lists/*

RUN mkdir -p /usr/share/casacore/data && \
    ln -s /usr/share/casacore /var/lib/casacore && \
    wget -qO - https://www.astron.nl/iers/WSRT_Measures.ztar | \
      tar -C /usr/share/casacore/data -xzf -

RUN mkdir -p /src && \
    cd /src && \
    git clone --depth 1 --branch "${CASACORE_TAG}" \
      https://github.com/casacore/casacore.git && \
    cmake -S casacore -B casacore/build \
      -DBUILD_TESTING=OFF \
      -DBUILD_PYTHON3=OFF \
      -DDATA_DIR=/usr/share/casacore/data && \
    cmake --build casacore/build --parallel "$(nproc)" && \
    cmake --install casacore/build && \
    rm -rf casacore && \
    \
    git clone --depth 1 --branch "${WSCLEAN_TAG}" --recurse-submodules --shallow-submodules \
      https://gitlab.com/aroffringa/wsclean.git && \
    cmake -S wsclean -B wsclean/build && \
    cmake --build wsclean/build --parallel "$(nproc)" && \
    cmake --install wsclean/build && \
    rm -rf wsclean && \
    wsclean --version && \
    \
    git clone --depth 1 --branch "${AOFLAGGER_TAG}" --recurse-submodules --shallow-submodules \
      https://gitlab.com/aroffringa/aoflagger.git && \
    cmake -S aoflagger -B aoflagger/build -DENABLE_GUI=OFF && \
    cmake --build aoflagger/build --parallel "$(nproc)" && \
    cmake --install aoflagger/build && \
    cd / && \
    if ! python3 -c "import importlib.util, sys; sys.exit(importlib.util.find_spec('aoflagger') is None)"; then \
      site_packages="$(python3 -c 'import sysconfig; print(sysconfig.get_paths()["platlib"])')" && \
      mkdir -p "${site_packages}" && \
      cp /src/aoflagger/build/python/aoflagger*.so "${site_packages}/"; \
    fi && \
    python3 -c "import aoflagger; print(aoflagger.__file__)" && \
    aoflagger --version && \
    rm -rf aoflagger

FROM ubuntu:24.04

ENV \
    CASASITECONFIG=/usr/local/etc/casasiteconfig.py \
    DEBIAN_FRONTEND=noninteractive \
    PYTHONNOUSERSITE=1

COPY --from=native-builder /usr/local /usr/local
COPY --from=native-builder /usr/share/casacore /usr/share/casacore

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      ca-certificates \
      libblas-dev \
      libboost-filesystem-dev \
      libboost-program-options-dev \
      libboost-python-dev \
      libcfitsio-dev \
      libfftw3-dev \
      libgsl-dev \
      libhdf5-dev \
      liblapack-dev \
      liblua5.3-dev \
      libopenmpi-dev \
      libpng-dev \
      libpython3-dev \
      libreadline-dev \
      python3 \
      python3-dev \
      python3-pip \
      wcslib-dev && \
    rm -rf /var/lib/apt/lists/* && \
    ln -s /usr/share/casacore /var/lib/casacore && \
    ldconfig && \
    python3 -m pip install --no-cache-dir --no-compile --break-system-packages \
      numpy \
      aplpy \
      astropy \
      astroquery \
      cmasher \
      ipython \
      matplotlib \
      scipy \
      reproject \
      pyregion \
      protobuf==3.20 \
      'shadems>=0.5.4,<0.6.0' \
      'dask>=2024.4.1,<2024.5.0' \
      casaconfig==1.4.0 \
      casatools==6.7.2.42 \
      casatasks==6.7.2.42 \
      casashell==6.7.2.42 \
      casatestutils==6.7.2.42 \
      casatablebrowser==0.0.39 \
      casalogger==1.0.23 \
      casampi==0.5.9

COPY pyproject.toml README.md /opt/eMERLIN_CASA_pipeline/
COPY src /opt/eMERLIN_CASA_pipeline/src
COPY .git /opt/eMERLIN_CASA_pipeline/.git

RUN apt-get update && \
    apt-get install -y --no-install-recommends git && \
    rm -rf /var/lib/apt/lists/* && \
    cd /opt/eMERLIN_CASA_pipeline && \
    python3 -m pip install --no-cache-dir --no-compile --break-system-packages \
      --no-deps . && \
    apt-get purge -y --auto-remove git git-man && \
    mkdir -p /root/.casa/data /usr/local/etc && \
    printf '%s\n' \
      'measurespath = "/root/.casa/data"' \
      'datapath = ["/root/.casa/data", "/usr/share/casacore/data"]' \
      'data_auto_update = False' \
      'measures_auto_update = False' \
      'nologger = True' \
      'nogui = True' \
      > "${CASASITECONFIG}" && \
    rm -rf /opt/eMERLIN_CASA_pipeline && \
    cd / && \
    python3 --version && \
    wsclean --version && \
    aoflagger --version && \
    python3 -c "import aoflagger; print(aoflagger.__file__)" && \
    python3 -c "import importlib.util; assert importlib.util.find_spec('casadata') is None" && \
    python3 -c "import eMCP; import casatasks; import casatools" && \
    emcp -h

WORKDIR /work
