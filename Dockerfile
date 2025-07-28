FROM amigahub/wsclean-dysco:v1

ENV LC_ALL=C.UTF-8 \
    QT_QPA_PLATFORM=offscreen

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        python3 \
        python3-pip \
        libgfortran5 \
        libstdc++6 \
        libqt5core5a \
        libqt5gui5 \
        libqt5widgets5 \
        xvfb && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* && \
    pip3 install --upgrade pip && \
    pip3 install git+https://github.com/e-merlin/eMERLIN_CASA_pipeline.git@casa6 && \
    mkdir -p /root/.casa/data && \
    emcp -l

# Extract plotms AppImage and modify the calls  
RUN PLOTMS_DIR=$(find /usr/local -name "casaplotms-x86_64.AppImage" -exec dirname {} \; | head -1) && \
    cd "$PLOTMS_DIR" && \
    ./casaplotms-x86_64.AppImage --appimage-extract && \
    rm casaplotms-x86_64.AppImage && \
    find squashfs-root -type d | xargs chmod 775 && \
    chmod +x squashfs-root/AppRun && \
    find /usr/local -name "*.py" -exec grep -l "casaplotms-x86_64.AppImage" {} \; | \
        xargs sed -i 's/casaplotms-x86_64.AppImage/squashfs-root\/AppRun/g'
