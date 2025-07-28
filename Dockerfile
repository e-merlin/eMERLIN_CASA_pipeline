FROM amigahub/wsclean-dysco:v1

RUN apt-get update && \
    apt-get install -y --no-install-recommends libfuse2 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*


RUN pip3 install --upgrade pip
RUN pip3 install git+https://github.com/e-merlin/eMERLIN_CASA_pipeline.git@casa6
RUN mkdir -p /root/.casa/data

RUN emcp -h
