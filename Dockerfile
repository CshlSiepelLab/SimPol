FROM continuumio/miniconda3

COPY . /SimPolv2

RUN apt-get update && \
    apt-get -y install build-essential cmake gdb libhdf5-dev && \
    /opt/conda/bin/conda env create -f /SimPolv2/environment.yml && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

SHELL ["/bin/bash", "-c"]
ENV PATH /opt/conda/envs/$(head -n 1 environment.yml | cut -d' ' -f2)/bin:$PATH
