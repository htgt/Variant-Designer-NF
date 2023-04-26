FROM ubuntu:20.04

# install dependencies
RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive \
    apt-get install -y \
        zlib1g-dev \
        libbz2-dev \
        libssl-dev \
        liblzma-dev \
        libcurl4-openssl-dev \
        procps \
        wget \
        gcc \
        make && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# install bcftools
RUN wget https://github.com/samtools/bcftools/releases/download/1.17/bcftools-1.17.tar.bz2 && \
    tar xjf bcftools-1.17.tar.bz2 && \
    cd bcftools-1.17 && \
    ./configure && \
    make && \
    make install
