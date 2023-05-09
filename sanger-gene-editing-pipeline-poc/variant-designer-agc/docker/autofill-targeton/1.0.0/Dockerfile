FROM ubuntu:23.04

# Install required packages
RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive \
    apt-get install -y \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libfontconfig1-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libfreetype6-dev \
        libpng-dev \
        libtiff5-dev \
        libjpeg-dev \
        r-base \
        procps && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R -e "install.packages(c('validate', 'glue', 'here', 'vroom', 'cli'), repos='http://cran.rstudio.com/')" && \
    R -e "install.packages(c('tidyverse'), repos='http://cran.rstudio.com/')" && \
    R -e "if (!require('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install('Biostrings')"
