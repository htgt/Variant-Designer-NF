FROM ubuntu:20.04

# Install required packages
RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive \
    apt-get install -y \
        git \
        perl \
        libbio-perl-perl \
        build-essential \
        cpanminus && \
        procps && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install Perl packages
RUN cpanm Getopt::Long && \
    cpanm Bio::Seq && \
    cpanm Bio::SeqIO && \
    cpanm Bio::Tools::IUPAC && \
    cpanm Bio::Tools::SeqPattern

# clone Dispr repository and set the PATH environment variable
# to include the Dispr executable
RUN git clone https://github.com/douglasgscofield/dispr.git
ENV PATH="/dispr:${PATH}"
