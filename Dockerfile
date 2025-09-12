# Phasing T2T Project Docker Container
# This container includes all dependencies for running the genomic variant phasing pipeline scripts

FROM ubuntu:22.04

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=UTC

# Set working directory
WORKDIR /phasing_T2T_project

# Install system dependencies
RUN apt-get update && apt-get install -y \
    # Core system tools
    # bash \
    coreutils \
    findutils \
    grep \
    sed \
    wget \
    curl \
    git \
    build-essential \
    cmake \
    pkg-config \
    # Compression and file handling
    gzip \
    apt-utils

RUN apt-get install -y \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    # libcurl4-openssl-dev \
    # libcurl4-gnutls-dev \
    libssl-dev \
    # Python and R dependencies
    python3 \
    python3-pip \
    python3-dev \
    r-base \
    r-base-dev \
    # Java for GATK
    openjdk-11-jdk \
    # Additional libraries
    libhts-dev \
    libncurses5-dev \
    libtbb-dev \
    autoconf \
    automake \
    libtool \
    # Parallel processing
    parallel \
    # Docker support (for scripts that use Docker)
    docker.io \
    && rm -rf /var/lib/apt/lists/*

# Install bcftools from source (version 1.22) with LiftOver plugin
RUN wget https://github.com/samtools/bcftools/releases/download/1.22/bcftools-1.22.tar.bz2 && \
    tar -xjf bcftools-1.22.tar.bz2
RUN cd bcftools-1.22 && \
    wget -P plugins http://raw.githubusercontent.com/freeseek/score/master/liftover.c && \
    ./configure --prefix=/usr/local && \
    make -j24 && \
    make install
RUN bcftools +liftover --help

# Install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.22/samtools-1.22.tar.bz2 && \
    tar -xjf samtools-1.22.tar.bz2 && \
    cd samtools-1.22 && \
    ./configure --prefix=/usr/local && \
    make -j24 && \
    make install && \
    cd .. && \
    rm -rf samtools-1.22*

# Install htslib
RUN wget https://github.com/samtools/htslib/releases/download/1.22/htslib-1.22.tar.bz2 && \
    tar -xjf htslib-1.22.tar.bz2 && \
    cd htslib-1.22 && \
    ./configure --prefix=/usr/local && \
    make -j24 && \
    make install && \
    cd .. && \
    rm -rf htslib-1.22*

# Install Python packages
RUN pip3 install --no-cache-dir \
    polars \
    pandas \
    numpy \
    pysam \
    cyvcf2 \
    intervaltree \
    tqdm \
    pyliftover \
    biopython \
    argparse \
    matplotlib \
    seaborn \
    scipy \
    scikit-learn

# Install R packages
RUN R -e "install.packages(c('karyoploteR', 'arrow', 'tidyverse', 'GenomicRanges', 'rtracklayer', 'magick'), repos='https://cran.rstudio.com/', dependencies=TRUE)"

# Install BiocManager and Bioconductor packages
RUN R -e "install.packages('BiocManager', repos='https://cran.rstudio.com/')" && \
    R -e "BiocManager::install(c('GenomicRanges', 'rtracklayer', 'karyoploteR'))"

# Install PlotTools R package (if available on CRAN, otherwise skip)
RUN R -e "tryCatch(install.packages('PlotTools', repos='https://cran.rstudio.com/'), error=function(e) cat('PlotTools not available on CRAN\n'))"

# Set up environment variables
ENV PATH="/usr/local/bin:${PATH}"
ENV LD_LIBRARY_PATH="/usr/local/lib:${LD_LIBRARY_PATH}"

# Default command
CMD ["/bin/bash"]

# Add labels for documentation
LABEL description="Docker container for T2T genomic variant phasing pipeline"
LABEL version="1.0"
LABEL maintainer="Phasing T2T Project"

# Health check to verify bcftools installation
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD bcftools --version || exit 1