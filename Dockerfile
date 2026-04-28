# Phasing T2T Project Docker Container
# Runs the phasing pipeline test cases (chr22_test, chr15_test).
# Includes bcftools/samtools/htslib, Python, and the R stack needed for Figure 6.
# No Java or GATK.
#
# ── Volume mount requirements ──
# The following MUST be bind-mounted at runtime (too large for the image).
# They are intentionally not tracked in git and are excluded from the image
# via .dockerignore. Supply them as bind mounts to the paths listed below.
#
# Reference genomes:
#   -v <chm13_fa>:/phasing_T2T_project/resources/chm13v2.0.fa.gz:ro  (+.fai, +.gzi)
#   -v <grch38_fa>:/phasing_T2T_project/resources/GRCh38_full_analysis_set_plus_decoy_hla.fa.gz:ro (+.fai, +.gzi)
#
# Unphased variant calls:
#   -v <unphased_vcfs>:/phasing_T2T_project/unphased_variant_calls:ro
#
# Pangenome VCFs:
#   -v <hprc_chm13_vcf>:/phasing_T2T_project/resources/hprc-v1.1-mc-chm13.vcfbub.a100k.wave.vcf.gz:ro  (+.tbi)
#   -v <hprc_grch38_vcf>:/phasing_T2T_project/resources/hprc-v1.1-mc-grch38.vcfbub.a100k.wave.vcf.gz:ro (+.tbi)
#
# SGDP ground truth:
#   -v <sgdp_grch38_dir>:/phasing_T2T_project/resources/SGDP_variation/grch38:ro
#   -v <sgdp_t2t_dir>:/phasing_T2T_project/resources/SGDP_variation/t2t:ro
#
# Precomputed GRCh38 panels (read-write: assess_imputation.sh writes filtered BCFs here):
#   -v <grch38_panels>:/phasing_T2T_project/phased_panels/grch38
#
# Persistent working directories:
#   -v <output_dir>/working_directories:/phasing_T2T_project/working_directories
#   -v <output_dir>/intermediate_data:/phasing_T2T_project/intermediate_data
#   -v <output_dir>/imputation_statistics:/phasing_T2T_project/imputation_statistics
#   -v <output_dir>/SHAPEIT5_switch_output:/phasing_T2T_project/SHAPEIT5_switch_output
#
# See scripts/run_docker_smoke_test.sh plus docker.env.example for a complete example.

FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=UTC

WORKDIR /build

# ── System dependencies ──
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    autoconf \
    automake \
    libtool \
    pkg-config \
    wget \
    curl \
    git \
    gzip \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libssl-dev \
    libncurses5-dev \
    libtbb-dev \
    libcurl4-openssl-dev \
    libcairo2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libjpeg-dev \
    libpng-dev \
    libtiff5-dev \
    libxml2-dev \
    python3 \
    python3-pip \
    python3-dev \
    r-base \
    r-base-dev \
    cmake \
    fonts-liberation2 \
    parallel \
    && rm -rf /var/lib/apt/lists/*

# ── Build htslib 1.22 ──
RUN wget -q https://github.com/samtools/htslib/releases/download/1.22/htslib-1.22.tar.bz2 && \
    tar -xjf htslib-1.22.tar.bz2 && \
    cd htslib-1.22 && \
    ./configure --prefix=/usr/local && \
    make -j24 && \
    make install && \
    cd .. && rm -rf htslib-1.22*

# ── Build samtools 1.22 ──
RUN wget -q https://github.com/samtools/samtools/releases/download/1.22/samtools-1.22.tar.bz2 && \
    tar -xjf samtools-1.22.tar.bz2 && \
    cd samtools-1.22 && \
    ./configure --prefix=/usr/local && \
    make -j24 && \
    make install && \
    cd .. && rm -rf samtools-1.22*

# ── Build bcftools 1.22 with liftover plugin ──
RUN wget -q https://github.com/samtools/bcftools/releases/download/1.22/bcftools-1.22.tar.bz2 && \
    tar -xjf bcftools-1.22.tar.bz2 && \
    cd bcftools-1.22 && \
    wget -q -P plugins https://raw.githubusercontent.com/freeseek/score/master/liftover.c && \
    ./configure --prefix=/usr/local && \
    make -j24 && \
    make install && \
    cd .. && rm -rf bcftools-1.22*

# Verify liftover plugin is available
RUN bcftools +liftover --help 2>&1 | head -1

# ── Python packages ──
COPY requirements.txt /tmp/requirements.txt
RUN pip3 install --no-cache-dir -r /tmp/requirements.txt && rm /tmp/requirements.txt

# ── R packages for Figure 6 generation ──
RUN R -q -e "install.packages(c('BiocManager', 'arrow', 'PlotTools', 'tidyverse', 'svglite'), repos='https://cloud.r-project.org')" && \
    R -q -e "BiocManager::install(c('GenomicRanges', 'rtracklayer', 'karyoploteR'), ask=FALSE, update=FALSE)"

# ── Set up project directory ──
WORKDIR /phasing_T2T_project

ENV PATH="/usr/local/bin:${PATH}"
ENV LD_LIBRARY_PATH="/usr/local/lib:${LD_LIBRARY_PATH}"

# ── COPY static binaries (SHAPEIT5 etc.) ──
COPY bin/ bin/

# ── COPY core pipeline scripts ──
COPY scripts/ scripts/

# ── COPY small resource files ──
COPY resources/ resources/

# ── COPY Jupyter notebooks for figure generation ──
COPY notebooks/ notebooks/

# ── Create mount-point directories for bind-mounted data ──
# These directories are empty in the image; fill them with -v at runtime.
RUN mkdir -p \
    unphased_variant_calls/t2t \
    unphased_variant_calls/grch38 \
    phased_panels/grch38 \
    working_directories \
    intermediate_data \
    imputation_statistics \
    SHAPEIT5_switch_output \
    figures \
    tables \
    notebook_runs \
    resources/SGDP_variation/t2t \
    resources/SGDP_variation/grch38

LABEL description="T2T genomic variant phasing pipeline with Python and R analysis support"
LABEL version="2.0"

CMD ["/bin/bash"]
