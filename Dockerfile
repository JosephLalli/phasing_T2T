# Phasing T2T Project Docker Container
# Runs the phasing pipeline and includes the Python/R support needed for Figure 5.
# The published image is used by ./run_pipeline.sh and scripts/run_make_plots.sh.
# External inputs and generated outputs stay in the bind-mounted repository.

FROM python:3.11-bookworm@sha256:a8f8fbe1a0edc9e4dddafa64ba73f7e04be7be5ebc23f332362e779e0a2e4e52

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=UTC

WORKDIR /build

# ── System dependencies ──
# Microsoft core fonts are installed from Debian's contrib package so SVG/PDF
# generation resolves Arial consistently. The EULA is accepted noninteractively.
RUN sed -i 's/^Components: main$/Components: main contrib non-free non-free-firmware/' /etc/apt/sources.list.d/debian.sources && \
    echo 'msttcorefonts msttcorefonts/accepted-mscorefonts-eula select true' | debconf-set-selections && \
    apt-get update && apt-get install -y --no-install-recommends \
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
    libtiff-dev \
    libxml2-dev \
    python3-dev \
    r-base \
    r-base-dev \
    cmake \
    fontconfig \
    fonts-liberation2 \
    ttf-mscorefonts-installer \
    librsvg2-bin \
    poppler-utils \
    parallel \
    && fc-cache -f \
    && test "$(fc-match -f '%{family}' Arial)" = 'Arial' \
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
RUN python -m pip install --no-cache-dir -r /tmp/requirements.txt && rm /tmp/requirements.txt

# ── R packages for Figure 5 generation ──
# fs can build its bundled libuv without another system dependency.
RUN USE_BUNDLED_LIBUV=1 R -q -e "options(Ncpus=8); install.packages(c('BiocManager', 'arrow', 'PlotTools', 'tidyverse', 'svglite'), repos='https://cloud.r-project.org')" && \
    USE_BUNDLED_LIBUV=1 R -q -e "options(Ncpus=8); BiocManager::install(c('GenomicRanges', 'rtracklayer', 'karyoploteR'), ask=FALSE, update=FALSE); stopifnot(all(vapply(c('arrow', 'PlotTools', 'tidyverse', 'svglite', 'GenomicRanges', 'rtracklayer', 'karyoploteR'), requireNamespace, logical(1), quietly=TRUE)))"

# The resource downloader decompresses the GRCh38 reference with zstd.
RUN apt-get update && apt-get install -y --no-install-recommends zstd && \
    rm -rf /var/lib/apt/lists/*

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
LABEL version="2.1"

CMD ["/bin/bash"]
