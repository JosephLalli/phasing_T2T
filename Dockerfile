# Phasing T2T Project Docker Container
# Runs the phasing pipeline test cases (chr22_test, chr15_test).
# No R, Java, or GATK -- only bcftools/samtools/htslib + Python.
#
# ── Volume mount requirements ──
# The following MUST be bind-mounted at runtime (too large for the image):
#
#   /phasing_T2T_project/resources/chm13v2.0.fa.gz          (+.fai, +.gzi)
#   /phasing_T2T_project/resources/GRCh38_full_analysis_set_plus_decoy_hla.fa.gz (+.fai, +.gzi)
#   /phasing_T2T_project/unphased_variant_calls/            (t2t/ and grch38/ subdirs)
#   /phasing_T2T_project/resources/hprc-v1.1-mc-chm13.vcfbub.a100k.wave.vcf.gz (+.tbi)
#   /phasing_T2T_project/resources/hgsvc3-*                 (all hgsvc3 files)
#   /phasing_T2T_project/resources/SGDP_variation/          (t2t/ and grch38/ subdirs)
#   /phasing_T2T_project/phased_panels/                     (grch38/ subdir with precomputed GRCh38 panels)
#
# Example invocation:
#   docker run -v /data/refs:/phasing_T2T_project/resources \
#              -v /data/vcf:/phasing_T2T_project/unphased_variant_calls \
#              -v /data/panels:/phasing_T2T_project/phased_panels \
#              phasing-t2t bash -c "cd scripts && ./create_and_assess_haplotype_panels.sh chr22_test 12 test CHM13v2.0"

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
    python3 \
    python3-pip \
    python3-dev \
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
    scipy \
    scikit-learn \
    matplotlib \
    seaborn

# ── Set up project directory ──
WORKDIR /phasing_T2T_project

ENV PATH="/usr/local/bin:${PATH}"
ENV LD_LIBRARY_PATH="/usr/local/lib:${LD_LIBRARY_PATH}"

# ── COPY static binaries (SHAPEIT5 etc.) ──
COPY bin/ bin/

# ── COPY core pipeline scripts ──
COPY scripts/create_and_assess_haplotype_panels.sh scripts/
COPY scripts/run_and_assess_imputation.sh scripts/
COPY scripts/assess_imputation.sh scripts/
COPY scripts/parameters.sh scripts/
COPY scripts/find_trio_singletons.sh scripts/
COPY scripts/liftover_panel.sh scripts/
COPY scripts/liftover_indels.py scripts/
COPY scripts/get_discordant_multiallelic_sites.py scripts/
COPY scripts/generate_regions_for_rare_phasing.py scripts/
COPY scripts/create_syn_nonsyn_bins.sh scripts/
COPY scripts/calc_genomewide_imputation_statistics_full.sh scripts/
COPY scripts/create_summary_phasing_dataframes_polars_regional.py scripts/

# ── COPY small resource files ──
COPY resources/sample_subsets/ resources/sample_subsets/
COPY resources/pedigrees/ resources/pedigrees/
COPY resources/recombination_maps/t2t_native_scaled_maps/ resources/recombination_maps/t2t_native_scaled_maps/
COPY resources/recombination_maps/grch38/ resources/recombination_maps/grch38/
COPY resources/chm13v2-syntenic_to_hg38.bed resources/
COPY resources/hg38.GCA_009914755.4.synNet.summary.bed.gz resources/
COPY resources/chm13v2.0_cytobands_allchrs.bed resources/
COPY resources/grch38_cytobands_allchrs.bed resources/
COPY resources/chm13v2-grch38.chain resources/
COPY resources/grch38-chm13v2.chain resources/
COPY resources/grch38-chm13v2.sort.vcf.gz resources/
COPY resources/grch38-chm13v2.sort.vcf.gz.tbi resources/
COPY resources/chm13v2-grch38.sort.vcf.gz resources/
COPY resources/chm13v2-grch38.sort.vcf.gz.tbi resources/
COPY resources/regions.txt resources/
COPY resources/1000G_omni2.5.hg38.t2t-chm13-v2.0.biallelic.vcf.gz resources/
COPY resources/1000G_omni2.5.hg38.t2t-chm13-v2.0.biallelic.vcf.gz.tbi resources/
COPY resources/1000_genomes_meta.tsv resources/

LABEL description="T2T genomic variant phasing pipeline (no R/Java/GATK)"
LABEL version="2.0"

CMD ["/bin/bash"]
