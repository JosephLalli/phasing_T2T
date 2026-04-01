#!/usr/bin/env bash
# Download large external resources required by the phasing pipeline.
# Run from the repository root:  bash scripts/utility/download_resources.sh [OPTIONS]
#
# Files that already exist are skipped (wget -nc).  Re-run safely at any time.
#
# Options:
#   --chromosomes CHR_LIST   Comma-separated chromosomes to download (default: all)
#                            Example: --chromosomes 15,22  (for test regions)
#   --skip-unphased          Skip downloading unphased variant calls
#   --skip-grch38-panels     Skip downloading GRCh38 phased panels
#   --help                   Show this help
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
RESOURCES="${REPO_ROOT}/resources"

# ── Parse arguments ───────────────────────────────────────────────────────────
CHROMOSOMES=""
SKIP_UNPHASED=false
SKIP_GRCH38_PANELS=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --chromosomes)
            CHROMOSOMES="$2"
            shift 2
            ;;
        --skip-unphased)
            SKIP_UNPHASED=true
            shift
            ;;
        --skip-grch38-panels)
            SKIP_GRCH38_PANELS=true
            shift
            ;;
        --help|-h)
            sed -n '2,/^set /{ /^#/s/^# \?//p }' "$0"
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            exit 2
            ;;
    esac
done

# Build chromosome list
if [[ -n "${CHROMOSOMES}" ]]; then
    IFS=',' read -ra CHR_ARRAY <<< "${CHROMOSOMES}"
else
    CHR_ARRAY=()
    for i in $(seq 1 22); do CHR_ARRAY+=("$i"); done
    CHR_ARRAY+=("X")
fi

# ── Create directory structure ────────────────────────────────────────────────
mkdir -p \
    "${RESOURCES}" \
    "${RESOURCES}/SGDP_variation/t2t" \
    "${RESOURCES}/SGDP_variation/grch38" \
    "${REPO_ROOT}/unphased_variant_calls/t2t" \
    "${REPO_ROOT}/unphased_variant_calls/grch38" \
    "${REPO_ROOT}/phased_panels/grch38"

# ── Helper ────────────────────────────────────────────────────────────────────
fetch() {
    local dest_dir="$1" url="$2"
    local filename
    filename=$(basename "${url}")
    if [[ -f "${dest_dir}/${filename}" ]]; then
        echo "    Already exists: ${filename}"
        return 0
    fi
    wget -nc -P "${dest_dir}" "${url}" || {
        echo "WARNING: failed to download ${url}" >&2
    }
}

# Download a file and rename it to a different local name
fetch_as() {
    local dest_path="$1" url="$2"
    if [[ -f "${dest_path}" ]]; then
        echo "    Already exists: $(basename "${dest_path}")"
        return 0
    fi
    local tmpfile="${dest_path}.downloading"
    wget -q -O "${tmpfile}" "${url}" && mv "${tmpfile}" "${dest_path}" || {
        rm -f "${tmpfile}"
        echo "WARNING: failed to download ${url}" >&2
    }
}

# ── Reference genomes ─────────────────────────────────────────────────────────
echo "==> Downloading reference genomes..."

fetch "${RESOURCES}" \
    "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/GCA_009914755.4/chm13v2.0.fa.gz"
fetch "${RESOURCES}" \
    "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/GCA_009914755.4/chm13v2.0.fa.gz.fai"

# GRCh38 FASTA requires bgzip on the fly; only download if missing.
if [[ ! -f "${RESOURCES}/GRCh38_full_analysis_set_plus_decoy_hla.fa.gz" ]]; then
    echo "    Downloading and bgzipping GRCh38 reference (this may take a while)..."
    curl -f "https://42basepairs.com/download/s3/1000genomes/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa" \
        | bgzip > "${RESOURCES}/GRCh38_full_analysis_set_plus_decoy_hla.fa.gz"
else
    echo "    GRCh38 FASTA already present, skipping."
fi
# Index reference FASTAs (creates .fa.gz.fai and .fa.gz.gzi)
# Note: the upstream .fa.fai from 42basepairs is for uncompressed FASTA and
# cannot be used with the bgzipped version. samtools faidx must re-index.
if command -v samtools &>/dev/null; then
    for fasta in "${RESOURCES}/chm13v2.0.fa.gz" "${RESOURCES}/GRCh38_full_analysis_set_plus_decoy_hla.fa.gz"; do
        if [[ -f "${fasta}" && ! -f "${fasta}.fai" ]]; then
            echo "    Indexing $(basename "${fasta}") (generates .fai and .gzi)..."
            samtools faidx "${fasta}"
        fi
    done
else
    echo "WARNING: samtools not found. Reference FASTAs will need manual indexing:" >&2
    echo "  samtools faidx resources/chm13v2.0.fa.gz" >&2
    echo "  samtools faidx resources/GRCh38_full_analysis_set_plus_decoy_hla.fa.gz" >&2
fi

# ── Unphased variant calls ───────────────────────────────────────────────────
if [[ "${SKIP_UNPHASED}" == false ]]; then
    echo ""
    echo "==> Downloading unphased variant calls..."

    T2T_UNPHASED_BASE="https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/variants/1000_Genomes_Project/chm13v2.0/all_samples_3202"
    GRCH38_UNPHASED_BASE="https://42basepairs.com/download/s3/1000genomes/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot"

    for chr in "${CHR_ARRAY[@]}"; do
        echo "  chr${chr}..."
        # T2T-CHM13 calls
        fetch "${REPO_ROOT}/unphased_variant_calls/t2t" \
            "${T2T_UNPHASED_BASE}/1KGP.CHM13v2.0.chr${chr}.recalibrated.snp_indel.pass.vcf.gz"
        fetch "${REPO_ROOT}/unphased_variant_calls/t2t" \
            "${T2T_UNPHASED_BASE}/1KGP.CHM13v2.0.chr${chr}.recalibrated.snp_indel.pass.vcf.gz.tbi"

        # GRCh38 calls
        fetch "${REPO_ROOT}/unphased_variant_calls/grch38" \
            "${GRCH38_UNPHASED_BASE}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.recalibrated_variants.vcf.gz"
        fetch "${REPO_ROOT}/unphased_variant_calls/grch38" \
            "${GRCH38_UNPHASED_BASE}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.recalibrated_variants.vcf.gz.tbi"
    done
else
    echo ""
    echo "==> Skipping unphased variant calls (--skip-unphased)"
fi

# ── GRCh38 phased reference panels (Byrska-Bishop et al. 2022) ──────────────
# The upstream files from 1KGP use a different naming convention than the pipeline
# expects. This section downloads and renames them to match.
if [[ "${SKIP_GRCH38_PANELS}" == false ]]; then
    echo ""
    echo "==> Downloading GRCh38 phased panels (Byrska-Bishop et al. 2022)..."

    GRCH38_PANEL_BASE="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV"
    GRCH38_PANEL_DIR="${REPO_ROOT}/phased_panels/grch38"

    for chr in "${CHR_ARRAY[@]}"; do
        upstream_name="1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
        local_name="1KGP.GRCh38.chr${chr}.recalibrated.snp_indel.pass.phased.3202.vcf.gz"

        echo "  chr${chr}: ${upstream_name} -> ${local_name}"
        fetch_as "${GRCH38_PANEL_DIR}/${local_name}" \
            "${GRCH38_PANEL_BASE}/${upstream_name}"

        # Download tabix index and rename to match
        fetch_as "${GRCH38_PANEL_DIR}/${local_name}.tbi" \
            "${GRCH38_PANEL_BASE}/${upstream_name}.tbi"
    done
else
    echo ""
    echo "==> Skipping GRCh38 phased panels (--skip-grch38-panels)"
fi

# ── HPRC v1.1 pangenome VCFs ─────────────────────────────────────────────────
echo ""
echo "==> Downloading HPRC v1.1 pangenome VCFs..."

HPRC_BASE="https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus"

fetch "${RESOURCES}" "${HPRC_BASE}/hprc-v1.1-mc-chm13/hprc-v1.1-mc-chm13.vcfbub.a100k.wave.vcf.gz"
fetch "${RESOURCES}" "${HPRC_BASE}/hprc-v1.1-mc-chm13/hprc-v1.1-mc-chm13.vcfbub.a100k.wave.vcf.gz.tbi"

fetch "${RESOURCES}" "${HPRC_BASE}/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.vcfbub.a100k.wave.vcf.gz"
fetch "${RESOURCES}" "${HPRC_BASE}/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.vcfbub.a100k.wave.vcf.gz.tbi"

# ── HGSVC3 pangenome VCFs ────────────────────────────────────────────────────
echo ""
echo "==> Downloading HGSVC3 pangenome VCFs..."

HGSVC3_BASE="https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/scratch/2024_02_26_minigraph_cactus_hgsvc3"
HGSVC3_HPRC_BASE="https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/scratch/2024_02_23_minigraph_cactus_hgsvc3_hprc"

# HGSVC3-only (CHM13 coordinates)
fetch "${RESOURCES}" "${HGSVC3_BASE}/hgsvc3-2024-02-23-mc-chm13-vcfbub.a100k.wave.norm.vcf.gz"
fetch "${RESOURCES}" "${HGSVC3_BASE}/hgsvc3-2024-02-23-mc-chm13-vcfbub.a100k.wave.norm.vcf.gz.tbi"

# HGSVC3-only (GRCh38 coordinates)
fetch "${RESOURCES}" "${HGSVC3_BASE}/hgsvc3-2024-02-23-mc-chm13.GRCh38-vcfbub.a100k.wave.norm.vcf.gz"
fetch "${RESOURCES}" "${HGSVC3_BASE}/hgsvc3-2024-02-23-mc-chm13.GRCh38-vcfbub.a100k.wave.norm.vcf.gz.tbi"

# HGSVC3+HPRC combined (CHM13 coordinates)
fetch "${RESOURCES}" "${HGSVC3_HPRC_BASE}/hgsvc3-hprc-2024-02-23-mc-chm13-vcfbub.a100k.wave.norm.vcf.gz"
fetch "${RESOURCES}" "${HGSVC3_HPRC_BASE}/hgsvc3-hprc-2024-02-23-mc-chm13-vcfbub.a100k.wave.norm.vcf.gz.tbi"

# HGSVC3+HPRC combined (GRCh38 coordinates)
fetch "${RESOURCES}" "${HGSVC3_HPRC_BASE}/hgsvc3-hprc-2024-02-23-mc-chm13.GRCh38-vcfbub.a100k.wave.norm.vcf.gz"
fetch "${RESOURCES}" "${HGSVC3_HPRC_BASE}/hgsvc3-hprc-2024-02-23-mc-chm13.GRCh38-vcfbub.a100k.wave.norm.vcf.gz.tbi"

# ── SGDP ground truth (T2T coordinates) ───────────────────────────────────────
echo ""
echo "==> Downloading SGDP ground truth VCFs (T2T-CHM13 coordinates)..."

SGDP_T2T_BASE="https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/variants/SGDP/chm13v2.0"
SGDP_T2T_DIR="${RESOURCES}/SGDP_variation/t2t"

for chr in "${CHR_ARRAY[@]}"; do
    fetch "${SGDP_T2T_DIR}" "${SGDP_T2T_BASE}/SGDP.CHM13v2.0.chr${chr}.recalibrated.snp_indel.pass.vcf.gz"
    fetch "${SGDP_T2T_DIR}" "${SGDP_T2T_BASE}/SGDP.CHM13v2.0.chr${chr}.recalibrated.snp_indel.pass.vcf.gz.tbi"
done

# ── SGDP ground truth (GRCh38 coordinates) ───────────────────────────────────
echo ""
echo "==> SGDP ground truth (GRCh38 coordinates)"
echo "   The GRCh38-coordinate SGDP VCFs must be downloaded manually."
echo "   Obtain them from Zenodo: [TBD]"
echo "   (Originally sourced from https://anvil.terra.bio/#workspaces/anvil-datastorage/AnVIL_T2T_CHRY/data)"
echo "   Place the per-chromosome VCFs at:"
echo "     ${RESOURCES}/SGDP_variation/grch38/"
echo "   Expected naming: chr{N}.recalibrated.snp_indel.pass.vcf.gz (with .tbi indexes)"

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "==> Download complete."
echo ""
echo "Directory contents:"
echo "  resources/               - Reference genomes, pangenome VCFs, chain files"
echo "  resources/SGDP_variation - SGDP ground truth for imputation validation"
echo "  unphased_variant_calls/  - Raw 1KGP variant calls (input to phasing)"
echo "  phased_panels/grch38/    - GRCh38 phased panels (Byrska-Bishop et al. 2022)"
echo ""
echo "Verify with:"
echo "  ls -lh resources/*.fa.gz resources/*.vcf.gz"
echo "  ls -lh unphased_variant_calls/t2t/*.vcf.gz | head -3"
echo "  ls -lh phased_panels/grch38/*.vcf.gz | head -3"
echo "  ls -lh resources/SGDP_variation/t2t/*.vcf.gz | head -3"
echo ""
echo "Remaining manual steps:"
echo "  1. Download SGDP GRCh38 VCFs (see above)"
echo "  2. Run the pipeline: ./scripts/create_and_assess_haplotype_panels.sh chr22_test 12 my_suffix CHM13v2.0"
