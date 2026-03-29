#!/usr/bin/env bash
# Download large external resources required by the phasing pipeline.
# Run from the repository root:  bash scripts/utility/download_resources.sh
#
# Files that already exist are skipped (wget -nc).  Re-run safely at any time.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
RESOURCES="${REPO_ROOT}/resources"

mkdir -p "${RESOURCES}" "${RESOURCES}/SGDP_variation"

# ── Helper ────────────────────────────────────────────────────────────────────
fetch() {
    local dest_dir="$1" url="$2"
    wget -nc -P "${dest_dir}" "${url}" || {
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
fetch "${RESOURCES}" \
    "https://42basepairs.com/download/s3/1000genomes/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai"

# ── HPRC v1.1 pangenome VCFs ─────────────────────────────────────────────────
echo "==> Downloading HPRC v1.1 pangenome VCFs..."

HPRC_BASE="https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus"

fetch "${RESOURCES}" "${HPRC_BASE}/hprc-v1.1-mc-chm13/hprc-v1.1-mc-chm13.vcfbub.a100k.wave.vcf.gz"
fetch "${RESOURCES}" "${HPRC_BASE}/hprc-v1.1-mc-chm13/hprc-v1.1-mc-chm13.vcfbub.a100k.wave.vcf.gz.tbi"

fetch "${RESOURCES}" "${HPRC_BASE}/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.vcfbub.a100k.wave.vcf.gz"
fetch "${RESOURCES}" "${HPRC_BASE}/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.vcfbub.a100k.wave.vcf.gz.tbi"

# ── HGSVC3 pangenome VCFs ────────────────────────────────────────────────────
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
echo "==> Downloading SGDP ground truth VCFs (T2T-CHM13 coordinates)..."

SGDP_T2T_BASE="https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/variants/SGDP/chm13v2.0"
SGDP_T2T_DIR="${RESOURCES}/SGDP_variation/t2t"
mkdir -p "${SGDP_T2T_DIR}"

for chr in {1..22} X; do
    fetch "${SGDP_T2T_DIR}" "${SGDP_T2T_BASE}/SGDP.CHM13v2.0.chr${chr}.recalibrated.snp_indel.pass.vcf.gz"
    fetch "${SGDP_T2T_DIR}" "${SGDP_T2T_BASE}/SGDP.CHM13v2.0.chr${chr}.recalibrated.snp_indel.pass.vcf.gz.tbi"
done

# ── SGDP ground truth (GRCh38 coordinates) ───────────────────────────────────
echo ""
echo "==> SGDP ground truth (GRCh38 coordinates)"
echo "   The GRCh38-coordinate SGDP VCFs must be downloaded manually."
echo "   Obtain them from Zenodo: [TBD]"
echo "   (Originally sourced from https://anvil.terra.bio/#workspaces/anvil-datastorage/AnVIL_T2T_CHRY/data)"
echo "   Place (or symlink) the per-chromosome VCFs at:"
echo "     resources/SGDP_variation/grch38/"
echo "   Expected naming: chr{N}.recalibrated.snp_indel.pass.vcf.gz"
echo "   See resources/SGDP_variation/README.md for details."

echo ""
echo "==> Done.  Verify with: ls -lh resources/*.vcf.gz resources/*.fa.gz resources/SGDP_variation/t2t/*.vcf.gz"
