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
#   --test                   Test mode: download only chr15/chr22 test-region data
#                            Uses bcftools to stream only the relevant regions from
#                            remote VCFs, drastically reducing download size.
#   --help                   Show this help
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
RESOURCES="${REPO_ROOT}/resources"

# ── Parse arguments ───────────────────────────────────────────────────────────
CHROMOSOMES=""
SKIP_UNPHASED=false
SKIP_GRCH38_PANELS=false
TEST_MODE=false

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
        --test)
            TEST_MODE=true
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

# Apply --test defaults
if [[ "${TEST_MODE}" == true ]]; then
    [[ -z "${CHROMOSOMES}" ]] && CHROMOSOMES="15,22"
    echo "==> Test mode: downloading only test-region data for chr15/chr22"
    echo "    (using bcftools to stream subsets from remote VCFs)"
fi

# Build chromosome list
if [[ -n "${CHROMOSOMES}" ]]; then
    IFS=',' read -ra CHR_ARRAY <<< "${CHROMOSOMES}"
else
    CHR_ARRAY=()
    for i in $(seq 1 22); do CHR_ARRAY+=("$i"); done
    CHR_ARRAY+=("X")
fi

# ── Test-region coordinates ──────────────────────────────────────────────────
# These must match the regions defined in create_and_assess_haplotype_panels.sh
declare -A T2T_TEST_REGIONS GRCH38_TEST_REGIONS
T2T_TEST_REGIONS=(
    [15]="chr15:17904139-25242443"
    [22]="chr22:18671427-25461594"
)
GRCH38_TEST_REGIONS=(
    [15]="chr15:20000000-27500000"
    [22]="chr22:18000000-25000000"
)

# ── Create directory structure ────────────────────────────────────────────────
mkdir -p \
    "${RESOURCES}" \
    "${RESOURCES}/SGDP_variation/t2t" \
    "${RESOURCES}/SGDP_variation/grch38" \
    "${REPO_ROOT}/unphased_variant_calls/t2t" \
    "${REPO_ROOT}/unphased_variant_calls/grch38" \
    "${REPO_ROOT}/phased_panels/grch38"

# Clean up stale .downloading files from previous interrupted runs
find "${RESOURCES}" "${REPO_ROOT}/unphased_variant_calls" "${REPO_ROOT}/phased_panels/grch38" \
    -name '*.downloading' -type f -delete 2>/dev/null || true

# ── Helpers ──────────────────────────────────────────────────────────────────
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

# Stream a region from a remote VCF using bcftools, saving a bgzipped subset.
# Usage: fetch_region <dest.vcf.gz> <remote_url.vcf.gz> <region>
# Requires the remote server to support range requests and a .tbi index to exist.
fetch_region() {
    local dest_path="$1" url="$2" region="$3"
    if [[ -f "${dest_path}" ]]; then
        echo "    Already exists: $(basename "${dest_path}")"
        return 0
    fi
    local tmpfile="${dest_path}.downloading"
    bcftools view -r "${region}" "${url}" -Oz -o "${tmpfile}" && \
        bcftools index -t "${tmpfile}" && \
        mv "${tmpfile}" "${dest_path}" && \
        mv "${tmpfile}.tbi" "${dest_path}.tbi" || {
            rm -f "${tmpfile}" "${tmpfile}.tbi"
            echo "WARNING: failed to extract region ${region} from ${url}" >&2
        }
}

# Ensure a .fa.gz file is bgzipped (required by samtools faidx).
# The upstream CHM13 .fa.gz is regular gzip; convert to bgzip if needed.
ensure_bgzip() {
    local fasta="$1"
    [[ -f "${fasta}" ]] || return 0
    local flag
    flag=$(od -A n -t x1 -j 3 -N 1 "${fasta}" | tr -d ' ')
    if [[ "${flag}" != "04" ]]; then
        echo "    Re-compressing $(basename "${fasta}") as bgzip..."
        local tmp="${fasta}.rebgzip"
        gunzip -c "${fasta}" | bgzip -@ 2 > "${tmp}" && mv "${tmp}" "${fasta}"
    fi
}

# Collect background PIDs and wait for all to succeed.
BGPIDS=()
bg_wait() {
    local failed=0
    for pid in "${BGPIDS[@]}"; do
        wait "${pid}" || ((failed++))
    done
    BGPIDS=()
    if [[ ${failed} -gt 0 ]]; then
        echo "WARNING: ${failed} background job(s) failed in this section" >&2
    fi
}

# ── Reference genomes ─────────────────────────────────────────────────────────
echo "==> Downloading reference genomes..."

fetch "${RESOURCES}" \
    "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/GCA_009914755.4/chm13v2.0.fa.gz"

# GRCh38 FASTA: download zstd-compressed version, stream-decompress, and bgzip.
if [[ ! -f "${RESOURCES}/GRCh38_full_analysis_set_plus_decoy_hla.fa.gz" ]]; then
    echo "    Downloading, decompressing, and bgzipping GRCh38 reference..."
    wget -qO- "https://www.dropbox.com/s/xyggouv3tnamh0j/GRCh38_full_analysis_set_plus_decoy_hla.fa.zst?dl=1" \
        | zstd -dq \
        | bgzip -@ 2 > "${RESOURCES}/GRCh38_full_analysis_set_plus_decoy_hla.fa.gz"
else
    echo "    GRCh38 FASTA already present, skipping."
fi

ensure_bgzip "${RESOURCES}/chm13v2.0.fa.gz"
ensure_bgzip "${RESOURCES}/GRCh38_full_analysis_set_plus_decoy_hla.fa.gz"

# Index reference FASTAs (creates .fa.gz.fai and .fa.gz.gzi)
for fasta in "${RESOURCES}/chm13v2.0.fa.gz" "${RESOURCES}/GRCh38_full_analysis_set_plus_decoy_hla.fa.gz"; do
    if [[ -f "${fasta}" && ! -f "${fasta}.fai" ]]; then
        echo "    Indexing $(basename "${fasta}") (generates .fai and .gzi)..."
        samtools faidx "${fasta}"
    fi
done

# ── Unphased variant calls ───────────────────────────────────────────────────
if [[ "${SKIP_UNPHASED}" == false ]]; then
    echo ""
    echo "==> Downloading unphased variant calls..."

    T2T_UNPHASED_BASE="https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/variants/1000_Genomes_Project/chm13v2.0/all_samples_3202"
    # bcftools can use s3:// directly; wget needs https://
    GRCH38_UNPHASED_S3="s3://1000genomes/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot"
    GRCH38_UNPHASED_HTTPS="https://s3.amazonaws.com/1000genomes/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot"

    for chr in "${CHR_ARRAY[@]}"; do
        echo "  chr${chr}..."
        if [[ "${TEST_MODE}" == true && -n "${T2T_TEST_REGIONS[${chr}]:-}" ]]; then
            fetch_region "${REPO_ROOT}/unphased_variant_calls/t2t/1KGP.CHM13v2.0.chr${chr}.recalibrated.snp_indel.pass.vcf.gz" \
                "${T2T_UNPHASED_BASE}/1KGP.CHM13v2.0.chr${chr}.recalibrated.snp_indel.pass.vcf.gz" \
                "${T2T_TEST_REGIONS[${chr}]}" &
            BGPIDS+=($!)

            fetch_region "${REPO_ROOT}/unphased_variant_calls/grch38/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.recalibrated_variants.vcf.gz" \
                "${GRCH38_UNPHASED_S3}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.recalibrated_variants.vcf.gz" \
                "${GRCH38_TEST_REGIONS[${chr}]}" &
            BGPIDS+=($!)
        else
            fetch "${REPO_ROOT}/unphased_variant_calls/t2t" \
                "${T2T_UNPHASED_BASE}/1KGP.CHM13v2.0.chr${chr}.recalibrated.snp_indel.pass.vcf.gz" &
            BGPIDS+=($!)
            fetch "${REPO_ROOT}/unphased_variant_calls/t2t" \
                "${T2T_UNPHASED_BASE}/1KGP.CHM13v2.0.chr${chr}.recalibrated.snp_indel.pass.vcf.gz.tbi" &
            BGPIDS+=($!)

            fetch "${REPO_ROOT}/unphased_variant_calls/grch38" \
                "${GRCH38_UNPHASED_HTTPS}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.recalibrated_variants.vcf.gz" &
            BGPIDS+=($!)
            fetch "${REPO_ROOT}/unphased_variant_calls/grch38" \
                "${GRCH38_UNPHASED_HTTPS}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.recalibrated_variants.vcf.gz.tbi" &
            BGPIDS+=($!)
        fi
    done
    bg_wait
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

        if [[ "${TEST_MODE}" == true && -n "${GRCH38_TEST_REGIONS[${chr}]:-}" ]]; then
            echo "  chr${chr}: streaming test region -> ${local_name}"
            fetch_region "${GRCH38_PANEL_DIR}/${local_name}" \
                "${GRCH38_PANEL_BASE}/${upstream_name}" \
                "${GRCH38_TEST_REGIONS[${chr}]}" &
            BGPIDS+=($!)
        else
            echo "  chr${chr}: ${upstream_name} -> ${local_name}"
            fetch_as "${GRCH38_PANEL_DIR}/${local_name}" \
                "${GRCH38_PANEL_BASE}/${upstream_name}" &
            BGPIDS+=($!)
            fetch_as "${GRCH38_PANEL_DIR}/${local_name}.tbi" \
                "${GRCH38_PANEL_BASE}/${upstream_name}.tbi" &
            BGPIDS+=($!)
        fi
    done
    bg_wait

    # Create biallelic 2504-sample BCF panels for imputation comparison.
    # The imputation script expects these as the GRCh38 native reference baseline.
    UNRELATED_SAMPLES="${REPO_ROOT}/resources/sample_subsets/unrelated_samples.txt"
    echo "  Creating biallelic 2504-member BCF panels..."
    for chr in "${CHR_ARRAY[@]}"; do
        local_name="1KGP.GRCh38.chr${chr}.recalibrated.snp_indel.pass.phased"
        source_vcf="${GRCH38_PANEL_DIR}/${local_name}.3202.vcf.gz"
        dest_bcf="${GRCH38_PANEL_DIR}/${local_name}.biallelic.2504.bcf"
        if [[ -f "${dest_bcf}" && -f "${dest_bcf}.csi" ]]; then
            echo "    Already exists: $(basename "${dest_bcf}")"
            continue
        fi
        if [[ ! -f "${source_vcf}" ]]; then
            echo "    WARNING: source panel missing, skipping: $(basename "${source_vcf}")" >&2
            continue
        fi
        echo "    chr${chr}: subsetting to 2504 unrelated, biallelic -> $(basename "${dest_bcf}")"
        bcftools view -Ou --threads 4 -S "${UNRELATED_SAMPLES}" --force-samples "${source_vcf}" \
        | bcftools view -Ou --threads 2 -m2 -M2 -c 1:minor - \
        | bcftools annotate -Ou -x INFO/MAC,INFO/AN,INFO/AC,INFO/MAF - \
        | bcftools +fill-tags -Ou --threads 4 - -- -t AN,AC,MAF,MAC:1=MAC \
        | bcftools annotate -Ou --threads 2 -x ^INFO/MAF,^INFO/MAC,^INFO/AN,^FORMAT/GT - \
        | bcftools view --threads 4 -Ob - > "${dest_bcf}" \
        && bcftools index --threads 4 "${dest_bcf}" &
        BGPIDS+=($!)
    done
    bg_wait
else
    echo ""
    echo "==> Skipping GRCh38 phased panels (--skip-grch38-panels)"
fi

# ── Pangenome VCFs ───────────────────────────────────────────────────────────
echo ""
echo "==> Downloading pangenome VCFs (HPRC v1.1 + HGSVC3)..."

HPRC_BASE="https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus"
HGSVC3_BASE="https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/scratch/2024_02_26_minigraph_cactus_hgsvc3"
HGSVC3_HPRC_BASE="https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/scratch/2024_02_23_minigraph_cactus_hgsvc3_hprc"

if [[ "${TEST_MODE}" == true ]]; then
    # Build combined region strings for all test chromosomes
    t2t_regions="" grch38_regions=""
    for chr in "${CHR_ARRAY[@]}"; do
        [[ -n "${T2T_TEST_REGIONS[${chr}]:-}" ]] && t2t_regions="${t2t_regions:+${t2t_regions},}${T2T_TEST_REGIONS[${chr}]}"
        [[ -n "${GRCH38_TEST_REGIONS[${chr}]:-}" ]] && grch38_regions="${grch38_regions:+${grch38_regions},}${GRCH38_TEST_REGIONS[${chr}]}"
    done

    # All 6 pangenome VCFs in parallel
    fetch_region "${RESOURCES}/hprc-v1.1-mc-chm13.vcfbub.a100k.wave.vcf.gz" \
        "${HPRC_BASE}/hprc-v1.1-mc-chm13/hprc-v1.1-mc-chm13.vcfbub.a100k.wave.vcf.gz" \
        "${t2t_regions}" &
    BGPIDS+=($!)
    fetch_region "${RESOURCES}/hprc-v1.1-mc-grch38.vcfbub.a100k.wave.vcf.gz" \
        "${HPRC_BASE}/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.vcfbub.a100k.wave.vcf.gz" \
        "${grch38_regions}" &
    BGPIDS+=($!)
    fetch_region "${RESOURCES}/hgsvc3-2024-02-23-mc-chm13-vcfbub.a100k.wave.norm.vcf.gz" \
        "${HGSVC3_BASE}/hgsvc3-2024-02-23-mc-chm13-vcfbub.a100k.wave.norm.vcf.gz" \
        "${t2t_regions}" &
    BGPIDS+=($!)
    fetch_region "${RESOURCES}/hgsvc3-2024-02-23-mc-chm13.GRCh38-vcfbub.a100k.wave.norm.vcf.gz" \
        "${HGSVC3_BASE}/hgsvc3-2024-02-23-mc-chm13.GRCh38-vcfbub.a100k.wave.norm.vcf.gz" \
        "${grch38_regions}" &
    BGPIDS+=($!)
    fetch_region "${RESOURCES}/hgsvc3-hprc-2024-02-23-mc-chm13-vcfbub.a100k.wave.norm.vcf.gz" \
        "${HGSVC3_HPRC_BASE}/hgsvc3-hprc-2024-02-23-mc-chm13-vcfbub.a100k.wave.norm.vcf.gz" \
        "${t2t_regions}" &
    BGPIDS+=($!)
    fetch_region "${RESOURCES}/hgsvc3-hprc-2024-02-23-mc-chm13.GRCh38-vcfbub.a100k.wave.norm.vcf.gz" \
        "${HGSVC3_HPRC_BASE}/hgsvc3-hprc-2024-02-23-mc-chm13.GRCh38-vcfbub.a100k.wave.norm.vcf.gz" \
        "${grch38_regions}" &
    BGPIDS+=($!)
    bg_wait
else
    # HPRC v1.1
    fetch "${RESOURCES}" "${HPRC_BASE}/hprc-v1.1-mc-chm13/hprc-v1.1-mc-chm13.vcfbub.a100k.wave.vcf.gz" &
    BGPIDS+=($!)
    fetch "${RESOURCES}" "${HPRC_BASE}/hprc-v1.1-mc-chm13/hprc-v1.1-mc-chm13.vcfbub.a100k.wave.vcf.gz.tbi" &
    BGPIDS+=($!)
    fetch "${RESOURCES}" "${HPRC_BASE}/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.vcfbub.a100k.wave.vcf.gz" &
    BGPIDS+=($!)
    fetch "${RESOURCES}" "${HPRC_BASE}/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.vcfbub.a100k.wave.vcf.gz.tbi" &
    BGPIDS+=($!)

    # HGSVC3-only
    fetch "${RESOURCES}" "${HGSVC3_BASE}/hgsvc3-2024-02-23-mc-chm13-vcfbub.a100k.wave.norm.vcf.gz" &
    BGPIDS+=($!)
    fetch "${RESOURCES}" "${HGSVC3_BASE}/hgsvc3-2024-02-23-mc-chm13-vcfbub.a100k.wave.norm.vcf.gz.tbi" &
    BGPIDS+=($!)
    fetch "${RESOURCES}" "${HGSVC3_BASE}/hgsvc3-2024-02-23-mc-chm13.GRCh38-vcfbub.a100k.wave.norm.vcf.gz" &
    BGPIDS+=($!)
    fetch "${RESOURCES}" "${HGSVC3_BASE}/hgsvc3-2024-02-23-mc-chm13.GRCh38-vcfbub.a100k.wave.norm.vcf.gz.tbi" &
    BGPIDS+=($!)

    # HGSVC3+HPRC combined
    fetch "${RESOURCES}" "${HGSVC3_HPRC_BASE}/hgsvc3-hprc-2024-02-23-mc-chm13-vcfbub.a100k.wave.norm.vcf.gz" &
    BGPIDS+=($!)
    fetch "${RESOURCES}" "${HGSVC3_HPRC_BASE}/hgsvc3-hprc-2024-02-23-mc-chm13-vcfbub.a100k.wave.norm.vcf.gz.tbi" &
    BGPIDS+=($!)
    fetch "${RESOURCES}" "${HGSVC3_HPRC_BASE}/hgsvc3-hprc-2024-02-23-mc-chm13.GRCh38-vcfbub.a100k.wave.norm.vcf.gz" &
    BGPIDS+=($!)
    fetch "${RESOURCES}" "${HGSVC3_HPRC_BASE}/hgsvc3-hprc-2024-02-23-mc-chm13.GRCh38-vcfbub.a100k.wave.norm.vcf.gz.tbi" &
    BGPIDS+=($!)
    bg_wait
fi

# ── SGDP ground truth (T2T coordinates) ───────────────────────────────────────
echo ""
echo "==> Downloading SGDP ground truth VCFs (T2T-CHM13 coordinates)..."

SGDP_T2T_BASE="https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/variants/SGDP/chm13v2.0"
SGDP_T2T_DIR="${RESOURCES}/SGDP_variation/t2t"

for chr in "${CHR_ARRAY[@]}"; do
    if [[ "${TEST_MODE}" == true && -n "${T2T_TEST_REGIONS[${chr}]:-}" ]]; then
        fetch_region "${SGDP_T2T_DIR}/SGDP.CHM13v2.0.chr${chr}.recalibrated.snp_indel.pass.vcf.gz" \
            "${SGDP_T2T_BASE}/SGDP.CHM13v2.0.chr${chr}.recalibrated.snp_indel.pass.vcf.gz" \
            "${T2T_TEST_REGIONS[${chr}]}" &
        BGPIDS+=($!)
    else
        fetch "${SGDP_T2T_DIR}" "${SGDP_T2T_BASE}/SGDP.CHM13v2.0.chr${chr}.recalibrated.snp_indel.pass.vcf.gz" &
        BGPIDS+=($!)
        fetch "${SGDP_T2T_DIR}" "${SGDP_T2T_BASE}/SGDP.CHM13v2.0.chr${chr}.recalibrated.snp_indel.pass.vcf.gz.tbi" &
        BGPIDS+=($!)
    fi
done
bg_wait

# ── SGDP ground truth (GRCh38 coordinates) ───────────────────────────────────
echo ""
echo "==> Downloading SGDP ground truth VCFs (GRCh38 coordinates)..."

SGDP_GRCH38_DIR="${RESOURCES}/SGDP_variation/grch38"
SGDP_GRCH38_ZENODO="https://zenodo.org/records/19371182/files"

if [[ "${TEST_MODE}" == true ]]; then
    # Test mode: download pre-subsetted test-region VCFs from Zenodo
    SGDP_GRCH38_TARBALL="${SGDP_GRCH38_DIR}/SGDP_GRCh38_test_regions.tar.gz"
    if [[ -f "${SGDP_GRCH38_DIR}/chr22.recalibrated.snp_indel.pass.vcf.gz" && \
          -f "${SGDP_GRCH38_DIR}/chr15.recalibrated.snp_indel.pass.vcf.gz" ]]; then
        echo "    SGDP GRCh38 test VCFs already present, skipping."
    else
        echo "    Downloading SGDP GRCh38 test-region tarball from Zenodo..."
        wget -q -O "${SGDP_GRCH38_TARBALL}" \
            "${SGDP_GRCH38_ZENODO}/SGDP_GRCh38_test_regions.tar.gz?download=1" && \
            tar -xzf "${SGDP_GRCH38_TARBALL}" -C "${SGDP_GRCH38_DIR}" && \
            rm -f "${SGDP_GRCH38_TARBALL}" || {
                rm -f "${SGDP_GRCH38_TARBALL}"
                echo "WARNING: failed to download SGDP GRCh38 test data from Zenodo" >&2
            }
    fi
    # Index any extracted VCFs that are missing a .tbi
    for vcf in "${SGDP_GRCH38_DIR}"/*.vcf.gz; do
        [[ -f "${vcf}" ]] || continue
        if [[ ! -f "${vcf}.tbi" ]]; then
            echo "    Indexing $(basename "${vcf}")..."
            bcftools index -t "${vcf}"
        fi
    done
else
    # Full mode: download whole-genome SGDP GRCh38 tarball from Zenodo
    SGDP_GRCH38_TARBALL="${SGDP_GRCH38_DIR}/SGDP_GRCh38_all_chromosomes.tar.gz"
    # Check if all requested chromosomes are already present
    all_present=true
    for chr in "${CHR_ARRAY[@]}"; do
        if [[ ! -f "${SGDP_GRCH38_DIR}/chr${chr}.recalibrated.snp_indel.pass.vcf.gz" ]]; then
            all_present=false
            break
        fi
    done
    if [[ "${all_present}" == true ]]; then
        echo "    SGDP GRCh38 VCFs already present, skipping."
    else
        echo "    Downloading SGDP GRCh38 whole-genome tarball from Zenodo..."
        wget -q -O "${SGDP_GRCH38_TARBALL}" \
            "${SGDP_GRCH38_ZENODO}/SGDP_GRCh38_all_chromosomes.tar.gz?download=1" && \
            tar -xzf "${SGDP_GRCH38_TARBALL}" -C "${SGDP_GRCH38_DIR}" && \
            rm -f "${SGDP_GRCH38_TARBALL}" || {
                rm -f "${SGDP_GRCH38_TARBALL}"
                echo "WARNING: failed to download SGDP GRCh38 data from Zenodo" >&2
            }
    fi
    # Index any extracted VCFs that are missing a .tbi
    for vcf in "${SGDP_GRCH38_DIR}"/*.vcf.gz; do
        [[ -f "${vcf}" ]] || continue
        if [[ ! -f "${vcf}.tbi" ]]; then
            echo "    Indexing $(basename "${vcf}")..."
            bcftools index -t "${vcf}"
        fi
    done
fi

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
echo "Next steps:"
echo "  Run the pipeline: ./scripts/create_and_assess_haplotype_panels.sh chr22_test 12 my_suffix CHM13v2.0"
echo ""
echo "For a quick test with minimal downloads:"
echo "  bash scripts/utility/download_resources.sh --test"
