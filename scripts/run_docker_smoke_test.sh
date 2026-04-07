#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

usage() {
    cat <<'EOF'
Usage:
  scripts/run_docker_smoke_test.sh [docker.env]

Description:
  Run the chr22_test and chr15_test smoke workflow for both CHM13v2.0
  and GRCh38 inside Docker, then aggregate the test outputs and execute
  the downstream analysis notebooks and Figure 6 generation.
  This validates the containerized software stack and repo entrypoints,
  but still requires externally provided biological inputs.

Environment:
  Copy docker.env.example to docker.env, then review the settings.
  If you used scripts/utility/download_resources.sh, the canonical repo-local
  paths in docker.env.example should work without edits.
EOF
}

ENV_FILE="${1:-docker.env}"
if [[ "${ENV_FILE}" == "-h" || "${ENV_FILE}" == "--help" ]]; then
    usage
    exit 0
fi

if [[ ! -f "${ENV_FILE}" ]]; then
    echo "ERROR: environment file not found: ${ENV_FILE}" >&2
    echo "Create it from docker.env.example." >&2
    echo "For the canonical test-data setup, run:" >&2
    echo "  bash scripts/utility/download_resources.sh --test" >&2
    exit 1
fi

# shellcheck disable=SC1090
source "${ENV_FILE}"

IMAGE_NAME="${IMAGE_NAME:-jlalli/phasing_t2t_dep_container:v2.0}"
CONTAINER_NAME="${CONTAINER_NAME:-phasing-t2t-smoke}"
NUM_THREADS="${NUM_THREADS:-12}"
RUN_SUFFIX="${RUN_SUFFIX:-smoke}"
RUN_NOTEBOOKS="${RUN_NOTEBOOKS:-1}"
OUTPUT_DIR="${OUTPUT_DIR:-${PROJECT_DIR}/docker_smoke_output}"

require_var() {
    local name="$1"
    if [[ -z "${!name:-}" ]]; then
        echo "ERROR: required variable ${name} is not set in ${ENV_FILE}" >&2
        exit 1
    fi
}

require_file() {
    local name="$1"
    local path="$2"
    if [[ ! -f "${path}" ]]; then
        echo "ERROR: ${name} does not exist: ${path}" >&2
        exit 1
    fi
}

require_dir() {
    local name="$1"
    local path="$2"
    if [[ ! -d "${path}" ]]; then
        echo "ERROR: ${name} does not exist: ${path}" >&2
        exit 1
    fi
}

require_var CHM13_FA
require_var GRCH38_FA
require_var UNPHASED_T2T_DIR
require_var UNPHASED_GRCH38_DIR
require_var GRCH38_PANELS_DIR
require_var SGDP_GRCH38_DIR
require_var SGDP_T2T_DIR
require_var HPRC_CHM13_VCF
require_var HPRC_GRCH38_VCF
require_var HGSVC3_CHM13_VCF
require_var HGSVC3_CHM13_GRCH38_VCF
require_var HGSVC3_HPRC_CHM13_VCF
require_var HGSVC3_HPRC_CHM13_GRCH38_VCF

require_file CHM13_FA "${CHM13_FA}"
require_file CHM13_FA_FAI "${CHM13_FA}.fai"
require_file CHM13_FA_GZI "${CHM13_FA}.gzi"
require_file GRCH38_FA "${GRCH38_FA}"
require_file GRCH38_FA_FAI "${GRCH38_FA}.fai"
require_file GRCH38_FA_GZI "${GRCH38_FA}.gzi"

require_dir UNPHASED_T2T_DIR "${UNPHASED_T2T_DIR}"
require_dir UNPHASED_GRCH38_DIR "${UNPHASED_GRCH38_DIR}"
require_dir GRCH38_PANELS_DIR "${GRCH38_PANELS_DIR}"
require_dir SGDP_GRCH38_DIR "${SGDP_GRCH38_DIR}"
require_dir SGDP_T2T_DIR "${SGDP_T2T_DIR}"

require_file HPRC_CHM13_VCF "${HPRC_CHM13_VCF}"
require_file HPRC_CHM13_TBI "${HPRC_CHM13_VCF}.tbi"
require_file HPRC_GRCH38_VCF "${HPRC_GRCH38_VCF}"
require_file HPRC_GRCH38_TBI "${HPRC_GRCH38_VCF}.tbi"
require_file HGSVC3_CHM13_VCF "${HGSVC3_CHM13_VCF}"
require_file HGSVC3_CHM13_TBI "${HGSVC3_CHM13_VCF}.tbi"
require_file HGSVC3_CHM13_GRCH38_VCF "${HGSVC3_CHM13_GRCH38_VCF}"
require_file HGSVC3_CHM13_GRCH38_TBI "${HGSVC3_CHM13_GRCH38_VCF}.tbi"
require_file HGSVC3_HPRC_CHM13_VCF "${HGSVC3_HPRC_CHM13_VCF}"
require_file HGSVC3_HPRC_CHM13_TBI "${HGSVC3_HPRC_CHM13_VCF}.tbi"
require_file HGSVC3_HPRC_CHM13_GRCH38_VCF "${HGSVC3_HPRC_CHM13_GRCH38_VCF}"
require_file HGSVC3_HPRC_CHM13_GRCH38_TBI "${HGSVC3_HPRC_CHM13_GRCH38_VCF}.tbi"

mkdir -p \
    "${OUTPUT_DIR}/working_directories" \
    "${OUTPUT_DIR}/intermediate_data" \
    "${OUTPUT_DIR}/imputation_statistics" \
    "${OUTPUT_DIR}/SHAPEIT5_switch_output" \
    "${OUTPUT_DIR}/phased_panels" \
    "${OUTPUT_DIR}/figures" \
    "${OUTPUT_DIR}/tables" \
    "${OUTPUT_DIR}/notebook_runs"

docker run --rm \
    --name "${CONTAINER_NAME}" \
    -e NUM_THREADS="${NUM_THREADS}" \
    -e RUN_SUFFIX="${RUN_SUFFIX}" \
    -e RUN_NOTEBOOKS="${RUN_NOTEBOOKS}" \
    -v "${CHM13_FA}:/phasing_T2T_project/resources/chm13v2.0.fa.gz:ro" \
    -v "${CHM13_FA}.fai:/phasing_T2T_project/resources/chm13v2.0.fa.gz.fai:ro" \
    -v "${CHM13_FA}.gzi:/phasing_T2T_project/resources/chm13v2.0.fa.gz.gzi:ro" \
    -v "${GRCH38_FA}:/phasing_T2T_project/resources/GRCh38_full_analysis_set_plus_decoy_hla.fa.gz:ro" \
    -v "${GRCH38_FA}.fai:/phasing_T2T_project/resources/GRCh38_full_analysis_set_plus_decoy_hla.fa.gz.fai:ro" \
    -v "${GRCH38_FA}.gzi:/phasing_T2T_project/resources/GRCh38_full_analysis_set_plus_decoy_hla.fa.gz.gzi:ro" \
    -v "${UNPHASED_T2T_DIR}:/phasing_T2T_project/unphased_variant_calls/t2t:ro" \
    -v "${UNPHASED_GRCH38_DIR}:/phasing_T2T_project/unphased_variant_calls/grch38:ro" \
    -v "${HPRC_CHM13_VCF}:/phasing_T2T_project/resources/hprc-v1.1-mc-chm13.vcfbub.a100k.wave.vcf.gz:ro" \
    -v "${HPRC_CHM13_VCF}.tbi:/phasing_T2T_project/resources/hprc-v1.1-mc-chm13.vcfbub.a100k.wave.vcf.gz.tbi:ro" \
    -v "${HPRC_GRCH38_VCF}:/phasing_T2T_project/resources/hprc-v1.1-mc-grch38.vcfbub.a100k.wave.vcf.gz:ro" \
    -v "${HPRC_GRCH38_VCF}.tbi:/phasing_T2T_project/resources/hprc-v1.1-mc-grch38.vcfbub.a100k.wave.vcf.gz.tbi:ro" \
    -v "${HGSVC3_CHM13_VCF}:/phasing_T2T_project/resources/hgsvc3-2024-02-23-mc-chm13-vcfbub.a100k.wave.norm.vcf.gz:ro" \
    -v "${HGSVC3_CHM13_VCF}.tbi:/phasing_T2T_project/resources/hgsvc3-2024-02-23-mc-chm13-vcfbub.a100k.wave.norm.vcf.gz.tbi:ro" \
    -v "${HGSVC3_CHM13_GRCH38_VCF}:/phasing_T2T_project/resources/hgsvc3-2024-02-23-mc-chm13.GRCh38-vcfbub.a100k.wave.norm.vcf.gz:ro" \
    -v "${HGSVC3_CHM13_GRCH38_VCF}.tbi:/phasing_T2T_project/resources/hgsvc3-2024-02-23-mc-chm13.GRCh38-vcfbub.a100k.wave.norm.vcf.gz.tbi:ro" \
    -v "${HGSVC3_HPRC_CHM13_VCF}:/phasing_T2T_project/resources/hgsvc3-hprc-2024-02-23-mc-chm13-vcfbub.a100k.wave.norm.vcf.gz:ro" \
    -v "${HGSVC3_HPRC_CHM13_VCF}.tbi:/phasing_T2T_project/resources/hgsvc3-hprc-2024-02-23-mc-chm13-vcfbub.a100k.wave.norm.vcf.gz.tbi:ro" \
    -v "${HGSVC3_HPRC_CHM13_GRCH38_VCF}:/phasing_T2T_project/resources/hgsvc3-hprc-2024-02-23-mc-chm13.GRCh38-vcfbub.a100k.wave.norm.vcf.gz:ro" \
    -v "${HGSVC3_HPRC_CHM13_GRCH38_VCF}.tbi:/phasing_T2T_project/resources/hgsvc3-hprc-2024-02-23-mc-chm13.GRCh38-vcfbub.a100k.wave.norm.vcf.gz.tbi:ro" \
    -v "${SGDP_GRCH38_DIR}:/phasing_T2T_project/resources/SGDP_variation/grch38:ro" \
    -v "${SGDP_T2T_DIR}:/phasing_T2T_project/resources/SGDP_variation/t2t:ro" \
    -v "${OUTPUT_DIR}/working_directories:/phasing_T2T_project/working_directories" \
    -v "${OUTPUT_DIR}/intermediate_data:/phasing_T2T_project/intermediate_data" \
    -v "${OUTPUT_DIR}/imputation_statistics:/phasing_T2T_project/imputation_statistics" \
    -v "${OUTPUT_DIR}/SHAPEIT5_switch_output:/phasing_T2T_project/SHAPEIT5_switch_output" \
    -v "${OUTPUT_DIR}/phased_panels:/phasing_T2T_project/phased_panels" \
    -v "${GRCH38_PANELS_DIR}:/phasing_T2T_project/phased_panels/grch38" \
    -v "${OUTPUT_DIR}/figures:/phasing_T2T_project/figures" \
    -v "${OUTPUT_DIR}/tables:/phasing_T2T_project/tables" \
    -v "${OUTPUT_DIR}/notebook_runs:/phasing_T2T_project/notebook_runs" \
    "${IMAGE_NAME}" \
    bash -lc '
        set -euo pipefail

        make_run_suffix() {
            local base_suffix="$1"
            local genome="$2"
            if [[ -n "${base_suffix}" ]]; then
                printf "%s_%s\n" "${base_suffix}" "${genome}"
            else
                printf "%s\n" "${genome}"
            fi
        }

        CHM13_RUN_SUFFIX="$(make_run_suffix "${RUN_SUFFIX}" "CHM13v2.0")"
        GRCH38_RUN_SUFFIX="$(make_run_suffix "${RUN_SUFFIX}" "GRCh38")"

        run_panel_test() {
            local chrom="$1"
            local genome="$2"
            local run_suffix="$3"
            echo "=== Stage 1: ${genome} ${chrom} ==="
            ./scripts/create_and_assess_haplotype_panels.sh "${chrom}" "${NUM_THREADS}" "${run_suffix}" "${genome}"
        }

        run_panel_test chr22_test CHM13v2.0 "${CHM13_RUN_SUFFIX}" &
        run_panel_test chr15_test CHM13v2.0 "${CHM13_RUN_SUFFIX}" &
        run_panel_test chr22_test GRCh38 "${GRCH38_RUN_SUFFIX}" &
        run_panel_test chr15_test GRCh38 "${GRCH38_RUN_SUFFIX}" &

        wait

        IMPUTATION_RESULTS_DIR="./imputation_statistics/imputation_results_${CHM13_RUN_SUFFIX}"

        echo "=== Stage 2: syntenic/nonsyntenic bin generation ==="
        ./scripts/create_syn_nonsyn_bins.sh CHM13v2.0 "${CHM13_RUN_SUFFIX}" "${IMPUTATION_RESULTS_DIR}" true
        ./scripts/create_syn_nonsyn_bins.sh GRCh38 "${GRCH38_RUN_SUFFIX}" "${IMPUTATION_RESULTS_DIR}" true

        echo "=== Stage 3: genome-wide imputation aggregation ==="
        ./scripts/calc_genomewide_imputation_statistics_full.sh \
            "${GRCH38_RUN_SUFFIX}" \
            "${CHM13_RUN_SUFFIX}" \
            "${NUM_THREADS}" \
            true

        echo "=== Stage 4: summary parquet generation ==="
        python3 ./scripts/analysis/create_summary_phasing_dataframes_polars_regional.py \
            --CHM13_run_suffix "${CHM13_RUN_SUFFIX}" \
            --GRCh38_run_suffix "${GRCH38_RUN_SUFFIX}" \
            --test

        if [[ "${RUN_NOTEBOOKS}" == "1" ]]; then
            echo "=== Stage 5: notebook execution ==="
            (
                cd notebooks
                for nb in calc_figures_for_paper.ipynb calc_per_variant_figures_for_paper.ipynb make_plots.ipynb; do
                    echo "Executing ${nb}"
                    jupyter nbconvert \
                        --to notebook \
                        --execute \
                        --ExecutePreprocessor.timeout=600 \
                        --output-dir /phasing_T2T_project/notebook_runs \
                        "${nb}"
                done
            )

            echo "=== Stage 6: Figure 6 generation ==="
            Rscript ./scripts/figure6/Figure_6_script.R
            python3 ./scripts/figure6/stitch_svgs.py --batch ./figures/figure6
            python3 ./scripts/figure6/stitch_svgs.py --grid ./figures/figure6
        fi

        echo "=== Smoke test complete ==="
        ls -1 /phasing_T2T_project/intermediate_data/*.parquet 2>/dev/null || true
    '
