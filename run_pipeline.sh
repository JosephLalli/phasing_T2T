#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: ./run_pipeline.sh CHROM THREADS RUN_SUFFIX GENOME

Run scripts/create_and_assess_haplotype_panels.sh in the main Docker image.
GENOME must be CHM13v2.0 or GRCh38. Set PHASING_T2T_IMAGE to override the
default image tag (jlalli/phasing_t2t:v2.1).
EOF
}

if [[ $# -eq 1 && ( $1 == -h || $1 == --help ) ]]; then
    usage
    exit 0
fi

if [[ $# -ne 4 ]]; then
    usage >&2
    exit 2
fi

chrom=$1
threads=$2
run_suffix=$3
genome=$4

if [[ ! $threads =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: THREADS must be a positive integer." >&2
    exit 2
fi

if [[ $genome != CHM13v2.0 && $genome != GRCh38 ]]; then
    echo "ERROR: GENOME must be CHM13v2.0 or GRCh38." >&2
    exit 2
fi

project_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd -P)"
docker_args=(
    run --rm
    --user "$(id -u):$(id -g)"
    --cpus "$threads"
    --mount "type=bind,src=$project_dir,dst=/phasing_T2T_project"
    --workdir /phasing_T2T_project
)

if [[ -v PHASING_T2T_RUN_PROFILE ]]; then
    docker_args+=(--env "PHASING_T2T_RUN_PROFILE=$PHASING_T2T_RUN_PROFILE")
fi

exec docker "${docker_args[@]}" "${PHASING_T2T_IMAGE:-jlalli/phasing_t2t:v2.1}" \
    bash scripts/create_and_assess_haplotype_panels.sh \
    "$chrom" "$threads" "$run_suffix" "$genome"
