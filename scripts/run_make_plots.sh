#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
INPUT_DIR=
OUTPUT_DIR=

usage() {
    echo 'Usage: scripts/run_make_plots.sh --inputs INPUTDIR --output OUTPUTDIR'
    echo 'Run make_plots.ipynb in the main Docker image with the supplied figure inputs.'
    echo 'Set PHASING_T2T_IMAGE to use another tag (default: jlalli/phasing_t2t:v2.1).'
}

while (($#)); do
    case "$1" in
        --inputs) INPUT_DIR=${2:?missing INPUTDIR}; shift 2 ;;
        --output) OUTPUT_DIR=${2:?missing OUTPUTDIR}; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done
[[ -n "$INPUT_DIR" && -n "$OUTPUT_DIR" ]] || { usage >&2; exit 2; }
INPUT_DIR="$(cd -- "$INPUT_DIR" && pwd -P)"
for directory in intermediate_data_whole_genome imputation_statistics_whole_genome; do
    [[ -d "$INPUT_DIR/$directory" ]] || { echo "ERROR: missing input directory: $directory" >&2; exit 2; }
done
(cd -- "$INPUT_DIR" && sha256sum --check --quiet SHA256SUMS.txt)

if [[ -e "$OUTPUT_DIR" ]] && [[ -n "$(ls -A -- "$OUTPUT_DIR")" ]]; then
    echo "ERROR: output directory must be empty: $OUTPUT_DIR" >&2
    exit 2
fi
mkdir -p -- "$OUTPUT_DIR"
OUTPUT_DIR="$(cd -- "$OUTPUT_DIR" && pwd -P)"
for directory in intermediate_data_whole_genome imputation_statistics_whole_genome figures_whole_genome tables_whole_genome; do
    [[ ! -L "$PROJECT_DIR/$directory" ]] || { echo "ERROR: mount destination is a symlink: $directory" >&2; exit 2; }
    mkdir -p -- "$PROJECT_DIR/$directory"
done
mkdir -p -- "$OUTPUT_DIR/figures_whole_genome" "$OUTPUT_DIR/tables_whole_genome" "$OUTPUT_DIR/notebook_runs"

exec docker run --rm --network none --read-only --tmpfs /tmp:rw,size=8g \
    --user "$(id -u):$(id -g)" --cpus "${PHASING_T2T_THREADS:-8}" --memory "${PHASING_T2T_MEMORY:-64g}" \
    -e MPLBACKEND=Agg -e MPLCONFIGDIR=/tmp/matplotlib -e XDG_CACHE_HOME=/tmp/cache \
    -e JUPYTER_DATA_DIR=/tmp/jupyter-data -e JUPYTER_CONFIG_DIR=/tmp/jupyter-config \
    -e JUPYTER_RUNTIME_DIR=/tmp/jupyter-runtime -e IPYTHONDIR=/tmp/ipython \
    -e PHASING_T2T_RUN_PROFILE=whole_genome -e NUM_THREADS="${PHASING_T2T_THREADS:-8}" \
    -e OPENBLAS_NUM_THREADS="${PHASING_T2T_THREADS:-8}" -e MKL_NUM_THREADS="${PHASING_T2T_THREADS:-8}" \
    -e OMP_NUM_THREADS="${PHASING_T2T_THREADS:-8}" -e POLARS_MAX_THREADS="${PHASING_T2T_THREADS:-8}" \
    --mount "type=bind,src=$PROJECT_DIR,dst=/phasing_T2T_project,readonly" \
    --mount "type=bind,src=$INPUT_DIR/intermediate_data_whole_genome,dst=/phasing_T2T_project/intermediate_data_whole_genome,readonly" \
    --mount "type=bind,src=$INPUT_DIR/imputation_statistics_whole_genome,dst=/phasing_T2T_project/imputation_statistics_whole_genome,readonly" \
    --mount "type=bind,src=$OUTPUT_DIR/figures_whole_genome,dst=/phasing_T2T_project/figures_whole_genome" \
    --mount "type=bind,src=$OUTPUT_DIR/tables_whole_genome,dst=/phasing_T2T_project/tables_whole_genome" \
    --mount "type=bind,src=$OUTPUT_DIR/notebook_runs,dst=/notebook_runs" \
    --workdir /phasing_T2T_project "${PHASING_T2T_IMAGE:-jlalli/phasing_t2t:v2.1}" \
    jupyter nbconvert --to notebook --execute notebooks/notebooks_whole_genome/make_plots.ipynb \
        --output /notebook_runs/make_plots.executed.ipynb --ExecutePreprocessor.timeout=-1
