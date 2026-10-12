#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"

usage() {
    echo 'Usage: bash scripts/run_make_plots.sh'
    echo 'Run make_plots.ipynb in Docker using inputs and output folders in the repository.'
    echo 'Set PHASING_T2T_IMAGE to use another tag (default: jlalli/phasing_t2t:v2.1).'
}

if [[ $# -eq 1 && ( $1 == -h || $1 == --help ) ]]; then
    usage
    exit 0
fi
[[ $# -eq 0 ]] || { usage >&2; exit 2; }

for directory in intermediate_data_whole_genome imputation_statistics_whole_genome; do
    [[ -d "$PROJECT_DIR/$directory" ]] || { echo "ERROR: missing input directory: $directory" >&2; exit 2; }
done
(cd -- "$PROJECT_DIR" && sha256sum --check --quiet SHA256SUMS.txt)

for directory in intermediate_data_whole_genome imputation_statistics_whole_genome figures_whole_genome tables_whole_genome; do
    [[ ! -L "$PROJECT_DIR/$directory" ]] || { echo "ERROR: mount destination is a symlink: $directory" >&2; exit 2; }
    mkdir -p -- "$PROJECT_DIR/$directory"
done
mkdir -p -- "$PROJECT_DIR/notebook_runs"

exec docker run --rm --network none --read-only --tmpfs /tmp:rw,size=8g \
    --user "$(id -u):$(id -g)" --cpus "${PHASING_T2T_THREADS:-8}" --memory "${PHASING_T2T_MEMORY:-64g}" \
    -e MPLBACKEND=Agg -e MPLCONFIGDIR=/tmp/matplotlib -e XDG_CACHE_HOME=/tmp/cache \
    -e JUPYTER_DATA_DIR=/tmp/jupyter-data -e JUPYTER_CONFIG_DIR=/tmp/jupyter-config \
    -e JUPYTER_RUNTIME_DIR=/tmp/jupyter-runtime -e IPYTHONDIR=/tmp/ipython \
    -e PHASING_T2T_RUN_PROFILE=whole_genome -e NUM_THREADS="${PHASING_T2T_THREADS:-8}" \
    -e OPENBLAS_NUM_THREADS="${PHASING_T2T_THREADS:-8}" -e MKL_NUM_THREADS="${PHASING_T2T_THREADS:-8}" \
    -e OMP_NUM_THREADS="${PHASING_T2T_THREADS:-8}" -e POLARS_MAX_THREADS="${PHASING_T2T_THREADS:-8}" \
    --mount "type=bind,src=$PROJECT_DIR,dst=/phasing_T2T_project" \
    --workdir /phasing_T2T_project "${PHASING_T2T_IMAGE:-jlalli/phasing_t2t:v2.1}" \
    jupyter nbconvert --to notebook --execute notebooks/notebooks_whole_genome/make_plots.ipynb \
        --output /phasing_T2T_project/notebook_runs/make_plots.executed.ipynb --ExecutePreprocessor.timeout=-1
