#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
IMAGE_TAG="phasing-t2t-figures:reported"
CPUS=8
MEMORY=64g
BUILD=0
INPUT_DIR=
OUTPUT_DIR=

usage() {
    cat <<'EOF'
Usage: scripts/run_publication_notebook.sh --inputs INPUTDIR --output OUTPUTDIR [--image TAG] [--build] [--cpus N] [--memory SIZE]

Execute the canonical whole-genome publication notebook in an isolated Docker
container. INPUTDIR must be the extracted NG-TR68126R notebook-input release;
OUTPUTDIR must be new or contain no prior figure, table, or run receipt data.
EOF
}

while (($#)); do
    case "$1" in
        --inputs) INPUT_DIR=${2:?missing INPUTDIR}; shift 2 ;;
        --output) OUTPUT_DIR=${2:?missing OUTPUTDIR}; shift 2 ;;
        --image) IMAGE_TAG=${2:?missing TAG}; shift 2 ;;
        --build) BUILD=1; shift ;;
        --cpus) CPUS=${2:?missing N}; shift 2 ;;
        --memory) MEMORY=${2:?missing SIZE}; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done
[[ "$CPUS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: --cpus must be a positive integer" >&2; exit 2; }
[[ -n "$INPUT_DIR" && -n "$OUTPUT_DIR" ]] || { usage >&2; exit 2; }
[[ -d "$INPUT_DIR" ]] || { echo "ERROR: input directory does not exist: $INPUT_DIR" >&2; exit 2; }
INPUT_DIR="$(cd "$INPUT_DIR" && pwd -P)"
OUTPUT_DIR="$(mkdir -p "$OUTPUT_DIR" && cd "$OUTPUT_DIR" && pwd -P)"
for required in INPUT_MANIFEST.json SHA256SUMS.txt README.md intermediate_data_whole_genome imputation_statistics_whole_genome; do
    [[ -e "$INPUT_DIR/$required" ]] || { echo "ERROR: missing input release member: $required" >&2; exit 2; }
done
for output in figures_whole_genome tables_whole_genome scratch; do
    if [[ -e "$OUTPUT_DIR/$output" ]] && find "$OUTPUT_DIR/$output" -mindepth 1 -print -quit | grep -q .; then
        echo "ERROR: output already contains data: $OUTPUT_DIR/$output" >&2; exit 2
    fi
    mkdir -p "$OUTPUT_DIR/$output"
done

BUILD_DIR=
cleanup() { [[ -z "$BUILD_DIR" ]] || rm -rf -- "$BUILD_DIR"; }
trap cleanup EXIT
if (( BUILD )) || ! docker image inspect "$IMAGE_TAG" >/dev/null 2>&1; then
    BUILD_DIR="$(mktemp -d)"
    cp "$PROJECT_DIR/Dockerfile.reported-versions" "$PROJECT_DIR/requirements.figures.lock.txt" "$BUILD_DIR/"
    docker build -f "$BUILD_DIR/Dockerfile.reported-versions" -t "$IMAGE_TAG" "$BUILD_DIR"
fi
IMAGE_ID="$(docker image inspect --format '{{.Id}}' "$IMAGE_TAG")"
printf '%q ' "$0" --inputs "$INPUT_DIR" --output "$OUTPUT_DIR" --image "$IMAGE_TAG" --cpus "$CPUS" --memory "$MEMORY" > "$OUTPUT_DIR/scratch/run_command.txt"
printf '\n' >> "$OUTPUT_DIR/scratch/run_command.txt"
STARTED="$(date +%s)"
set +e
docker run --rm --network none --read-only --tmpfs /tmp:rw,size=8g \
    --user "$(id -u):$(id -g)" --cpus "$CPUS" --memory "$MEMORY" \
    -e MPLBACKEND=Agg -e MPLCONFIGDIR=/tmp/matplotlib -e XDG_CACHE_HOME=/tmp/cache \
    -e JUPYTER_DATA_DIR=/tmp/jupyter-data -e JUPYTER_CONFIG_DIR=/tmp/jupyter-config \
    -e JUPYTER_RUNTIME_DIR=/tmp/jupyter-runtime -e IPYTHONDIR=/tmp/ipython \
    -e PHASING_T2T_RUN_PROFILE=whole_genome -e NUM_THREADS="$CPUS" \
    -e OPENBLAS_NUM_THREADS="$CPUS" -e MKL_NUM_THREADS="$CPUS" -e OMP_NUM_THREADS="$CPUS" -e POLARS_MAX_THREADS="$CPUS" \
    -e IMAGE_ID="$IMAGE_ID" -e STARTED="$STARTED" -e RUN_CPUS="$CPUS" -e RUN_MEMORY="$MEMORY" \
    -e MOUNT_REPO_SOURCE="$PROJECT_DIR" -e MOUNT_INPUT_SOURCE="$INPUT_DIR" \
    -e MOUNT_INTERMEDIATE_SOURCE="$INPUT_DIR/intermediate_data_whole_genome" \
    -e MOUNT_IMPUTATION_SOURCE="$INPUT_DIR/imputation_statistics_whole_genome" \
    -e MOUNT_FIGURES_SOURCE="$OUTPUT_DIR/figures_whole_genome" \
    -e MOUNT_TABLES_SOURCE="$OUTPUT_DIR/tables_whole_genome" -e MOUNT_SCRATCH_SOURCE="$OUTPUT_DIR/scratch" \
    -v "$PROJECT_DIR:/work/repo:ro" -v "$INPUT_DIR:/inputs:ro" \
    -v "$INPUT_DIR/intermediate_data_whole_genome:/work/repo/intermediate_data_whole_genome:ro" \
    -v "$INPUT_DIR/imputation_statistics_whole_genome:/work/repo/imputation_statistics_whole_genome:ro" \
    -v "$OUTPUT_DIR/figures_whole_genome:/work/repo/figures_whole_genome" \
    -v "$OUTPUT_DIR/tables_whole_genome:/work/repo/tables_whole_genome" \
    -v "$OUTPUT_DIR/scratch:/scratch" -w /work/repo "$IMAGE_TAG" bash -lc '
        set +e
        python scripts/analysis/verify_notebook_reproduction.py --preflight --write-run-recipe --scratch /scratch --image-id "$IMAGE_ID"
        preflight_status=$?
        run_status=$preflight_status
        if [ "$preflight_status" -ne 0 ]; then exit "$preflight_status"; fi
        if [ "$preflight_status" -eq 0 ]; then
            jupyter nbconvert --to notebook --execute notebooks/notebooks_whole_genome/make_plots.ipynb --output /scratch/make_plots.executed.ipynb --ExecutePreprocessor.timeout=-1
            run_status=$?
        fi
        python scripts/analysis/verify_notebook_reproduction.py --verify --executed /scratch/make_plots.executed.ipynb --figures /work/repo/figures_whole_genome --tables /work/repo/tables_whole_genome --scratch /scratch --exit-status "$run_status" --started "$STARTED" --image-id "$IMAGE_ID"
        exit $?
    ' 2>&1 | tee "$OUTPUT_DIR/scratch/run.log"
PIPE_STATUS=${PIPESTATUS[0]}
set -e
printf '%s\n' "$PIPE_STATUS" > "$OUTPUT_DIR/scratch/docker_exit_status.txt"
exit "$PIPE_STATUS"
