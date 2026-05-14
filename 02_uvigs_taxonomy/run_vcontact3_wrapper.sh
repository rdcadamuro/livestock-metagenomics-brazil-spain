#!/usr/bin/env bash
set -euo pipefail

# USER: set your base path here

ENV_PATH="/mnt/nvme2/conda_envs/vcontact3"
# USER: set your base path here
DB_PATH="/mnt/nvme2/vcontact3_db"
THREADS="${THREADS:-24}"

usage() {
    cat <<'EOF'
Usage:
  run_vcontact3_wrapper.sh nucleotide <input_fasta> <output_dir> [extra vcontact3 args...]
  run_vcontact3_wrapper.sh proteins <proteins_faa> <gene2genome.tsv> <output_dir> [len_nucleotide.tsv] [extra vcontact3 args...]

Examples:
  run_vcontact3_wrapper.sh nucleotide bins.fna results_vcontact --db-domain prokaryotes --exports graphml cytoscape
  run_vcontact3_wrapper.sh proteins proteins.faa gene2genome.tsv results_vcontact lengths.tsv --db-domain prokaryotes

Notes:
  - Uses the environment at /mnt/nvme2/conda_envs/vcontact3
  - Uses the database directory at /mnt/nvme2/vcontact3_db
  - Default threads comes from THREADS env var or 24
EOF
}

if [[ $# -lt 1 ]]; then
    usage
    exit 1
fi

MODE="$1"
shift

if [[ ! -x "${ENV_PATH}/bin/vcontact3" ]]; then
    echo "vcontact3 not found in ${ENV_PATH}" >&2
    exit 1
fi

if [[ ! -e "${DB_PATH}" ]]; then
    echo "Database path not found: ${DB_PATH}" >&2
    exit 1
fi

BASE_CMD=(
    /home/rafael/miniforge3/bin/conda run -p "${ENV_PATH}" vcontact3 run
    --threads "${THREADS}"
    --db-path "${DB_PATH}"
    --db-domain prokaryotes
    --distance-metric SqRoot
    --max-iterations 3
    --reduce-memory
    --force-overwrite
)

case "${MODE}" in
    nucleotide)
        if [[ $# -lt 2 ]]; then
            usage
            exit 1
        fi
        INPUT_FASTA="$1"
        OUTPUT_DIR="$2"
        shift 2
        exec "${BASE_CMD[@]}" --nucleotide "${INPUT_FASTA}" --output "${OUTPUT_DIR}" "$@"
        ;;
    proteins)
        if [[ $# -lt 3 ]]; then
            usage
            exit 1
        fi
        PROTEINS="$1"
        GENE2GENOME="$2"
        OUTPUT_DIR="$3"
        shift 3

        OPTS=(--proteins "${PROTEINS}" --gene2genome "${GENE2GENOME}" --output "${OUTPUT_DIR}")
        if [[ $# -ge 1 && -f "$1" && ( "$1" == *.tsv || "$1" == *.parquet ) ]]; then
            OPTS+=(--len-nucleotide "$1")
            shift
        fi
        exec "${BASE_CMD[@]}" "${OPTS[@]}" "$@"
        ;;
    *)
        usage
        exit 1
        ;;
esac
