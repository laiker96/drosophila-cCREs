#!/usr/bin/env bash
#SBATCH --job-name=atlas-to-review
#SBATCH --partition=cpu
#SBATCH --exclude=nodo3,nodo4
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=7-00:00:00
#SBATCH --propagate=NONE
#SBATCH --output=slurm-%x-%j.out

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-/grupos/Lucai/ANTARKA/drosophila-cCREs}"
HOST_SHORT="$(hostname -s)"
OPEN_FILE_LIMIT="${OPEN_FILE_LIMIT:-65536}"
MIN_OPEN_FILE_LIMIT="${MIN_OPEN_FILE_LIMIT:-8192}"

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo "Refusing to download or process outside a SLURM allocation." >&2
    exit 97
fi
case "$HOST_SHORT" in
    nodo3|nodo4|cranex*|cluster|cluster.*|login*)
        echo "Refusing to download or process on forbidden host $HOST_SHORT." >&2
        exit 97
        ;;
esac

raise_open_file_limit() {
    local soft_before hard_limit target soft_after
    soft_before="$(ulimit -Sn)"
    hard_limit="$(ulimit -Hn)"
    target="$OPEN_FILE_LIMIT"
    if [[ "$hard_limit" != unlimited ]] && (( target > hard_limit )); then
        target="$hard_limit"
    fi
    if [[ "$soft_before" != unlimited ]] && (( soft_before < target )); then
        if ! ulimit -Sn "$target"; then
            echo "Could not raise the open-file limit from $soft_before to $target." >&2
            exit 2
        fi
    fi
    soft_after="$(ulimit -Sn)"
    echo "open_file_limit soft_before=$soft_before hard=$hard_limit applied=$soft_after"
    if [[ "$soft_after" != unlimited ]] && (( soft_after < MIN_OPEN_FILE_LIMIT )); then
        echo "Open-file limit $soft_after is below required minimum $MIN_OPEN_FILE_LIMIT." >&2
        echo "Ask the cluster administrator to raise the compute-node hard NOFILE limit." >&2
        exit 2
    fi
}

SAMPLE_SHEET="${SAMPLE_SHEET:-$REPO_ROOT/resources/atlas_samples_ip_only.tsv}"
PROJECT="${PROJECT:-drosophila-atlas}"
QC_RUN_ID="${QC_RUN_ID:-qc-v1}"
GENOME="${GENOME:-dm6}"
RAW_ROOT="${RAW_ROOT:-$REPO_ROOT/data/raw/drosophila-atlas}"
MANIFEST="${MANIFEST:-$RAW_ROOT/download_manifest.tsv}"
CONFIG_ROOT="${CONFIG_ROOT:-$REPO_ROOT/configs/cluster}"
TOTAL_CORES="${TOTAL_CORES:-${SLURM_CPUS_PER_TASK:-8}}"
TOTAL_MEMORY_MB="${TOTAL_MEMORY_MB:-${SLURM_MEM_PER_NODE:-16384}}"
MAX_JOBS="${MAX_JOBS:-$TOTAL_CORES}"
MAX_THREADS="${MAX_THREADS:-$TOTAL_CORES}"
DOWNLOAD_BACKEND="${DOWNLOAD_BACKEND:-auto}"
RESOLVE_JOBS="${RESOLVE_JOBS:-8}"
FILE_JOBS="${FILE_JOBS:-4}"
CONNECTIONS="${CONNECTIONS:-2}"
CHECKSUM_RETRIES="${CHECKSUM_RETRIES:-5}"
SRA_JOBS="${SRA_JOBS:-2}"
DOWNLOAD_THREADS="${DOWNLOAD_THREADS:-${SLURM_CPUS_PER_TASK:-8}}"

for identifier in "$PROJECT" "$QC_RUN_ID"; do
    if [[ ! "$identifier" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]]; then
        echo "Invalid identifier: $identifier" >&2
        exit 2
    fi
done
for value in "$TOTAL_CORES" "$TOTAL_MEMORY_MB" "$MAX_JOBS" "$MAX_THREADS" "$RESOLVE_JOBS" \
    "$FILE_JOBS" "$CONNECTIONS" "$CHECKSUM_RETRIES" "$SRA_JOBS" "$DOWNLOAD_THREADS" \
    "$OPEN_FILE_LIMIT" "$MIN_OPEN_FILE_LIMIT"; do
    if [[ ! "$value" =~ ^[1-9][0-9]*$ ]]; then
        echo "Expected a positive integer, got: $value" >&2
        exit 2
    fi
done
if (( MIN_OPEN_FILE_LIMIT > OPEN_FILE_LIMIT )); then
    echo "MIN_OPEN_FILE_LIMIT cannot exceed OPEN_FILE_LIMIT." >&2
    exit 2
fi
if [[ "$GENOME" != dm6 && "$GENOME" != hg38 ]]; then
    echo "GENOME must be dm6 or hg38." >&2
    exit 2
fi
if [[ "$DOWNLOAD_BACKEND" != auto && "$DOWNLOAD_BACKEND" != ena \
    && "$DOWNLOAD_BACKEND" != sra ]]; then
    echo "DOWNLOAD_BACKEND must be auto, ena, or sra." >&2
    exit 2
fi
if [[ ! -f "$SAMPLE_SHEET" ]]; then
    echo "Sample sheet does not exist: $SAMPLE_SHEET" >&2
    exit 2
fi
if [[ ! -x "$REPO_ROOT/.venv/bin/python" ]] \
    || [[ ! -x "$REPO_ROOT/.venv/bin/snakemake" ]]; then
    echo "The repository environment is missing. Submit slurm/install_environment.sbatch first." >&2
    exit 2
fi

cd "$REPO_ROOT"
mkdir -p logs/slurm "$RAW_ROOT" "$CONFIG_ROOT" .cache
exec > >(tee -a "$REPO_ROOT/logs/slurm/run-to-review-${SLURM_JOB_ID}.log") 2>&1

raise_open_file_limit

exec 9>"$REPO_ROOT/.atlas-pipeline.lock"
if ! flock -n 9; then
    echo "Another atlas pipeline job holds $REPO_ROOT/.atlas-pipeline.lock" >&2
    exit 75
fi

TMP_PARENT="${SLURM_TMPDIR:-${TMPDIR:-/tmp}}"
if [[ ! -d "$TMP_PARENT" ]]; then
    echo "Temporary directory does not exist: $TMP_PARENT" >&2
    exit 2
fi
JOB_TMP="$(mktemp -d "$TMP_PARENT/atlas-to-review-${SLURM_JOB_ID}.XXXXXX")"
cleanup() {
    rm -rf -- "$JOB_TMP"
}
trap cleanup EXIT

export TMPDIR="$JOB_TMP"
export MAMBA_ROOT_PREFIX="$REPO_ROOT/.micromamba"
export CONDA_PKGS_DIRS="$REPO_ROOT/.micromamba/pkgs"
export CONDA_ENVS_PATH="$REPO_ROOT/.micromamba/envs"
export CONDARC="$REPO_ROOT/.condarc"
export XDG_CACHE_HOME="$REPO_ROOT/.cache"
export PYTHONNOUSERSITE=1
export PATH="$REPO_ROOT/.venv/bin:$PATH"

PROFILE_DIR="$JOB_TMP/local-profile"
mkdir -p "$PROFILE_DIR"
cat > "$PROFILE_DIR/config.yaml" <<EOF
executor: local
software-deployment-method:
  - conda
conda-prefix: "$REPO_ROOT/.snakemake/conda"
conda-frontend: conda
resources:
  mem_mb: $TOTAL_MEMORY_MB
default-resources:
  mem_mb: 2000
latency-wait: 180
printshellcmds: true
show-failed-logs: true
rerun-incomplete: true
EOF

echo "[$(date --iso-8601=seconds)] host=$HOST_SHORT job=$SLURM_JOB_ID"
echo "repo=$REPO_ROOT project=$PROJECT run=$QC_RUN_ID genome=$GENOME"
echo "execution=single-node jobs=$MAX_JOBS cores=$TOTAL_CORES mem_mb=$TOTAL_MEMORY_MB nofile=$(ulimit -Sn)"

PIPELINE_COMMAND=(
    "$REPO_ROOT/.venv/bin/python"
    "$REPO_ROOT/src/run_pipeline.py"
    "$SAMPLE_SHEET"
    --project "$PROJECT"
    --run-id "$QC_RUN_ID"
    --genome "$GENOME"
    --until-stage qc
    --config-dir "$CONFIG_ROOT"
    --reference-root "$REPO_ROOT/references"
    --output-dir "$RAW_ROOT"
    --manifest "$MANIFEST"
    --backend "$DOWNLOAD_BACKEND"
    --resolve-jobs "$RESOLVE_JOBS"
    --file-jobs "$FILE_JOBS"
    --connections "$CONNECTIONS"
    --checksum-retries "$CHECKSUM_RETRIES"
    --sra-jobs "$SRA_JOBS"
    --threads "$DOWNLOAD_THREADS"
    --workflow-profile "$PROFILE_DIR"
    --cores "$TOTAL_CORES"
    --jobs "$MAX_JOBS"
    --max-threads "$MAX_THREADS"
    --snakemake-arg=--retries
    --snakemake-arg=2
)
if [[ -f "$MANIFEST" ]]; then
    echo "Reusing completed download manifest: $MANIFEST"
    PIPELINE_COMMAND+=(--skip-download)
fi
printf 'Running:'
printf ' %q' "${PIPELINE_COMMAND[@]}"
printf '\n'
"${PIPELINE_COMMAND[@]}"

ATAC_ROOT="$REPO_ROOT/results/$PROJECT.atac.$GENOME/$QC_RUN_ID"
H3_ROOT="$REPO_ROOT/results/$PROJECT.chip_histone.$GENOME/$QC_RUN_ID"
cat <<EOF

[$(date --iso-8601=seconds)] QC is complete. The pipeline is stopped for manual review.

ATAC review table:
  $ATAC_ROOT/qc/library-review.tsv
H3K27ac review table:
  $H3_ROOT/qc/library-review.tsv

Copy both tables outside results before editing:
  mkdir -p "$REPO_ROOT/data/reviewed"
  cp "$ATAC_ROOT/qc/library-review.tsv" "$REPO_ROOT/data/reviewed/atlas-atac.library-review.tsv"
  cp "$H3_ROOT/qc/library-review.tsv" "$REPO_ROOT/data/reviewed/atlas-h3k27ac.library-review.tsv"

Edit only qc_decision, estimated_fragment_length_bp, and notes. Use pass or fail
for every qc_decision; failed rows require notes. Only single-end H3K27ac rows
may contain estimated_fragment_length_bp. Then submit run_after_review.sbatch.
EOF
