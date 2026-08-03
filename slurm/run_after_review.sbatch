#!/usr/bin/env bash
#SBATCH --job-name=atlas-after-review
#SBATCH --partition=cpu
#SBATCH --exclude=nodo3
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=7-00:00:00
#SBATCH --propagate=NONE
#SBATCH --output=slurm-%x-%j.out

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-/grupos/Lucai/ANTARKA/drosophila-cCREs}"
HOST_SHORT="$(hostname -s)"
OPEN_FILE_LIMIT="${OPEN_FILE_LIMIT:-65536}"
MIN_OPEN_FILE_LIMIT="${MIN_OPEN_FILE_LIMIT:-8192}"

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo "Refusing to process outside a SLURM allocation." >&2
    exit 97
fi
case "$HOST_SHORT" in
    nodo3|cranex*|cluster|cluster.*|login*)
        echo "Refusing to process on forbidden host $HOST_SHORT." >&2
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
MASTER_RUN_ID="${MASTER_RUN_ID:-master-v1}"
FINAL_RUN_ID="${FINAL_RUN_ID:-catalog-v1}"
GENOME="${GENOME:-dm6}"
CONFIG_ROOT="${CONFIG_ROOT:-$REPO_ROOT/configs/cluster}"
TOTAL_CORES="${TOTAL_CORES:-${SLURM_CPUS_PER_TASK:-8}}"
TOTAL_MEMORY_MB="${TOTAL_MEMORY_MB:-${SLURM_MEM_PER_NODE:-32768}}"
MAX_JOBS="${MAX_JOBS:-$TOTAL_CORES}"
MAX_THREADS="${MAX_THREADS:-$TOTAL_CORES}"
REVIEW_ROOT="${REVIEW_ROOT:-$REPO_ROOT/data/reviewed}"
ATAC_REVIEW="${ATAC_REVIEW:-$REVIEW_ROOT/atlas-atac.library-review.tsv}"
H3_REVIEW="${H3_REVIEW:-$REVIEW_ROOT/atlas-h3k27ac.library-review.tsv}"
ATAC_REVIEWED_MANIFEST="${ATAC_REVIEWED_MANIFEST:-$REVIEW_ROOT/atlas-atac.reviewed.tsv}"
H3_REVIEWED_MANIFEST="${H3_REVIEWED_MANIFEST:-$REVIEW_ROOT/atlas-h3k27ac.reviewed.tsv}"

for identifier in "$PROJECT" "$QC_RUN_ID" "$MASTER_RUN_ID" "$FINAL_RUN_ID"; do
    if [[ ! "$identifier" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]]; then
        echo "Invalid identifier: $identifier" >&2
        exit 2
    fi
done
for value in "$TOTAL_CORES" "$TOTAL_MEMORY_MB" "$MAX_JOBS" "$MAX_THREADS" "$OPEN_FILE_LIMIT" \
    "$MIN_OPEN_FILE_LIMIT"; do
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
if [[ "$MASTER_RUN_ID" == "$FINAL_RUN_ID" ]]; then
    echo "MASTER_RUN_ID and FINAL_RUN_ID must name different result namespaces." >&2
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

ATAC_QC_ROOT="$REPO_ROOT/results/$PROJECT.atac.$GENOME/$QC_RUN_ID"
H3_QC_ROOT="$REPO_ROOT/results/$PROJECT.chip_histone.$GENOME/$QC_RUN_ID"
ATAC_FINAL_MANIFEST="$ATAC_QC_ROOT/provenance/manifests/final-bams.tsv"
H3_FINAL_MANIFEST="$H3_QC_ROOT/provenance/manifests/final-bams.tsv"
ATAC_QC_CHECKPOINT="$ATAC_QC_ROOT/provenance/manifests/qc.checkpoint.json"
MASTER_ROOT="$REPO_ROOT/results/$PROJECT/$MASTER_RUN_ID"
MASTER_MANIFEST="$MASTER_ROOT/provenance/manifests/master-dhs.tsv"
FINAL_ROOT="$REPO_ROOT/results/$PROJECT/$FINAL_RUN_ID"

for required_file in "$ATAC_REVIEW" "$H3_REVIEW" "$ATAC_FINAL_MANIFEST" \
    "$H3_FINAL_MANIFEST" "$ATAC_QC_CHECKPOINT"; do
    if [[ ! -f "$required_file" ]]; then
        echo "Required reviewed or QC file does not exist: $required_file" >&2
        exit 2
    fi
done

cd "$REPO_ROOT"
mkdir -p logs/slurm "$CONFIG_ROOT" "$REVIEW_ROOT" .cache
exec > >(tee -a "$REPO_ROOT/logs/slurm/run-after-review-${SLURM_JOB_ID}.log") 2>&1

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
JOB_TMP="$(mktemp -d "$TMP_PARENT/atlas-after-review-${SLURM_JOB_ID}.XXXXXX")"
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
echo "repo=$REPO_ROOT project=$PROJECT master_run=$MASTER_RUN_ID final_run=$FINAL_RUN_ID"
echo "execution=single-node jobs=$MAX_JOBS cores=$TOTAL_CORES mem_mb=$TOTAL_MEMORY_MB nofile=$(ulimit -Sn)"

echo "Applying the completed manual review tables."
"$REPO_ROOT/.venv/bin/python" "$REPO_ROOT/src/review_final_bam_manifest.py" \
    "$ATAC_FINAL_MANIFEST" \
    --review-table "$ATAC_REVIEW" \
    --output "$ATAC_REVIEWED_MANIFEST"
"$REPO_ROOT/.venv/bin/python" "$REPO_ROOT/src/review_final_bam_manifest.py" \
    "$H3_FINAL_MANIFEST" \
    --review-table "$H3_REVIEW" \
    --output "$H3_REVIEWED_MANIFEST"

MASTER_COMMAND=(
    "$REPO_ROOT/.venv/bin/python"
    "$REPO_ROOT/src/run_pipeline.py"
    "$SAMPLE_SHEET"
    --from-stage qc
    --checkpoint-manifest "$ATAC_QC_CHECKPOINT"
    --final-bam-manifest "$ATAC_REVIEWED_MANIFEST"
    --project "$PROJECT"
    --run-id "$MASTER_RUN_ID"
    --genome "$GENOME"
    --until-stage master
    --config-dir "$CONFIG_ROOT"
    --reference-root "$REPO_ROOT/references"
    --workflow-profile "$PROFILE_DIR"
    --cores "$TOTAL_CORES"
    --jobs "$MAX_JOBS"
    --max-threads "$MAX_THREADS"
    --snakemake-arg=--retries
    --snakemake-arg=2
)
printf 'Running:'
printf ' %q' "${MASTER_COMMAND[@]}"
printf '\n'
"${MASTER_COMMAND[@]}"

if [[ ! -f "$MASTER_MANIFEST" ]]; then
    echo "Master DHS manifest was not produced: $MASTER_MANIFEST" >&2
    exit 2
fi

FINAL_COMMAND=(
    "$REPO_ROOT/.venv/bin/python"
    "$REPO_ROOT/src/run_pipeline.py"
    "$SAMPLE_SHEET"
    --from-stage master
    --master-manifest "$MASTER_MANIFEST"
    --activity-bam-manifest "$ATAC_REVIEWED_MANIFEST"
    --activity-bam-manifest "$H3_REVIEWED_MANIFEST"
    --report-source-root "$ATAC_QC_ROOT"
    --report-source-root "$H3_QC_ROOT"
    --report-source-root "$MASTER_ROOT"
    --project "$PROJECT"
    --run-id "$FINAL_RUN_ID"
    --genome "$GENOME"
    --until-stage report
    --config-dir "$CONFIG_ROOT"
    --reference-root "$REPO_ROOT/references"
    --workflow-profile "$PROFILE_DIR"
    --cores "$TOTAL_CORES"
    --jobs "$MAX_JOBS"
    --max-threads "$MAX_THREADS"
    --snakemake-arg=--retries
    --snakemake-arg=2
)
printf 'Running:'
printf ' %q' "${FINAL_COMMAND[@]}"
printf '\n'
"${FINAL_COMMAND[@]}"

cat <<EOF

[$(date --iso-8601=seconds)] post-review pipeline complete.
Master DHS manifest:
  $MASTER_MANIFEST
Regulatory-element tables and IGV sessions:
  $FINAL_ROOT/activity/catalog
Candidate-gene link tables:
  $FINAL_ROOT/activity/links
Integrated report:
  $FINAL_ROOT/activity/report/integrated_qc_report.html
EOF
