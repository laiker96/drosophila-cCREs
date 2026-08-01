#!/usr/bin/env bash
#SBATCH --job-name=atlas-env
#SBATCH --partition=cpu
#SBATCH --exclude=nodo3
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --propagate=NONE
#SBATCH --output=slurm-%x-%j.out

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-/grupos/Lucai/ANTARKA/drosophila-cCREs}"
HOST_SHORT="$(hostname -s)"
OPEN_FILE_LIMIT="${OPEN_FILE_LIMIT:-65536}"
MIN_OPEN_FILE_LIMIT="${MIN_OPEN_FILE_LIMIT:-8192}"

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo "Refusing to install outside a SLURM allocation." >&2
    exit 97
fi
case "$HOST_SHORT" in
    nodo3|cranex*|cluster|cluster.*|login*)
        echo "Refusing to install on forbidden host $HOST_SHORT." >&2
        exit 97
        ;;
esac

for value in "$OPEN_FILE_LIMIT" "$MIN_OPEN_FILE_LIMIT"; do
    if [[ ! "$value" =~ ^[1-9][0-9]*$ ]]; then
        echo "Open-file limits must be positive integers, got: $value" >&2
        exit 2
    fi
done
if (( MIN_OPEN_FILE_LIMIT > OPEN_FILE_LIMIT )); then
    echo "MIN_OPEN_FILE_LIMIT cannot exceed OPEN_FILE_LIMIT." >&2
    exit 2
fi

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

for command_name in curl flock install sha256sum tar; do
    if ! command -v "$command_name" >/dev/null 2>&1; then
        echo "Required bootstrap command is unavailable: $command_name" >&2
        exit 127
    fi
done
if [[ ! -f "$REPO_ROOT/environment.yml" ]] || [[ ! -d "$REPO_ROOT/workflow/envs" ]]; then
    echo "Repository not found or incomplete: $REPO_ROOT" >&2
    exit 2
fi

cd "$REPO_ROOT"
mkdir -p logs/slurm .cache .micromamba/bin .micromamba/pkgs \
    .micromamba/envs .snakemake/conda
exec > >(tee -a "$REPO_ROOT/logs/slurm/install-environment-${SLURM_JOB_ID}.log") 2>&1

raise_open_file_limit

exec 9>"$REPO_ROOT/.install-environment.lock"
if ! flock -n 9; then
    echo "Another environment installer holds $REPO_ROOT/.install-environment.lock" >&2
    exit 75
fi

TMP_PARENT="${SLURM_TMPDIR:-${TMPDIR:-/tmp}}"
if [[ ! -d "$TMP_PARENT" ]]; then
    echo "Temporary directory does not exist: $TMP_PARENT" >&2
    exit 2
fi
JOB_TMP="$(mktemp -d "$TMP_PARENT/atlas-env-${SLURM_JOB_ID}.XXXXXX")"
cleanup() {
    rm -rf -- "$JOB_TMP"
}
trap cleanup EXIT

export MAMBA_ROOT_PREFIX="$REPO_ROOT/.micromamba"
export CONDA_PKGS_DIRS="$REPO_ROOT/.micromamba/pkgs"
export CONDA_ENVS_PATH="$REPO_ROOT/.micromamba/envs"
export CONDARC="$REPO_ROOT/.condarc"
export XDG_CACHE_HOME="$REPO_ROOT/.cache"
MICROMAMBA="$REPO_ROOT/.micromamba/bin/micromamba"
BUNDLE_ROOT="$REPO_ROOT/.cluster-bootstrap/atlas-offline"
BUNDLED_MICROMAMBA="$BUNDLE_ROOT/tool/bin/micromamba"
PACKED_ORCHESTRATION="$BUNDLE_ROOT/orchestration.tar.gz"
PACKED_CONDA_ENVS="$BUNDLE_ROOT/conda-envs"
BUNDLE_CHECKSUMS="$BUNDLE_ROOT/offline-bundle.sha256"

install_packed_prefix() {
    local archive="$1"
    local target="$2"
    local ready_file="$3"
    local temporary="${target}.tmp-${SLURM_JOB_ID}"
    if [[ -f "$ready_file" ]] && [[ -d "$target/conda-meta" ]]; then
        echo "Packed environment already ready: $target"
        return
    fi
    if [[ -e "$target" ]]; then
        echo "Incomplete packed environment already exists: $target" >&2
        exit 2
    fi
    if [[ -e "$temporary" ]]; then
        echo "Stale staging path exists: $temporary" >&2
        exit 2
    fi
    mkdir -p "$temporary"
    tar -xzf "$archive" -C "$temporary"
    if [[ ! -d "$temporary/conda-meta" ]]; then
        echo "Packed environment lacks conda-meta: $archive" >&2
        exit 2
    fi
    mv "$temporary" "$target"
    if [[ -x "$target/bin/conda-unpack" ]]; then
        env PATH="$target/bin:$PATH" "$target/bin/conda-unpack"
    fi
    touch "$ready_file"
}

echo "[$(date --iso-8601=seconds)] host=$HOST_SHORT job=$SLURM_JOB_ID"
echo "repo=$REPO_ROOT"

OFFLINE_BUNDLE=false
if [[ -f "$BUNDLE_CHECKSUMS" ]] && [[ -f "$PACKED_ORCHESTRATION" ]]; then
    echo "Verifying the staged offline environment bundle."
    (
        cd "$BUNDLE_ROOT"
        sha256sum --check "$(basename "$BUNDLE_CHECKSUMS")"
    )
    OFFLINE_BUNDLE=true
    export CONDA_OFFLINE=true
fi

if [[ ! -x "$MICROMAMBA" ]] && [[ -x "$BUNDLED_MICROMAMBA" ]]; then
    install -m 0755 "$BUNDLED_MICROMAMBA" "$MICROMAMBA"
elif [[ ! -x "$MICROMAMBA" ]] && [[ "$OFFLINE_BUNDLE" != true ]]; then
    echo "Installing Micromamba inside the repository."
    curl --fail --location --retry 5 \
        https://micro.mamba.pm/api/micromamba/linux-64/latest \
        | tar -xj -C "$JOB_TMP" bin/micromamba
    install -m 0755 "$JOB_TMP/bin/micromamba" "$MICROMAMBA"
fi

if [[ -e "$REPO_ROOT/.venv" ]] && [[ ! -d "$REPO_ROOT/.venv/conda-meta" ]]; then
    echo "$REPO_ROOT/.venv exists but is not a Conda environment; move it aside first." >&2
    exit 2
fi
if [[ "$OFFLINE_BUNDLE" == true ]] && [[ ! -d "$REPO_ROOT/.venv/conda-meta" ]]; then
    echo "Installing the packed orchestration environment."
    install_packed_prefix \
        "$PACKED_ORCHESTRATION" \
        "$REPO_ROOT/.venv" \
        "$REPO_ROOT/.venv/.atlas_packed_ready"
elif [[ "$OFFLINE_BUNDLE" == true ]]; then
    echo "Reusing the existing packed orchestration environment."
elif [[ -d "$REPO_ROOT/.venv/conda-meta" ]]; then
    echo "Updating the repository orchestration environment."
    "$MICROMAMBA" install --yes --prefix "$REPO_ROOT/.venv" \
        --file "$REPO_ROOT/environment.yml" --channel-priority strict
else
    echo "Creating the repository orchestration environment."
    "$MICROMAMBA" create --yes --prefix "$REPO_ROOT/.venv" \
        --file "$REPO_ROOT/environment.yml" --channel-priority strict
fi

if [[ "$OFFLINE_BUNDLE" == true ]]; then
    if ! compgen -G "$PACKED_CONDA_ENVS/*.tar.gz" >/dev/null; then
        echo "Offline bundle contains no packed Snakemake environments." >&2
        exit 2
    fi
    for archive in "$PACKED_CONDA_ENVS"/*.tar.gz; do
        env_name="$(basename "$archive" .tar.gz)"
        target="$REPO_ROOT/.snakemake/conda/$env_name"
        ready_file="$REPO_ROOT/.snakemake/conda/${env_name}.env_setup_done"
        install_packed_prefix "$archive" "$target" "$ready_file"
    done
fi

mapfile -t RULE_ENVS < <(find "$REPO_ROOT/workflow/envs" -maxdepth 1 \
    -type f -name '*.yaml' -print | sort)
if (( ${#RULE_ENVS[@]} == 0 )); then
    echo "No Snakemake rule environments were found." >&2
    exit 2
fi

ENV_SNAKEFILE="$JOB_TMP/install-rule-envs.smk"
{
    echo "rule all:"
    echo "    input:"
    for index in "${!RULE_ENVS[@]}"; do
        printf '        "%s/env-%s.ready",\n' "$JOB_TMP" "$index"
    done
    for index in "${!RULE_ENVS[@]}"; do
        printf '\nrule env_%s:\n' "$index"
        printf '    output: "%s/env-%s.ready"\n' "$JOB_TMP" "$index"
        printf '    conda: "%s"\n' "${RULE_ENVS[$index]}"
        printf '    shell: "touch {output}"\n'
    done
} > "$ENV_SNAKEFILE"

echo "Pre-creating ${#RULE_ENVS[@]} Snakemake rule environments."
"$REPO_ROOT/.venv/bin/snakemake" \
    --snakefile "$ENV_SNAKEFILE" \
    --cores 1 \
    --software-deployment-method conda \
    --conda-prefix "$REPO_ROOT/.snakemake/conda" \
    --conda-frontend conda \
    --conda-create-envs-only \
    --printshellcmds

"$REPO_ROOT/.venv/bin/python" --version
"$REPO_ROOT/.venv/bin/conda" --version
"$REPO_ROOT/.venv/bin/snakemake" --version
echo "[$(date --iso-8601=seconds)] complete environment installed under $REPO_ROOT"
