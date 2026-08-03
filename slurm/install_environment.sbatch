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

for command_name in awk curl flock install sha256sum tar; do
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
PACKED_ENV_MAP="$BUNDLE_ROOT/conda-env-map.tsv"
BUNDLE_CHECKSUMS="$BUNDLE_ROOT/offline-bundle.sha256"

finalize_packed_prefix() {
    local target="$1"
    local ready_file="$2"
    local unpack_python
    if [[ ! -f "$target/bin/conda-unpack" ]]; then
        echo "Packed environment lacks conda-unpack: $target" >&2
        exit 2
    fi
    if [[ -x "$target/bin/python" ]]; then
        unpack_python="$target/bin/python"
    elif [[ -x "$REPO_ROOT/.venv/bin/python" ]]; then
        unpack_python="$REPO_ROOT/.venv/bin/python"
    else
        echo "No Python interpreter is available to finalize: $target" >&2
        exit 2
    fi
    "$unpack_python" "$target/bin/conda-unpack"
    touch "$ready_file"
}

install_packed_prefix() {
    local archive="$1"
    local target="$2"
    local ready_file="$3"
    local temporary="${target}.tmp-${SLURM_JOB_ID}"
    if [[ -d "$target/conda-meta" ]]; then
        if [[ -f "$ready_file" ]]; then
            echo "Packed environment already ready: $target"
            return
        fi
        echo "Finalizing an unpacked environment: $target"
        finalize_packed_prefix "$target" "$ready_file"
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
    finalize_packed_prefix "$target" "$ready_file"
}

install_conda_activation_wrapper() {
    local activate="$REPO_ROOT/.venv/bin/activate"
    local conda_sh="$REPO_ROOT/.venv/etc/profile.d/conda.sh"
    local temporary="${activate}.tmp-${SLURM_JOB_ID}"
    if [[ ! -f "$conda_sh" ]]; then
        echo "Conda shell initialization is missing: $conda_sh" >&2
        exit 2
    fi
    if [[ "$REPO_ROOT" == *"'"* ]] || [[ "$REPO_ROOT" == *$'\n'* ]]; then
        echo "Repository path cannot contain a quote or newline: $REPO_ROOT" >&2
        exit 2
    fi
    {
        printf '%s\n' '#!/bin/sh'
        printf "_CONDA_ROOT='%s'\n" "$REPO_ROOT/.venv"
        printf '%s\n' '. "$_CONDA_ROOT/etc/profile.d/conda.sh" || return $?'
        printf '%s\n' 'conda activate "$@"'
    } > "$temporary"
    chmod 0755 "$temporary"
    mv -f -- "$temporary" "$activate"
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
if [[ "$OFFLINE_BUNDLE" == true ]] && [[ ! -f "$REPO_ROOT/.venv/.atlas_packed_ready" ]]; then
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
install_conda_activation_wrapper
export PATH="$REPO_ROOT/.venv/bin:$PATH"

mapfile -t RULE_ENVS < <(find "$REPO_ROOT/workflow/envs" -maxdepth 1 \
    -type f -name '*.yaml' -print | sort)
if (( ${#RULE_ENVS[@]} == 0 )); then
    echo "No Snakemake rule environments were found." >&2
    exit 2
fi

ENV_SNAKEFILE="$JOB_TMP/install-rule-envs.smk"
ENV_LIST="$JOB_TMP/conda-envs.tsv"
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

if [[ "$OFFLINE_BUNDLE" == true ]]; then
    if [[ ! -f "$PACKED_ENV_MAP" ]]; then
        echo "Offline bundle lacks its Snakemake environment map: $PACKED_ENV_MAP" >&2
        exit 2
    fi
    declare -A PACKED_ARCHIVE_BY_ENV=()
    while IFS=$'\t' read -r mapped_env mapped_archive extra; do
        if [[ -z "$mapped_env" ]] || [[ -z "$mapped_archive" ]] || [[ -n "$extra" ]]; then
            echo "Invalid offline environment-map row: $mapped_env $mapped_archive $extra" >&2
            exit 2
        fi
        if [[ -n "${PACKED_ARCHIVE_BY_ENV[$mapped_env]+present}" ]]; then
            echo "Duplicate offline environment mapping: $mapped_env" >&2
            exit 2
        fi
        if [[ "$mapped_env" == */* ]]; then
            echo "Offline environment-map keys must be basenames: $mapped_env" >&2
            exit 2
        fi
        case "$mapped_archive" in
            conda-envs/*/*.tar.gz)
                echo "Invalid packed environment path: $mapped_archive" >&2
                exit 2
                ;;
            conda-envs/*.tar.gz) ;;
            *)
                echo "Invalid packed environment path: $mapped_archive" >&2
                exit 2
                ;;
        esac
        PACKED_ARCHIVE_BY_ENV["$mapped_env"]="$mapped_archive"
    done < "$PACKED_ENV_MAP"

    "$REPO_ROOT/.venv/bin/snakemake" \
        --snakefile "$ENV_SNAKEFILE" \
        --cores 1 \
        --software-deployment-method conda \
        --conda-prefix "$REPO_ROOT/.snakemake/conda" \
        --conda-frontend conda \
        --list-conda-envs \
        > "$ENV_LIST"

    installed_env_count=0
    while IFS=$'\t' read -r env_path target; do
        env_name="$(basename "$env_path")"
        mapped_archive="${PACKED_ARCHIVE_BY_ENV[$env_name]:-}"
        if [[ -z "$mapped_archive" ]]; then
            echo "No packed environment mapped for $env_name" >&2
            exit 2
        fi
        archive="$BUNDLE_ROOT/$mapped_archive"
        if [[ ! -f "$archive" ]]; then
            echo "Mapped packed environment does not exist: $archive" >&2
            exit 2
        fi
        if [[ "$target" != /* ]]; then
            target="$REPO_ROOT/$target"
        fi
        case "$target" in
            "$REPO_ROOT/.snakemake/conda/"*) ;;
            *)
                echo "Snakemake environment target is outside the repository: $target" >&2
                exit 2
                ;;
        esac
        ready_file="${target}.env_setup_done"
        install_packed_prefix "$archive" "$target" "$ready_file"
        install -m 0644 "$env_path" "${target}.yaml"
        (( installed_env_count += 1 ))
    done < <(awk -F '\t' 'NR > 1 {print $1 "\t" $3}' "$ENV_LIST")
    if (( installed_env_count != ${#RULE_ENVS[@]} )); then
        echo "Expected ${#RULE_ENVS[@]} packed environments, installed $installed_env_count" >&2
        exit 2
    fi
fi

echo "Pre-creating ${#RULE_ENVS[@]} Snakemake rule environments."
"$REPO_ROOT/.venv/bin/snakemake" \
    --snakefile "$ENV_SNAKEFILE" \
    --cores 1 \
    --software-deployment-method conda \
    --conda-prefix "$REPO_ROOT/.snakemake/conda" \
    --conda-frontend conda \
    --conda-create-envs-only \
    --printshellcmds

if [[ ! -s "$ENV_LIST" ]]; then
    "$REPO_ROOT/.venv/bin/snakemake" \
        --snakefile "$ENV_SNAKEFILE" \
        --cores 1 \
        --software-deployment-method conda \
        --conda-prefix "$REPO_ROOT/.snakemake/conda" \
        --conda-frontend conda \
        --list-conda-envs \
        > "$ENV_LIST"
fi
mapfile -t READ_QC_TARGETS < <(
    awk -F '\t' 'NR > 1 && $1 ~ /(^|\/)read_qc\.yaml$/ {print $3}' "$ENV_LIST"
)
if (( ${#READ_QC_TARGETS[@]} != 1 )); then
    echo "Expected one resolved read_qc environment, found ${#READ_QC_TARGETS[@]}." >&2
    exit 2
fi
READ_QC_TARGET="${READ_QC_TARGETS[0]}"
if [[ "$READ_QC_TARGET" != /* ]]; then
    READ_QC_TARGET="$REPO_ROOT/$READ_QC_TARGET"
fi
if ! ACTIVATION_REPORT="$(
    bash -c '
        set -euo pipefail
        source "$1" "$2"
        printf "%s\t%s\t%s\n" \
            "$CONDA_PREFIX" "$(command -v fastqc)" "$(command -v cutadapt)"
    ' bash "$REPO_ROOT/.venv/bin/activate" "$READ_QC_TARGET"
)"; then
    echo "The read_qc environment cannot be activated: $READ_QC_TARGET" >&2
    exit 2
fi
IFS=$'\t' read -r ACTIVATED_PREFIX FASTQC_PATH CUTADAPT_PATH <<< "$ACTIVATION_REPORT"
if [[ "$ACTIVATED_PREFIX" != "$READ_QC_TARGET" ]] \
    || [[ "$FASTQC_PATH" != "$READ_QC_TARGET/bin/fastqc" ]] \
    || [[ "$CUTADAPT_PATH" != "$READ_QC_TARGET/bin/cutadapt" ]]; then
    echo "The read_qc environment activated incorrectly: $ACTIVATION_REPORT" >&2
    exit 2
fi
echo "Verified read_qc activation: fastqc=$FASTQC_PATH cutadapt=$CUTADAPT_PATH"

"$REPO_ROOT/.venv/bin/python" --version
"$REPO_ROOT/.venv/bin/conda" --version
"$REPO_ROOT/.venv/bin/snakemake" --version
echo "[$(date --iso-8601=seconds)] complete environment installed under $REPO_ROOT"
