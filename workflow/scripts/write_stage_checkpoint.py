"""Write one deterministic, checksummed logical-stage checkpoint."""

import hashlib
import json
import os
from pathlib import Path
import tempfile


def sha256_file(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def relative_path(path, base):
    return Path(os.path.relpath(path.resolve(), base.resolve())).as_posix()


output = Path(str(snakemake.output.manifest))
log = Path(str(snakemake.log[0]))
output.parent.mkdir(parents=True, exist_ok=True)
log.parent.mkdir(parents=True, exist_ok=True)
artifact_paths = {
    str(label): Path(str(path))
    for label, path in dict(snakemake.params.artifacts).items()
}
missing = [path for path in artifact_paths.values() if not path.is_file()]
if missing:
    raise ValueError("Missing checkpoint artifacts: " + ", ".join(map(str, missing)))

payload = {
    "schema_version": 1,
    "stage": str(snakemake.wildcards.stage),
    "source_project": str(snakemake.params.source_project),
    "source_run_id": str(snakemake.params.source_run_id),
    "semantic_sha256": str(snakemake.params.semantic_sha256),
    "parameters": dict(snakemake.params.parameters),
    "artifacts": {
        label: {
            "path": relative_path(path, output.parent),
            "sha256": sha256_file(path),
        }
        for label, path in sorted(artifact_paths.items())
    },
}
content = json.dumps(payload, indent=2, sort_keys=True) + "\n"
if not output.is_file() or output.read_text(encoding="utf-8") != content:
    temporary_name = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            prefix=f".{output.name}.",
            dir=output.parent,
            delete=False,
        ) as handle:
            temporary_name = handle.name
            handle.write(content)
        os.replace(temporary_name, output)
    finally:
        if temporary_name and os.path.exists(temporary_name):
            os.unlink(temporary_name)
log.write_text(
    f"Exported {len(artifact_paths)} artifact(s) for {snakemake.wildcards.stage}\n",
    encoding="utf-8",
)
