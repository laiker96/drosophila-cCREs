"""Write a deterministic reusable manifest for a complete master-DHS bundle."""

import csv
import hashlib
import io
import json
import os
from pathlib import Path
import tempfile


ARTIFACT_FIELDS = (
    "master_bed",
    "summits_bed",
    "membership_tsv",
    "context_matrix_tsv",
    "stats_json",
)


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
metadata = dict(snakemake.params.metadata)
paths = {
    field: Path(str(getattr(snakemake.input, field)))
    for field in ARTIFACT_FIELDS
}
with paths["stats_json"].open(encoding="utf-8") as handle:
    stats = json.load(handle)
if stats.get("method") != metadata["method"]:
    raise ValueError("Master statistics method differs from export metadata")

fieldnames = ["genome", "method", "source_project", "source_run_id"]
fieldnames.extend(
    item
    for field in ARTIFACT_FIELDS
    for item in (field, f"{field}_sha256")
)
row = {
    key: metadata[key]
    for key in ("genome", "method", "source_project", "source_run_id")
}
for field, path in paths.items():
    if not path.is_file():
        raise ValueError(f"Missing master artifact: {path}")
    row[field] = relative_path(path, output.parent)
    row[f"{field}_sha256"] = (
        metadata.get(f"{field}_sha256") or sha256_file(path)
    )

buffer = io.StringIO()
writer = csv.DictWriter(
    buffer,
    fieldnames=fieldnames,
    delimiter="\t",
    lineterminator="\n",
)
writer.writeheader()
writer.writerow(row)
content = buffer.getvalue()
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
log.write_text(f"Exported master-DHS bundle to {output}\n", encoding="utf-8")
