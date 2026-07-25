"""Write a deterministic reusable manifest for completed final BAMs."""

import csv
import hashlib
import io
import os
from pathlib import Path
import tempfile


FIELDS = [
    "library_id",
    "assay",
    "context",
    "role",
    "layout",
    "bam",
    "bai",
    "genome",
    "filtering_contract",
    "bam_sha256",
    "bai_sha256",
    "qc_status",
    "source_project",
    "source_run_id",
    "notes",
]


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
rows = []
for raw in sorted(snakemake.params.rows, key=lambda item: item["library_id"]):
    row = dict(raw)
    bam = Path(str(row["bam"]))
    bai = Path(str(row["bai"]))
    if not bam.is_file() or not bai.is_file():
        raise ValueError(f"Incomplete final BAM pair for {row['library_id']}")
    row["bam_sha256"] = row.get("bam_sha256") or sha256_file(bam)
    row["bai_sha256"] = row.get("bai_sha256") or sha256_file(bai)
    row["bam"] = relative_path(bam, output.parent)
    row["bai"] = relative_path(bai, output.parent)
    rows.append(row)

buffer = io.StringIO()
writer = csv.DictWriter(
    buffer,
    fieldnames=FIELDS,
    delimiter="\t",
    lineterminator="\n",
)
writer.writeheader()
writer.writerows(rows)
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
log.write_text(
    f"Exported {len(rows)} final BAM artifact(s) to {output}\n",
    encoding="utf-8",
)
