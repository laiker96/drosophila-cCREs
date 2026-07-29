"""Validate one immutable external replicate-peak artifact."""

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


peak = Path(str(snakemake.input.peak))
expected = dict(snakemake.params.expected)
observed = sha256_file(peak)
if observed != expected["sha256"]:
    raise ValueError(f"QC replicate-peak SHA-256 mismatch: {peak}")
with peak.open(encoding="utf-8") as handle:
    records = sum(bool(line.strip()) for line in handle)

output = Path(str(snakemake.output.validation))
log = Path(str(snakemake.log[0]))
output.parent.mkdir(parents=True, exist_ok=True)
log.parent.mkdir(parents=True, exist_ok=True)
receipt = {
    "status": "ok",
    "library_id": str(snakemake.wildcards.sample),
    "method": expected["method"],
    "peak_sha256": observed,
    "record_count": records,
}
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
        json.dump(receipt, handle, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(temporary_name, output)
finally:
    if temporary_name and os.path.exists(temporary_name):
        os.unlink(temporary_name)
log.write_text(f"Validated external QC peak {peak}\n", encoding="utf-8")
