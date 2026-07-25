"""Validate one immutable external final BAM and write a deterministic receipt."""

import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def chrom_sizes(path: Path) -> list[tuple[str, int]]:
    result = []
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 2:
                raise ValueError(f"{path}:{line_number}: expected two columns")
            result.append((fields[0], int(fields[1])))
    if not result or len(result) != len({name for name, _size in result}):
        raise ValueError(f"Invalid chromosome sizes: {path}")
    return result


def bam_header(path: Path) -> tuple[str, list[tuple[str, int]]]:
    completed = subprocess.run(
        ["samtools", "view", "-H", str(path)],
        check=True,
        capture_output=True,
        text=True,
    )
    sort_order = ""
    sequences = []
    for line in completed.stdout.splitlines():
        fields = line.split("\t")
        if fields[0] == "@HD":
            sort_order = next(
                (field.removeprefix("SO:") for field in fields[1:] if field.startswith("SO:")),
                "",
            )
        elif fields[0] == "@SQ":
            name = next(
                field.removeprefix("SN:")
                for field in fields[1:]
                if field.startswith("SN:")
            )
            length = int(
                next(
                    field.removeprefix("LN:")
                    for field in fields[1:]
                    if field.startswith("LN:")
                )
            )
            sequences.append((name, length))
    return sort_order, sequences


bam = Path(str(snakemake.input.bam))
bai = Path(str(snakemake.input.bai))
sizes = Path(str(snakemake.input.chrom_sizes))
expected = dict(snakemake.params.expected)
output = Path(str(snakemake.output.validation))
log = Path(str(snakemake.log[0]))
output.parent.mkdir(parents=True, exist_ok=True)
log.parent.mkdir(parents=True, exist_ok=True)

with log.open("w", encoding="utf-8") as log_handle:
    subprocess.run(
        ["samtools", "quickcheck", "-v", str(bam)],
        check=True,
        stdout=log_handle,
        stderr=subprocess.STDOUT,
    )
    subprocess.run(
        ["samtools", "idxstats", "-X", str(bam), str(bai)],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=log_handle,
    )

observed_bam_sha256 = sha256_file(bam)
observed_bai_sha256 = sha256_file(bai)
if observed_bam_sha256 != expected["bam_sha256"]:
    raise ValueError(f"BAM SHA-256 mismatch for {bam}")
if observed_bai_sha256 != expected["bai_sha256"]:
    raise ValueError(f"BAI SHA-256 mismatch for {bai}")

sort_order, header_sequences = bam_header(bam)
reference_sequences = chrom_sizes(sizes)
if sort_order != "coordinate":
    raise ValueError(f"External BAM is not coordinate sorted: {bam}")
if header_sequences != reference_sequences:
    raise ValueError(f"External BAM sequence dictionary differs from {sizes}: {bam}")

receipt = {
    "status": "ok",
    "library_id": str(
        getattr(
            snakemake.wildcards,
            "sample",
            getattr(snakemake.wildcards, "library", ""),
        )
    ),
    "bam": str(bam.resolve()),
    "bai": str(bai.resolve()),
    "bam_sha256": observed_bam_sha256,
    "bai_sha256": observed_bai_sha256,
    "genome": expected["genome"],
    "filtering_contract": expected["filtering_contract"],
    "qc_status": expected["qc_status"],
    "sort_order": sort_order,
    "sequence_count": len(header_sequences),
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
