"""Validate an immutable master-DHS bundle and write a deterministic receipt."""

import csv
import hashlib
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


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def read_chrom_sizes(path: Path) -> dict[str, int]:
    result = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            chrom, size_text = line.rstrip("\n").split("\t")
            if chrom in result or int(size_text) < 1:
                raise ValueError(f"{path}:{line_number}: invalid chromosome")
            result[chrom] = int(size_text)
    if not result:
        raise ValueError(f"Empty chromosome sizes: {path}")
    return result


def read_bed(
    path: Path,
    *,
    sizes: dict[str, int],
    one_base: bool,
) -> list[tuple[str, int, int, str]]:
    order = {chrom: index for index, chrom in enumerate(sizes)}
    rows = []
    previous = None
    identifiers = set()
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 6:
                raise ValueError(f"{path}:{line_number}: expected strict BED6")
            chrom, start_text, end_text, identifier, score, strand = fields
            start, end = int(start_text), int(end_text)
            if (
                chrom not in sizes
                or start < 0
                or end <= start
                or end > sizes[chrom]
                or (one_base and end - start != 1)
                or identifier in identifiers
                or not score.isdigit()
                or not 0 <= int(score) <= 1000
                or strand not in {"+", "-", "."}
            ):
                raise ValueError(f"{path}:{line_number}: invalid BED6 record")
            key = (order[chrom], start, end, identifier)
            if previous is not None and key <= previous:
                raise ValueError(f"{path}:{line_number}: records are not strictly sorted")
            if rows and rows[-1][0] == chrom and start < rows[-1][2] and not one_base:
                raise ValueError(f"{path}:{line_number}: master intervals overlap")
            identifiers.add(identifier)
            rows.append((chrom, start, end, identifier))
            previous = key
    if not rows:
        raise ValueError(f"Master BED is empty: {path}")
    return rows


expected = dict(snakemake.params.expected)
paths = {field: Path(str(getattr(snakemake.input, field))) for field in ARTIFACT_FIELDS}
sizes = read_chrom_sizes(Path(str(snakemake.input.chrom_sizes)))
observed_hashes = {}
for field, path in paths.items():
    observed = sha256_file(path)
    if observed != expected[f"{field}_sha256"]:
        raise ValueError(f"SHA-256 mismatch for {field}: {path}")
    observed_hashes[f"{field}_sha256"] = observed

masters = read_bed(paths["master_bed"], sizes=sizes, one_base=False)
summits = read_bed(paths["summits_bed"], sizes=sizes, one_base=True)
master_by_id = {identifier: (chrom, start, end) for chrom, start, end, identifier in masters}
summit_by_id = {identifier: (chrom, start, end) for chrom, start, end, identifier in summits}
if list(master_by_id) != list(summit_by_id):
    raise ValueError("Master and summit BED identifiers or ordering differ")
for identifier, (chrom, start, end) in master_by_id.items():
    summit_chrom, summit, summit_end = summit_by_id[identifier]
    if summit_chrom != chrom or summit_end != summit + 1 or not start <= summit < end:
        raise ValueError(f"Invalid summit for {identifier}")

with paths["membership_tsv"].open(encoding="utf-8", newline="") as handle:
    membership = list(csv.DictReader(handle, delimiter="\t"))
membership_ids = {row.get("master_dhs_id", "") for row in membership}
if membership_ids != set(master_by_id):
    raise ValueError("Membership table does not cover exactly the master IDs")

with paths["context_matrix_tsv"].open(encoding="utf-8", newline="") as handle:
    matrix = list(csv.DictReader(handle, delimiter="\t"))
if [row.get("master_dhs_id", "") for row in matrix] != list(master_by_id):
    raise ValueError("Context matrix identifiers or ordering differ from the master BED")
for row in matrix:
    identifier = row["master_dhs_id"]
    chrom, start, end = master_by_id[identifier]
    summit = summit_by_id[identifier][1]
    if (
        row.get("chrom") != chrom
        or int(row.get("start", -1)) != start
        or int(row.get("end", -1)) != end
        or int(row.get("summit", -1)) != summit
    ):
        raise ValueError(f"Context matrix coordinates differ for {identifier}")

with paths["stats_json"].open(encoding="utf-8") as handle:
    stats = json.load(handle)
if stats.get("status") != "ok":
    raise ValueError("Master statistics status is not ok")
if int(stats.get("master_dhs_count", -1)) != len(masters):
    raise ValueError("Master statistics count differs from the BED")
if stats.get("method") != expected["method"]:
    raise ValueError("Master statistics method differs from the manifest")

receipt = {
    "status": "ok",
    "genome": expected["genome"],
    "method": expected["method"],
    "source_project": expected["source_project"],
    "source_run_id": expected["source_run_id"],
    "master_dhs_count": len(masters),
    **observed_hashes,
}
output = Path(str(snakemake.output.validation))
log = Path(str(snakemake.log[0]))
output.parent.mkdir(parents=True, exist_ok=True)
log.parent.mkdir(parents=True, exist_ok=True)
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
log.write_text("Validated immutable external master-DHS bundle\n", encoding="utf-8")
