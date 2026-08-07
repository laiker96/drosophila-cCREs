#!/usr/bin/env python3
"""Split one lane into deterministic, mate-preserving FASTQ chunks."""

import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import tempfile
from typing import BinaryIO


BUFFER_SIZE = 1024 * 1024


def _open_fastq(path: Path) -> BinaryIO:
    return gzip.open(path, "rb") if path.suffix == ".gz" else path.open("rb")


def _read_record(handle: BinaryIO, path: Path, record_number: int) -> tuple[bytes, ...] | None:
    header = handle.readline()
    if not header:
        return None
    record = (header, handle.readline(), handle.readline(), handle.readline())
    if any(not line for line in record[1:]):
        raise ValueError(f"{path}: incomplete FASTQ record {record_number}")
    if not record[0].startswith(b"@"):
        raise ValueError(f"{path}: record {record_number} header does not start with '@'")
    if not record[2].startswith(b"+"):
        raise ValueError(f"{path}: record {record_number} separator does not start with '+'")
    return record


def _read_name(record: tuple[bytes, ...]) -> bytes:
    name = record[0][1:].split(maxsplit=1)[0]
    if name.endswith((b"/1", b"/2")):
        name = name[:-2]
    return name


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(BUFFER_SIZE):
            digest.update(block)
    return digest.hexdigest()


def _gzip_writer(path: Path) -> gzip.GzipFile:
    return gzip.GzipFile(
        filename=str(path),
        mode="wb",
        compresslevel=1,
        mtime=0,
    )


def _replace_directory(staging: Path, destination: Path) -> None:
    backup: Path | None = None
    if destination.exists():
        backup = Path(
            tempfile.mkdtemp(prefix=f".{destination.name}.old.", dir=destination.parent)
        )
        backup.rmdir()
        os.replace(destination, backup)
    try:
        os.replace(staging, destination)
    except BaseException:
        if backup is not None and not destination.exists():
            os.replace(backup, destination)
        raise
    if backup is not None:
        shutil.rmtree(backup)


def split_fastqs(
    inputs: list[Path],
    output_dir: Path,
    *,
    unit: str,
    layout: str,
    records_per_chunk: int,
) -> dict[str, object]:
    if records_per_chunk < 1:
        raise ValueError("records_per_chunk must be positive")
    expected_inputs = 2 if layout == "paired" else 1
    if len(inputs) != expected_inputs:
        raise ValueError(f"{layout} layout requires {expected_inputs} FASTQ input(s)")

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(
        tempfile.mkdtemp(prefix=f".{output_dir.name}.partial.", dir=output_dir.parent)
    )
    chunks: list[dict[str, object]] = []
    handles = [_open_fastq(path) for path in inputs]
    writers: list[gzip.GzipFile] = []
    chunk_records = 0
    total_records = 0

    def close_chunk() -> None:
        nonlocal writers
        for writer in writers:
            writer.close()
        writers = []

    def open_chunk() -> None:
        nonlocal chunk_records, writers
        chunk_id = f"chunk{len(chunks) + 1:06d}"
        suffixes = ("R1", "R2") if layout == "paired" else ("SE",)
        names = [f"{chunk_id}_{suffix}.fastq.gz" for suffix in suffixes]
        writers = [_gzip_writer(staging / name) for name in names]
        chunks.append(
            {
                "id": chunk_id,
                "r1": names[0],
                "r2": names[1] if layout == "paired" else None,
                "records": 0,
            }
        )
        chunk_records = 0

    try:
        open_chunk()
        while True:
            record_number = total_records + 1
            records = [
                _read_record(handle, path, record_number)
                for handle, path in zip(handles, inputs, strict=True)
            ]
            if all(record is None for record in records):
                break
            if any(record is None for record in records):
                raise ValueError("Paired FASTQs contain different record counts")
            complete_records = [record for record in records if record is not None]
            if layout == "paired" and _read_name(complete_records[0]) != _read_name(
                complete_records[1]
            ):
                raise ValueError(
                    f"Paired FASTQ identifiers differ at record {record_number}"
                )
            if chunk_records == records_per_chunk:
                close_chunk()
                open_chunk()
            for writer, record in zip(writers, complete_records, strict=True):
                writer.writelines(record)
            chunk_records += 1
            total_records += 1
            chunks[-1]["records"] = chunk_records
        close_chunk()
        manifest: dict[str, object] = {
            "schema_version": 1,
            "unit": unit,
            "layout": layout,
            "records_per_chunk": records_per_chunk,
            "total_records": total_records,
            "inputs": [
                {
                    "path": str(path.resolve()),
                    "bytes": path.stat().st_size,
                    "sha256": _sha256(path),
                }
                for path in inputs
            ],
            "chunks": chunks,
        }
        manifest_path = staging / "manifest.json"
        manifest_path.write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        _replace_directory(staging, output_dir)
        print(
            f"split_fastq\tunit={unit}\trecords={total_records}\t"
            f"chunks={len(chunks)}\trecords_per_chunk={records_per_chunk}"
        )
        return manifest
    except BaseException:
        close_chunk()
        shutil.rmtree(staging, ignore_errors=True)
        raise
    finally:
        for handle in handles:
            handle.close()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--unit", required=True)
    parser.add_argument("--layout", choices=("single", "paired"), required=True)
    parser.add_argument("--records-per-chunk", type=int, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("fastq", nargs="+", type=Path)
    args = parser.parse_args()
    split_fastqs(
        args.fastq,
        args.output_dir,
        unit=args.unit,
        layout=args.layout,
        records_per_chunk=args.records_per_chunk,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
