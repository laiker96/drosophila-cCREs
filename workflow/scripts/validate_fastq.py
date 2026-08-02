#!/usr/bin/env python3
"""Strictly validate FASTQ structure and character ranges."""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path


SEQUENCE_CHARACTERS = frozenset(b"ACGTURYKMSWBDHVNacgturykmswbdhvn.*-")


def _strip_newline(line: bytes) -> bytes:
    if line.endswith(b"\n"):
        line = line[:-1]
    if line.endswith(b"\r"):
        line = line[:-1]
    return line


def validate_fastq(path: Path) -> int:
    opener = gzip.open if path.suffix == ".gz" else open
    records = 0
    with opener(path, "rb") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            records += 1
            sequence = handle.readline()
            separator = handle.readline()
            quality = handle.readline()
            if not sequence or not separator or not quality:
                raise ValueError(f"{path}: incomplete FASTQ record {records}")

            header = _strip_newline(header)
            sequence = _strip_newline(sequence)
            separator = _strip_newline(separator)
            quality = _strip_newline(quality)
            if not header.startswith(b"@"):
                raise ValueError(f"{path}: record {records} header does not start with '@'")
            if not separator.startswith(b"+"):
                raise ValueError(f"{path}: record {records} separator does not start with '+'")
            if not sequence:
                raise ValueError(f"{path}: record {records} has an empty sequence")
            invalid_sequence = set(sequence).difference(SEQUENCE_CHARACTERS)
            if invalid_sequence:
                values = ",".join(str(value) for value in sorted(invalid_sequence))
                raise ValueError(
                    f"{path}: record {records} contains invalid sequence byte(s): {values}"
                )
            if len(sequence) != len(quality):
                raise ValueError(
                    f"{path}: record {records} has sequence length {len(sequence)} "
                    f"but quality length {len(quality)}"
                )
            invalid_quality = sorted({value for value in quality if not 33 <= value <= 126})
            if invalid_quality:
                values = ",".join(str(value) for value in invalid_quality)
                raise ValueError(
                    f"{path}: record {records} contains invalid quality byte(s): {values}"
                )
    return records


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--paired", action="store_true")
    parser.add_argument("fastq", nargs="+", type=Path)
    args = parser.parse_args()
    if args.paired and len(args.fastq) != 2:
        parser.error("--paired requires exactly two FASTQ files")

    counts = []
    for path in args.fastq:
        count = validate_fastq(path)
        counts.append(count)
        print(f"validated\t{path}\trecords={count}")
    if args.paired and counts[0] != counts[1]:
        raise ValueError(
            f"paired FASTQs contain different record counts: {counts[0]} and {counts[1]}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
