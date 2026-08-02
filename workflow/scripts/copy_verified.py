#!/usr/bin/env python3
"""Copy a file and retry until the destination checksum matches the source."""

from __future__ import annotations

import argparse
import hashlib
import os
import time
from pathlib import Path


CHUNK_SIZE = 8 * 1024 * 1024


def _copy_and_hash(source: Path, destination: Path) -> str:
    digest = hashlib.sha256()
    with source.open("rb") as source_handle, destination.open("wb") as destination_handle:
        while chunk := source_handle.read(CHUNK_SIZE):
            digest.update(chunk)
            destination_handle.write(chunk)
        destination_handle.flush()
        os.fsync(destination_handle.fileno())
    return digest.hexdigest()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(CHUNK_SIZE):
            digest.update(chunk)
    return digest.hexdigest()


def copy_verified(source: Path, destination: Path, *, attempts: int = 5) -> int:
    if attempts < 1:
        raise ValueError("attempts must be positive")
    destination.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(1, attempts + 1):
        expected = _copy_and_hash(source, destination)
        actual = _sha256(destination)
        if actual == expected:
            print(
                f"verified_copy\t{source}\t{destination}\t"
                f"bytes={source.stat().st_size}\tattempt={attempt}\tsha256={expected}"
            )
            return attempt
        destination.unlink(missing_ok=True)
        print(
            f"checksum_mismatch\t{source}\t{destination}\t"
            f"attempt={attempt}\texpected={expected}\tactual={actual}",
            flush=True,
        )
        if attempt < attempts:
            time.sleep(min(2 ** (attempt - 1), 30))
    raise OSError(f"could not copy {source} to {destination} after {attempts} attempts")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--attempts", type=int, default=5)
    args = parser.parse_args()
    copy_verified(args.source, args.destination, attempts=args.attempts)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
