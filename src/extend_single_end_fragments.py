#!/usr/bin/env python3
"""Extend BED6 single-end reads to fixed, strand-aware fragment intervals."""

from __future__ import annotations

import argparse
from collections.abc import Iterable, Iterator
from pathlib import Path
import sys


def read_chrom_sizes(path: Path) -> dict[str, int]:
    sizes: dict[str, int] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 2:
                raise ValueError(f"{path}:{line_number}: expected two fields")
            chrom, size_text = fields
            size = int(size_text)
            if not chrom or size < 1 or chrom in sizes:
                raise ValueError(f"{path}:{line_number}: invalid chromosome size")
            sizes[chrom] = size
    if not sizes:
        raise ValueError(f"Chromosome sizes file is empty: {path}")
    return sizes


def extend_intervals(
    lines: Iterable[str],
    *,
    chrom_sizes: dict[str, int],
    fragment_length: int,
) -> Iterator[str]:
    """Yield clipped BED6 fragments anchored at each read's five-prime end."""

    if fragment_length < 1:
        raise ValueError("Fragment length must be positive")
    for line_number, line in enumerate(lines, start=1):
        if not line.strip():
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 6:
            raise ValueError(f"BED line {line_number}: expected six fields")
        chrom, start_text, end_text, name, score, strand = fields
        if chrom not in chrom_sizes:
            raise ValueError(f"BED line {line_number}: unknown chromosome {chrom!r}")
        try:
            start, end = int(start_text), int(end_text)
        except ValueError as error:
            raise ValueError(
                f"BED line {line_number}: start and end must be integers"
            ) from error
        chrom_size = chrom_sizes[chrom]
        if start < 0 or end <= start or end > chrom_size:
            raise ValueError(f"BED line {line_number}: invalid interval")
        if end - start > fragment_length:
            raise ValueError(
                f"BED line {line_number}: read is longer than fragment estimate"
            )
        if strand == "+":
            fragment_start = start
            fragment_end = min(chrom_size, start + fragment_length)
        elif strand == "-":
            fragment_start = max(0, end - fragment_length)
            fragment_end = end
        else:
            raise ValueError(f"BED line {line_number}: invalid strand {strand!r}")
        yield (
            f"{chrom}\t{fragment_start}\t{fragment_end}\t"
            f"{name}\t{score}\t{strand}\n"
        )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chrom-sizes", required=True, type=Path)
    parser.add_argument("--fragment-length", required=True, type=int)
    parser.add_argument("--count-output", required=True, type=Path)
    args = parser.parse_args()

    count = 0
    for interval in extend_intervals(
        sys.stdin,
        chrom_sizes=read_chrom_sizes(args.chrom_sizes),
        fragment_length=args.fragment_length,
    ):
        sys.stdout.write(interval)
        count += 1
    args.count_output.write_text(f"{count}\n", encoding="utf-8")
    print(
        f"Extended {count} reads to {args.fragment_length}-bp fragments",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
