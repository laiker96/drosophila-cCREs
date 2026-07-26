#!/usr/bin/env python3
"""Clip a bedGraph to declared chromosome bounds using an atomic output."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import tempfile


def read_chrom_sizes(path: Path) -> dict[str, int]:
    sizes: dict[str, int] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError(
                    f"{path}:{line_number}: expected chromosome and size"
                )
            chrom, size_text = fields[:2]
            size = int(size_text)
            if not chrom or size < 1 or chrom in sizes:
                raise ValueError(f"{path}:{line_number}: invalid chromosome size")
            sizes[chrom] = size
    if not sizes:
        raise ValueError(f"Chromosome sizes file is empty: {path}")
    return sizes


def clip_bedgraph(source: Path, destination: Path, sizes_path: Path) -> dict[str, int]:
    """Clip intervals to chromosome bounds and discard intervals outside them."""

    if source.resolve() == destination.resolve():
        raise ValueError("Input and output bedGraph paths must differ")
    sizes = read_chrom_sizes(sizes_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary_name: str | None = None
    stats = {
        "records": 0,
        "written": 0,
        "left_clipped": 0,
        "right_clipped": 0,
        "discarded_outside": 0,
        "maximum_right_clip": 0,
    }
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            prefix=f".{destination.name}.",
            dir=destination.parent,
            delete=False,
        ) as output:
            temporary_name = output.name
            with source.open(encoding="utf-8") as handle:
                for line_number, line in enumerate(handle, start=1):
                    if not line.strip() or line.startswith(("track", "browser", "#")):
                        continue
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) < 4:
                        raise ValueError(
                            f"{source}:{line_number}: expected four bedGraph fields"
                        )
                    chrom, start_text, end_text, value_text = fields[:4]
                    if chrom not in sizes:
                        raise ValueError(
                            f"{source}:{line_number}: unknown chromosome {chrom}"
                        )
                    try:
                        start, end = int(start_text), int(end_text)
                        float(value_text)
                    except ValueError as error:
                        raise ValueError(
                            f"{source}:{line_number}: invalid bedGraph value"
                        ) from error
                    if end <= start:
                        raise ValueError(
                            f"{source}:{line_number}: non-positive interval width"
                        )

                    stats["records"] += 1
                    chrom_size = sizes[chrom]
                    clipped_start = max(0, start)
                    clipped_end = min(chrom_size, end)
                    if clipped_start != start:
                        stats["left_clipped"] += 1
                    if clipped_end != end:
                        stats["right_clipped"] += 1
                        stats["maximum_right_clip"] = max(
                            stats["maximum_right_clip"], end - chrom_size
                        )
                    if clipped_end <= clipped_start:
                        stats["discarded_outside"] += 1
                        continue
                    output.write(
                        f"{chrom}\t{clipped_start}\t{clipped_end}\t{value_text}\n"
                    )
                    stats["written"] += 1
        os.replace(temporary_name, destination)
        temporary_name = None
    finally:
        if temporary_name:
            Path(temporary_name).unlink(missing_ok=True)
    return stats


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bedgraph", type=Path, required=True)
    parser.add_argument("--chrom-sizes", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    stats = clip_bedgraph(args.bedgraph, args.output, args.chrom_sizes)
    print(json.dumps(stats, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
