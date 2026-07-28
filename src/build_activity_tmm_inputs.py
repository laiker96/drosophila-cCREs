#!/usr/bin/env python3
"""Build a raw-count matrix and metadata for one activity TMM method."""

from __future__ import annotations

import argparse
from pathlib import Path

from short_read_processing.activity_tmm import TMM_METHODS, build_tmm_inputs


def keyed_path(value: str) -> tuple[str, Path]:
    key, separator, path = value.partition("=")
    if not separator or not key or not path:
        raise argparse.ArgumentTypeError("use LIBRARY_ID=PATH")
    return key, Path(path)


def unique_mapping(values: list[tuple[str, Path]], label: str) -> dict[str, Path]:
    result = {}
    for key, path in values:
        if key in result:
            raise ValueError(f"Duplicate {label} library ID: {key}")
        result[key] = path
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--method", choices=TMM_METHODS, required=True)
    parser.add_argument("--signal", action="append", type=keyed_path, required=True)
    parser.add_argument("--background-count", action="append", type=keyed_path, default=[])
    parser.add_argument("--output-counts", type=Path, required=True)
    parser.add_argument("--output-metadata", type=Path, required=True)
    args = parser.parse_args()
    metrics = build_tmm_inputs(
        method=args.method,
        signal_paths=unique_mapping(args.signal, "signal"),
        background_count_paths=(
            unique_mapping(args.background_count, "background-count")
            if args.background_count
            else None
        ),
        output_counts=args.output_counts,
        output_metadata=args.output_metadata,
    )
    print(
        f"Built {metrics['method']} matrix: {metrics['feature_count']} features, "
        f"{metrics['library_count']} libraries"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
