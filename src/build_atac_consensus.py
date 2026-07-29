#!/usr/bin/env python3
"""Retain pooled ATAC peaks supported by biological-replicate peak calls."""

from __future__ import annotations

import argparse
from pathlib import Path

from short_read_processing.consensus import build_condition_consensus


def replicate(value: str) -> tuple[str, Path]:
    sample, separator, path = value.partition("=")
    if not separator or not sample or not path:
        raise argparse.ArgumentTypeError("use SAMPLE_ID=PEAKS_BED")
    return sample, Path(path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--condition-id", required=True)
    parser.add_argument("--pooled-peaks", type=Path, required=True)
    parser.add_argument("--replicate", action="append", type=replicate, required=True)
    parser.add_argument("--minimum-replicates", type=int, default=2)
    parser.add_argument("--overlap-fraction", type=float, default=0.5)
    parser.add_argument("--output-bed", type=Path, required=True)
    parser.add_argument("--support-tsv", type=Path, required=True)
    parser.add_argument("--stats-json", type=Path, required=True)
    args = parser.parse_args()
    replicates = dict(args.replicate)
    if len(replicates) != len(args.replicate):
        parser.error("--replicate sample IDs must be unique")
    metrics = build_condition_consensus(
        condition_id=args.condition_id,
        pooled_peaks=args.pooled_peaks,
        replicate_peaks=replicates,
        output_bed=args.output_bed,
        support_tsv=args.support_tsv,
        stats_json=args.stats_json,
        minimum_replicates=args.minimum_replicates,
        overlap_fraction=args.overlap_fraction,
    )
    print(
        f"pooled={metrics['pooled_peaks']} "
        f"retained={metrics['retained_replicate_supported_peaks']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
