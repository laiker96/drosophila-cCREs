#!/usr/bin/env python3
"""Write summit-centered H3K27ac counting windows for a master DHS set."""

from __future__ import annotations

import argparse
from pathlib import Path

from short_read_processing.regulatory_elements import write_window_definitions


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--master-bed", type=Path, required=True)
    parser.add_argument("--summit-bed", type=Path, required=True)
    parser.add_argument("--chrom-sizes", type=Path, required=True)
    parser.add_argument("--output-table", type=Path, required=True)
    parser.add_argument("--output-bed", type=Path, required=True)
    args = parser.parse_args()
    metrics = write_window_definitions(
        master_bed=args.master_bed,
        summit_bed=args.summit_bed,
        chrom_sizes=args.chrom_sizes,
        output_table=args.output_table,
        output_bed=args.output_bed,
    )
    print(
        f"windows={metrics['window_count']} "
        f"countable={metrics['countable_window_count']} "
        f"zero_width={metrics['zero_width_window_count']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
