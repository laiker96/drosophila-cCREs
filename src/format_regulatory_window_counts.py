#!/usr/bin/env python3
"""Format one library's raw H3K27ac window-overlap counts."""

from __future__ import annotations

import argparse
from pathlib import Path

from short_read_processing.regulatory_elements import write_window_counts


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--window-table", type=Path, required=True)
    parser.add_argument("--coverage", type=Path, required=True)
    parser.add_argument("--total-units", type=Path, required=True)
    parser.add_argument("--library-id", required=True)
    parser.add_argument("--context", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    metrics = write_window_counts(
        window_table=args.window_table,
        coverage_path=args.coverage,
        total_units_path=args.total_units,
        library_id=args.library_id,
        context=args.context,
        output_path=args.output,
    )
    print(
        f"library={args.library_id} windows={metrics['window_count']} "
        f"overlaps={metrics['overlap_count_sum']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
