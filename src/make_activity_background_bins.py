#!/usr/bin/env python3
"""Create deterministic autosomal background bins for activity normalization."""

from __future__ import annotations

import argparse
from pathlib import Path

from short_read_processing.activity_tmm import (
    BACKGROUND_BIN_WIDTH,
    write_background_bins,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chrom-sizes", type=Path, required=True)
    parser.add_argument("--autosomes", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bin-width", type=int, default=BACKGROUND_BIN_WIDTH)
    args = parser.parse_args()
    count = write_background_bins(
        chrom_sizes_path=args.chrom_sizes,
        autosomes_path=args.autosomes,
        output_path=args.output,
        bin_width=args.bin_width,
    )
    print(f"Wrote {count} background bins")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
