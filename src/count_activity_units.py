#!/usr/bin/env python3
"""Count assay-unit BED records over a master DHS registry."""

from __future__ import annotations

import argparse
from pathlib import Path

from short_read_processing.activity import (
    count_units,
    read_master_elements,
    write_library_signal,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--master-bed", type=Path, required=True)
    parser.add_argument("--summit-bed", type=Path, required=True)
    parser.add_argument("--chrom-sizes", type=Path, required=True)
    parser.add_argument("--units-bed", type=Path, required=True)
    parser.add_argument("--total-units", type=Path, required=True)
    parser.add_argument("--library-id", required=True)
    parser.add_argument("--assay", choices=("atac", "h3k27ac"), required=True)
    parser.add_argument("--cohort", choices=("atlas", "reference"), required=True)
    parser.add_argument("--context", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    args = parser.parse_args()

    elements = read_master_elements(args.master_bed, args.summit_bed)
    chromosome_order = []
    with args.chrom_sizes.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 2 or not fields[0] or int(fields[1]) < 1:
                raise ValueError(
                    f"{args.chrom_sizes}:{line_number}: invalid chromosome sizes"
                )
            chromosome_order.append(fields[0])
    if not chromosome_order or len(chromosome_order) != len(set(chromosome_order)):
        raise ValueError(f"Invalid chromosome order: {args.chrom_sizes}")
    total_units = int(args.total_units.read_text(encoding="utf-8").strip())
    counts = count_units(
        elements,
        args.units_bed,
        expected_total=total_units,
        chromosome_order=chromosome_order,
    )
    write_library_signal(
        elements=elements,
        counts=counts,
        total_units=total_units,
        library_id=args.library_id,
        assay=args.assay,
        cohort=args.cohort,
        context=args.context,
        output=args.output,
        summary=args.summary,
    )
    print(
        f"library={args.library_id} assay={args.assay} "
        f"units={total_units} master_dhs={len(elements)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
