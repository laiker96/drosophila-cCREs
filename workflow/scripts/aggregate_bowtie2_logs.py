#!/usr/bin/env python3
"""Aggregate restartable chunk summaries into one Bowtie2-compatible lane log."""

import argparse
from pathlib import Path
import re


TOTAL_RE = re.compile(r"^\s*(\d+) reads; of these:$", re.MULTILINE)


def _one_count(pattern: str, text: str, path: Path) -> int:
    matches = re.findall(pattern, text, flags=re.MULTILINE)
    if len(matches) != 1:
        raise ValueError(f"{path}: expected one Bowtie2 count for {pattern!r}")
    return int(matches[0])


def bowtie2_counts(path: Path, layout: str) -> tuple[int, int, int, int]:
    text = path.read_text(encoding="utf-8")
    totals = TOTAL_RE.findall(text)
    if len(totals) != 1:
        raise ValueError(f"{path}: expected one Bowtie2 read total")
    concordant = "concordantly " if layout == "paired" else ""
    zero = _one_count(
        rf"^\s*(\d+) \([^\n]+\) aligned {concordant}0 times$", text, path
    )
    once = _one_count(
        rf"^\s*(\d+) \([^\n]+\) aligned {concordant}exactly 1 time$", text, path
    )
    multiple = _one_count(
        rf"^\s*(\d+) \([^\n]+\) aligned {concordant}>1 times$", text, path
    )
    total = int(totals[0])
    if zero + once + multiple != total:
        raise ValueError(f"{path}: Bowtie2 category counts do not sum to {total}")
    return total, zero, once, multiple


def _percentage(value: int, total: int) -> float:
    return 100.0 * value / total if total else 0.0


def aggregate(paths: list[Path], layout: str) -> str:
    total = zero = once = multiple = 0
    for path in paths:
        counts = bowtie2_counts(path, layout)
        total += counts[0]
        zero += counts[1]
        once += counts[2]
        multiple += counts[3]
    unit = "paired" if layout == "paired" else "unpaired"
    concordant = "concordantly " if layout == "paired" else ""
    aligned = once + multiple
    return "\n".join(
        [
            f"{total} reads; of these:",
            f"  {total} (100.00%) were {unit}; of these:",
            f"    {zero} ({_percentage(zero, total):.2f}%) "
            f"aligned {concordant}0 times",
            f"    {once} ({_percentage(once, total):.2f}%) "
            f"aligned {concordant}exactly 1 time",
            f"    {multiple} ({_percentage(multiple, total):.2f}%) "
            f"aligned {concordant}>1 times",
            f"{_percentage(aligned, total):.2f}% overall alignment rate",
            f"restartable chunks: {len(paths)}",
            "",
        ]
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--layout", choices=("single", "paired"), required=True)
    parser.add_argument("bowtie2_log", nargs="+", type=Path)
    args = parser.parse_args()
    print(aggregate(args.bowtie2_log, args.layout), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
