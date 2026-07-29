"""Validate ATAC context definitions and retain replicate-supported pooled peaks."""

from __future__ import annotations

import bisect
import json
import os
import re
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, TextIO

from .accessions import AcquisitionError


SAFE_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


@dataclass(frozen=True)
class ConditionSpec:
    condition_id: str
    label: str
    samples: tuple[str, ...]


@dataclass(frozen=True)
class Peak:
    chrom: str
    start: int
    end: int
    name: str
    score: int
    strand: str


def condition_specs(
    values: Any,
    *,
    sample_ids: Iterable[str],
    minimum_replicates: int,
) -> list[ConditionSpec]:
    if minimum_replicates < 2:
        raise AcquisitionError("ATAC minimum_replicates must be at least 2")
    if not isinstance(values, list) or not values:
        raise AcquisitionError("ATAC consensus conditions must be a non-empty list")

    expected = set(sample_ids)
    assigned: dict[str, str] = {}
    seen_conditions: set[str] = set()
    conditions: list[ConditionSpec] = []
    for index, value in enumerate(values, start=1):
        if not isinstance(value, dict):
            raise AcquisitionError(f"ATAC consensus condition {index} must be a mapping")
        condition = str(value.get("id", ""))
        label = str(value.get("label", ""))
        samples = value.get("samples")
        if not SAFE_ID_RE.fullmatch(condition):
            raise AcquisitionError(f"Invalid ATAC condition ID {condition!r}")
        if condition in seen_conditions:
            raise AcquisitionError(f"ATAC condition {condition!r} is defined more than once")
        if not label:
            raise AcquisitionError(f"ATAC condition {condition!r} has a blank label")
        if not isinstance(samples, list) or any(not isinstance(sample, str) for sample in samples):
            raise AcquisitionError(
                f"ATAC condition {condition!r} samples must be a list of library IDs"
            )
        if len(samples) < minimum_replicates:
            raise AcquisitionError(
                f"ATAC condition {condition!r} is below minimum_replicates"
            )
        if len(samples) != len(set(samples)):
            raise AcquisitionError(
                f"ATAC condition {condition!r} assigns a library more than once"
            )
        for sample in samples:
            if sample not in expected:
                raise AcquisitionError(
                    f"ATAC condition {condition!r} names unknown library {sample!r}"
                )
            if sample in assigned:
                raise AcquisitionError(
                    f"ATAC library {sample!r} is assigned to more than one condition"
                )
            assigned[sample] = condition
        seen_conditions.add(condition)
        conditions.append(ConditionSpec(condition, label, tuple(samples)))

    missing_samples = sorted(expected - set(assigned))
    if missing_samples:
        raise AcquisitionError(
            "ATAC consensus conditions do not assign: " + ", ".join(missing_samples)
        )
    return conditions


def read_peaks(path: Path) -> list[Peak]:
    peaks: list[Peak] = []
    previous: tuple[str, int] | None = None
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_number}: expected at least three columns")
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            if start < 0 or end <= start:
                raise ValueError(f"{path}:{line_number}: invalid interval")
            order = (chrom, start)
            if previous is not None and order < previous:
                raise ValueError(f"{path}:{line_number}: peaks are not sorted")
            previous = order
            name = fields[3] if len(fields) > 3 and fields[3] else f"peak_{len(peaks) + 1}"
            try:
                score = int(float(fields[4])) if len(fields) > 4 else 0
            except ValueError:
                score = 0
            strand = fields[5] if len(fields) > 5 and fields[5] in {"+", "-", "."} else "."
            peaks.append(Peak(chrom, start, end, name, min(1000, max(0, score)), strand))
    return peaks


class IntervalUnion:
    def __init__(self, intervals: Iterable[tuple[str, int, int]]):
        grouped: dict[str, list[tuple[int, int]]] = {}
        for chrom, start, end in intervals:
            grouped.setdefault(chrom, []).append((start, end))
        self.intervals: dict[str, list[tuple[int, int]]] = {}
        self.starts: dict[str, list[int]] = {}
        for chrom, values in grouped.items():
            merged: list[tuple[int, int]] = []
            for start, end in sorted(values):
                if not merged or start > merged[-1][1]:
                    merged.append((start, end))
                else:
                    merged[-1] = (merged[-1][0], max(merged[-1][1], end))
            self.intervals[chrom] = merged
            self.starts[chrom] = [start for start, _end in merged]

    def covered_bases(self, chrom: str, start: int, end: int) -> int:
        intervals = self.intervals.get(chrom, [])
        if not intervals:
            return 0
        index = max(0, bisect.bisect_right(self.starts[chrom], start) - 1)
        covered = 0
        while index < len(intervals) and intervals[index][0] < end:
            interval_start, interval_end = intervals[index]
            covered += max(0, min(end, interval_end) - max(start, interval_start))
            index += 1
        return covered


def _atomic_writer(path: Path) -> tuple[TextIO, Path]:
    path.parent.mkdir(parents=True, exist_ok=True)
    handle = tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        prefix=f".{path.name}.",
        dir=path.parent,
        delete=False,
    )
    return handle, Path(handle.name)


def _write_json_atomic(path: Path, value: dict[str, Any]) -> None:
    handle, temporary = _atomic_writer(path)
    try:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")
        handle.close()
        os.replace(temporary, path)
    finally:
        if not handle.closed:
            handle.close()
        temporary.unlink(missing_ok=True)


def build_condition_consensus(
    *,
    condition_id: str,
    pooled_peaks: Path,
    replicate_peaks: dict[str, Path],
    output_bed: Path,
    support_tsv: Path,
    stats_json: Path,
    minimum_replicates: int = 2,
    overlap_fraction: float = 0.5,
) -> dict[str, Any]:
    if not SAFE_ID_RE.fullmatch(condition_id):
        raise ValueError(f"Invalid condition ID: {condition_id!r}")
    if minimum_replicates < 2 or minimum_replicates > len(replicate_peaks):
        raise ValueError("minimum_replicates is inconsistent with replicate inputs")
    if not 0 < overlap_fraction <= 1:
        raise ValueError("overlap_fraction must be in (0, 1]")

    pooled = read_peaks(pooled_peaks)
    replicate_names = list(replicate_peaks)
    indexes = {
        sample: IntervalUnion(
            (peak.chrom, peak.start, peak.end) for peak in read_peaks(path)
        )
        for sample, path in replicate_peaks.items()
    }
    bed_handle, temporary_bed = _atomic_writer(output_bed)
    support_handle, temporary_support = _atomic_writer(support_tsv)
    retained_count = 0
    try:
        support_handle.write(
            "condition_id\tpooled_peak_id\tretained\tsupport_n\treplicate_n"
            "\tsupport_fraction\t"
            + "\t".join(f"{sample}_coverage_fraction" for sample in replicate_names)
            + "\n"
        )
        for peak in pooled:
            length = peak.end - peak.start
            fractions = {
                sample: indexes[sample].covered_bases(peak.chrom, peak.start, peak.end)
                / length
                for sample in replicate_names
            }
            supporting = [
                sample for sample in replicate_names if fractions[sample] >= overlap_fraction
            ]
            support_n = len(supporting)
            retained = support_n >= minimum_replicates
            support_fraction = support_n / len(replicate_names)
            support_handle.write(
                f"{condition_id}\t{peak.name}\t{int(retained)}\t{support_n}\t"
                f"{len(replicate_names)}\t{support_fraction:.6g}\t"
                + "\t".join(f"{fractions[sample]:.6g}" for sample in replicate_names)
                + "\n"
            )
            if not retained:
                continue
            retained_count += 1
            bed_handle.write(
                f"{peak.chrom}\t{peak.start}\t{peak.end}\t{peak.name}\t"
                f"{peak.score}\t{peak.strand}\t{condition_id}\t{support_n}\t"
                f"{len(replicate_names)}\t{support_fraction:.6g}\t"
                f"{','.join(supporting)}\tqpois\n"
            )
        bed_handle.close()
        support_handle.close()
        os.replace(temporary_bed, output_bed)
        os.replace(temporary_support, support_tsv)
    finally:
        if not bed_handle.closed:
            bed_handle.close()
        if not support_handle.closed:
            support_handle.close()
        temporary_bed.unlink(missing_ok=True)
        temporary_support.unlink(missing_ok=True)

    metrics: dict[str, Any] = {
        "status": "ok" if retained_count else "no_replicate_supported_peaks",
        "condition_id": condition_id,
        "peak_method": "qpois",
        "pooled_peaks": len(pooled),
        "retained_replicate_supported_peaks": retained_count,
        "replicate_n": len(replicate_names),
        "minimum_replicates": minimum_replicates,
        "minimum_pooled_peak_coverage_fraction": overlap_fraction,
        "replicates": replicate_names,
    }
    _write_json_atomic(stats_json, metrics)
    return metrics
