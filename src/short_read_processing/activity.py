"""Count and normalize ATAC/H3K27ac signal over a master-DHS registry."""

from __future__ import annotations

from bisect import bisect_right
from collections import defaultdict
import csv
from dataclasses import dataclass
import gzip
import hashlib
import io
import json
import math
import os
from pathlib import Path
import statistics
import tempfile
from typing import Any, Iterable


LIBRARY_SIGNAL_FIELDS = [
    "master_dhs_id",
    "chrom",
    "start",
    "end",
    "summit",
    "width_bp",
    "library_id",
    "assay",
    "cohort",
    "context",
    "raw_count",
    "total_units",
    "cpm_per_kb",
]
CONTEXT_SIGNAL_FIELDS = [
    "master_dhs_id",
    "chrom",
    "start",
    "end",
    "summit",
    "width_bp",
    "context",
    "assay",
    "library_n",
    "library_ids",
    "raw_count_sum",
    "raw_count_mean",
    "total_units_sum",
    "cpm_per_kb",
    "cpm_per_kb_sd",
]
QNORM_REFERENCE_FIELDS = [
    "master_dhs_id",
    "chrom",
    "start",
    "end",
    "summit",
    "width_bp",
    "assay",
    "reference_context",
    "library_n",
    "library_ids",
    "replicate_cpm_per_kb",
    "reference_cpm_per_kb",
]
FINAL_ACTIVITY_FIELDS = [
    "master_dhs_id",
    "chrom",
    "start",
    "end",
    "summit",
    "width_bp",
    "context",
    "atac_library_n",
    "atac_library_ids",
    "atac_raw_count_sum",
    "atac_raw_count_mean",
    "atac_total_units_sum",
    "atac_cpm_per_kb",
    "atac_cpm_per_kb_sd",
    "atac_qnorm",
    "h3k27ac_library_n",
    "h3k27ac_library_ids",
    "h3k27ac_raw_count_sum",
    "h3k27ac_raw_count_mean",
    "h3k27ac_total_units_sum",
    "h3k27ac_cpm_per_kb",
    "h3k27ac_cpm_per_kb_sd",
    "h3k27ac_qnorm",
    "activity",
]


@dataclass(frozen=True)
class MasterElement:
    chrom: str
    start: int
    end: int
    master_id: str
    summit: int

    @property
    def width(self) -> int:
        return self.end - self.start


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open(encoding="utf-8", newline="")


def read_master_elements(master_bed: Path, summit_bed: Path) -> list[MasterElement]:
    """Read a strict, sorted, non-overlapping BED6 registry and its summits."""

    masters = []
    identifiers = set()
    previous: tuple[str, int, int] | None = None
    completed_chromosomes = set()
    with master_bed.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 6:
                raise ValueError(f"{master_bed}:{line_number}: expected BED6")
            chrom, start_text, end_text, master_id, score, strand = fields
            start, end = int(start_text), int(end_text)
            if (
                not chrom
                or start < 0
                or end <= start
                or not master_id
                or master_id in identifiers
                or not score.isdigit()
                or strand not in {"+", "-", "."}
            ):
                raise ValueError(f"{master_bed}:{line_number}: invalid BED6 record")
            if previous is not None:
                if chrom != previous[0]:
                    completed_chromosomes.add(previous[0])
                if chrom in completed_chromosomes or (
                    chrom == previous[0] and start < previous[1]
                ):
                    raise ValueError(f"{master_bed}: records are not sorted")
                if chrom == previous[0] and start < previous[2]:
                    raise ValueError(f"{master_bed}: master elements overlap")
            identifiers.add(master_id)
            masters.append((chrom, start, end, master_id))
            previous = (chrom, start, end)
    if not masters:
        raise ValueError("Master DHS registry is empty")

    summits = {}
    summit_order = []
    with summit_bed.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 6:
                raise ValueError(f"{summit_bed}:{line_number}: expected BED6")
            chrom, start_text, end_text, master_id = fields[:4]
            start, end = int(start_text), int(end_text)
            if end != start + 1 or master_id in summits:
                raise ValueError(f"{summit_bed}:{line_number}: invalid summit")
            summits[master_id] = (chrom, start)
            summit_order.append(master_id)
    if summit_order != [item[3] for item in masters]:
        raise ValueError("Master and summit BED identifiers or order differ")

    result = []
    for chrom, start, end, master_id in masters:
        summit_chrom, summit = summits[master_id]
        if summit_chrom != chrom or not start <= summit < end:
            raise ValueError(f"Summit falls outside master element {master_id}")
        result.append(MasterElement(chrom, start, end, master_id, summit))
    return result


def cpm_per_kb(raw_count: int, total_units: int, width_bp: int) -> float:
    if raw_count < 0:
        raise ValueError("Raw count must be non-negative")
    if total_units <= 0:
        raise ValueError("Total units must be positive")
    if width_bp <= 0:
        raise ValueError("Element width must be positive")
    return raw_count * 1_000_000_000.0 / (total_units * width_bp)


def count_units(
    elements: list[MasterElement],
    unit_bed: Path,
    *,
    expected_total: int,
    chromosome_order: list[str] | None = None,
) -> list[int]:
    """Count sorted BED intervals over non-overlapping master elements."""

    if expected_total <= 0:
        raise ValueError("Library contains no usable assay units")
    by_chrom: dict[str, list[tuple[int, MasterElement]]] = defaultdict(list)
    for index, element in enumerate(elements):
        by_chrom[element.chrom].append((index, element))
    ordered_chromosomes = chromosome_order or list(by_chrom)
    chrom_order = {
        chrom: index for index, chrom in enumerate(ordered_chromosomes)
    }
    if len(chrom_order) != len(ordered_chromosomes) or any(
        chrom not in chrom_order for chrom in by_chrom
    ):
        raise ValueError("Chromosome order is duplicate or omits a master chromosome")
    ends = {
        chrom: [element.end for _index, element in values]
        for chrom, values in by_chrom.items()
    }
    counts = [0] * len(elements)
    observed_total = 0
    previous: tuple[int, int, int] | None = None
    with _open_text(unit_bed) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"{unit_bed}:{line_number}: expected BED3")
            chrom, start_text, end_text = fields[:3]
            start, end = int(start_text), int(end_text)
            if chrom not in chrom_order or start < 0 or end <= start:
                raise ValueError(f"{unit_bed}:{line_number}: invalid assay unit")
            key = (chrom_order[chrom], start, end)
            if previous is not None and key < previous:
                raise ValueError(f"{unit_bed}:{line_number}: units are not sorted")
            previous = key
            observed_total += 1
            if chrom not in by_chrom:
                continue
            values = by_chrom[chrom]
            position = bisect_right(ends[chrom], start)
            while position < len(values):
                element_index, element = values[position]
                if element.start >= end:
                    break
                counts[element_index] += 1
                position += 1
    if observed_total != expected_total:
        raise ValueError(
            f"Unit BED contains {observed_total} records; expected {expected_total}"
        )
    return counts


def _format_number(value: float | int) -> str:
    if isinstance(value, int):
        return str(value)
    if not math.isfinite(value):
        raise ValueError("Activity tables cannot contain NaN or infinity")
    return format(value, ".17g")


def _tsv_content(fields: list[str], rows: Iterable[dict[str, Any]]) -> str:
    buffer = io.StringIO()
    writer = csv.DictWriter(
        buffer,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    for raw in rows:
        row = {
            field: (
                _format_number(raw[field])
                if isinstance(raw.get(field), (int, float))
                else raw.get(field, "")
            )
            for field in fields
        }
        writer.writerow(row)
    return buffer.getvalue()


def _atomic_bytes_if_changed(path: Path, content: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.is_file() and path.read_bytes() == content:
        return
    temporary_name = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb",
            prefix=f".{path.name}.",
            dir=path.parent,
            delete=False,
        ) as handle:
            temporary_name = handle.name
            handle.write(content)
        os.replace(temporary_name, path)
    finally:
        if temporary_name and os.path.exists(temporary_name):
            os.unlink(temporary_name)


def write_deterministic_gzip(path: Path, content: str) -> None:
    _atomic_bytes_if_changed(
        path,
        gzip.compress(content.encode("utf-8"), compresslevel=6, mtime=0),
    )


def write_json_if_changed(path: Path, value: Any) -> None:
    content = (json.dumps(value, indent=2, sort_keys=True) + "\n").encode("utf-8")
    _atomic_bytes_if_changed(path, content)


def write_library_signal(
    *,
    elements: list[MasterElement],
    counts: list[int],
    total_units: int,
    library_id: str,
    assay: str,
    cohort: str,
    context: str,
    output: Path,
    summary: Path,
) -> None:
    if len(elements) != len(counts):
        raise ValueError("Element and count lengths differ")
    rows = []
    for element, raw_count in zip(elements, counts):
        rows.append(
            {
                "master_dhs_id": element.master_id,
                "chrom": element.chrom,
                "start": element.start,
                "end": element.end,
                "summit": element.summit,
                "width_bp": element.width,
                "library_id": library_id,
                "assay": assay,
                "cohort": cohort,
                "context": context,
                "raw_count": raw_count,
                "total_units": total_units,
                "cpm_per_kb": cpm_per_kb(
                    raw_count,
                    total_units,
                    element.width,
                ),
            }
        )
    write_deterministic_gzip(
        output,
        _tsv_content(LIBRARY_SIGNAL_FIELDS, rows),
    )
    write_json_if_changed(
        summary,
        {
            "status": "ok",
            "library_id": library_id,
            "assay": assay,
            "cohort": cohort,
            "context": context,
            "master_dhs_count": len(elements),
            "total_units": total_units,
            "element_overlap_count": sum(counts),
            "nonzero_element_count": sum(count > 0 for count in counts),
        },
    )


def tie_aware_quantile_normalize(
    target: list[float],
    reference: list[float],
) -> tuple[list[float], int]:
    if len(target) != len(reference) or not target:
        raise ValueError("Target and reference vectors must have equal nonzero length")
    if any(not math.isfinite(value) for value in target + reference):
        raise ValueError("Quantile-normalization vectors must be finite")
    ordered_target = sorted(range(len(target)), key=lambda index: (target[index], index))
    ordered_reference = sorted(reference)
    normalized = [0.0] * len(target)
    tie_groups = 0
    start = 0
    while start < len(target):
        end = start + 1
        value = target[ordered_target[start]]
        while end < len(target) and target[ordered_target[end]] == value:
            end += 1
        replacement = statistics.fmean(ordered_reference[start:end])
        for position in range(start, end):
            normalized[ordered_target[position]] = replacement
        if end - start > 1:
            tie_groups += 1
        start = end
    return normalized, tie_groups


def _read_library_signal(path: Path) -> list[dict[str, Any]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != LIBRARY_SIGNAL_FIELDS:
            raise ValueError(f"{path}: unexpected library-signal columns")
        rows = []
        for row in reader:
            rows.append(
                {
                    **row,
                    "start": int(row["start"]),
                    "end": int(row["end"]),
                    "summit": int(row["summit"]),
                    "width_bp": int(row["width_bp"]),
                    "raw_count": int(row["raw_count"]),
                    "total_units": int(row["total_units"]),
                    "cpm_per_kb": float(row["cpm_per_kb"]),
                }
            )
    if not rows:
        raise ValueError(f"{path}: library signal is empty")
    return rows


def _pearson(left: list[float], right: list[float]) -> float:
    if len(left) != len(right) or not left:
        raise ValueError("Correlation vectors differ")
    left_mean = statistics.fmean(left)
    right_mean = statistics.fmean(right)
    numerator = sum(
        (x - left_mean) * (y - right_mean)
        for x, y in zip(left, right)
    )
    left_ss = sum((value - left_mean) ** 2 for value in left)
    right_ss = sum((value - right_mean) ** 2 for value in right)
    if left_ss == 0 or right_ss == 0:
        return 0.0
    return numerator / math.sqrt(left_ss * right_ss)


def _aggregate_rows(
    rows_by_library: dict[str, list[dict[str, Any]]],
    library_ids: list[str],
) -> list[dict[str, Any]]:
    result = []
    for element_index in range(len(next(iter(rows_by_library.values())))):
        library_rows = [
            rows_by_library[library_id][element_index]
            for library_id in library_ids
        ]
        first = library_rows[0]
        raw_counts = [row["raw_count"] for row in library_rows]
        cpm_values = [row["cpm_per_kb"] for row in library_rows]
        result.append(
            {
                **{
                    field: first[field]
                    for field in (
                        "master_dhs_id",
                        "chrom",
                        "start",
                        "end",
                        "summit",
                        "width_bp",
                    )
                },
                "library_n": len(library_ids),
                "library_ids": ",".join(library_ids),
                "raw_count_sum": sum(raw_counts),
                "raw_count_mean": statistics.fmean(raw_counts),
                "total_units_sum": sum(
                    row["total_units"] for row in library_rows
                ),
                "cpm_per_kb": statistics.fmean(cpm_values),
                "cpm_per_kb_sd": (
                    statistics.stdev(cpm_values)
                    if len(cpm_values) > 1
                    else 0.0
                ),
            }
        )
    return result


def build_activity_outputs(
    *,
    signal_paths: dict[str, Path],
    atlas_contexts: list[str],
    reference_context: str,
    output_library_signal: Path,
    output_context_signal: Path,
    output_reference: Path,
    output_activity: Path,
    output_context_views: dict[str, Path],
    output_metrics: Path,
    output_provenance: Path,
    provenance: dict[str, Any],
) -> dict[str, Any]:
    """Aggregate library signals, quantile-normalize, and write final outputs."""

    rows_by_library = {
        library_id: _read_library_signal(path)
        for library_id, path in sorted(signal_paths.items())
    }
    first_rows = next(iter(rows_by_library.values()))
    master_ids = [row["master_dhs_id"] for row in first_rows]
    for library_id, rows in rows_by_library.items():
        if [row["master_dhs_id"] for row in rows] != master_ids:
            raise ValueError(f"Master IDs differ for library {library_id}")
    metadata = {
        library_id: {
            key: rows[0][key]
            for key in ("assay", "cohort", "context")
        }
        for library_id, rows in rows_by_library.items()
    }

    combined_rows = [
        rows_by_library[library_id][element_index]
        for element_index in range(len(master_ids))
        for library_id in sorted(
            rows_by_library,
            key=lambda item: (
                metadata[item]["cohort"],
                metadata[item]["assay"],
                metadata[item]["context"],
                item,
            ),
        )
    ]
    write_deterministic_gzip(
        output_library_signal,
        _tsv_content(LIBRARY_SIGNAL_FIELDS, combined_rows),
    )

    aggregates: dict[tuple[str, str], list[dict[str, Any]]] = {}
    context_rows = []
    for context in atlas_contexts:
        for assay in ("atac", "h3k27ac"):
            library_ids = sorted(
                library_id
                for library_id, values in metadata.items()
                if values
                == {"assay": assay, "cohort": "atlas", "context": context}
            )
            if not library_ids:
                raise ValueError(f"Atlas context {context!r} lacks {assay}")
            aggregate = _aggregate_rows(rows_by_library, library_ids)
            for row in aggregate:
                row["context"] = context
                row["assay"] = assay
            aggregates[(context, assay)] = aggregate
            context_rows.extend(aggregate)
    write_deterministic_gzip(
        output_context_signal,
        _tsv_content(CONTEXT_SIGNAL_FIELDS, context_rows),
    )

    reference_vectors = {}
    reference_rows = []
    reference_library_ids_by_assay = {}
    for assay in ("atac", "h3k27ac"):
        library_ids = sorted(
            library_id
            for library_id, values in metadata.items()
            if values
            == {
                "assay": assay,
                "cohort": "reference",
                "context": reference_context,
            }
        )
        if len(library_ids) < 2:
            raise ValueError(f"Reference requires two {assay} libraries")
        aggregate = _aggregate_rows(rows_by_library, library_ids)
        reference_vectors[assay] = [row["cpm_per_kb"] for row in aggregate]
        reference_library_ids_by_assay[assay] = library_ids
        for element_index, row in enumerate(aggregate):
            replicate_values = [
                rows_by_library[library_id][element_index]["cpm_per_kb"]
                for library_id in library_ids
            ]
            reference_rows.append(
                {
                    **{
                        field: row[field]
                        for field in (
                            "master_dhs_id",
                            "chrom",
                            "start",
                            "end",
                            "summit",
                            "width_bp",
                        )
                    },
                    "assay": assay,
                    "reference_context": reference_context,
                    "library_n": len(library_ids),
                    "library_ids": ",".join(library_ids),
                    "replicate_cpm_per_kb": ",".join(
                        _format_number(value) for value in replicate_values
                    ),
                    "reference_cpm_per_kb": row["cpm_per_kb"],
                }
            )
    write_deterministic_gzip(
        output_reference,
        _tsv_content(QNORM_REFERENCE_FIELDS, reference_rows),
    )

    qnorm = {}
    tie_counts = {}
    for context in atlas_contexts:
        for assay in ("atac", "h3k27ac"):
            values = [
                row["cpm_per_kb"] for row in aggregates[(context, assay)]
            ]
            normalized, ties = tie_aware_quantile_normalize(
                values,
                reference_vectors[assay],
            )
            qnorm[(context, assay)] = normalized
            tie_counts[f"{context}:{assay}"] = ties

    final_rows = []
    context_final_rows = {context: [] for context in atlas_contexts}
    for element_index, base in enumerate(first_rows):
        for context in atlas_contexts:
            assay_rows = {
                assay: aggregates[(context, assay)][element_index]
                for assay in ("atac", "h3k27ac")
            }
            atac_qnorm = qnorm[(context, "atac")][element_index]
            h3_qnorm = qnorm[(context, "h3k27ac")][element_index]
            row = {
                **{
                    field: base[field]
                    for field in (
                        "master_dhs_id",
                        "chrom",
                        "start",
                        "end",
                        "summit",
                        "width_bp",
                    )
                },
                "context": context,
            }
            for assay, prefix in (("atac", "atac"), ("h3k27ac", "h3k27ac")):
                values = assay_rows[assay]
                for source, suffix in (
                    ("library_n", "library_n"),
                    ("library_ids", "library_ids"),
                    ("raw_count_sum", "raw_count_sum"),
                    ("raw_count_mean", "raw_count_mean"),
                    ("total_units_sum", "total_units_sum"),
                    ("cpm_per_kb", "cpm_per_kb"),
                    ("cpm_per_kb_sd", "cpm_per_kb_sd"),
                ):
                    row[f"{prefix}_{suffix}"] = values[source]
            row["atac_qnorm"] = atac_qnorm
            row["h3k27ac_qnorm"] = h3_qnorm
            row["activity"] = math.sqrt(atac_qnorm * h3_qnorm)
            final_rows.append(row)
            context_final_rows[context].append(row)
    final_content = _tsv_content(FINAL_ACTIVITY_FIELDS, final_rows)
    write_deterministic_gzip(output_activity, final_content)
    for context, path in output_context_views.items():
        write_deterministic_gzip(
            path,
            _tsv_content(FINAL_ACTIVITY_FIELDS, context_final_rows[context]),
        )

    correlations = {}
    for cohort in ("atlas", "reference"):
        contexts = atlas_contexts if cohort == "atlas" else [reference_context]
        for context in contexts:
            for assay in ("atac", "h3k27ac"):
                library_ids = sorted(
                    library_id
                    for library_id, values in metadata.items()
                    if values
                    == {
                        "assay": assay,
                        "cohort": cohort,
                        "context": context,
                    }
                )
                for left_index, left in enumerate(library_ids):
                    for right in library_ids[left_index + 1 :]:
                        correlations[f"{cohort}:{context}:{assay}:{left}:{right}"] = (
                            _pearson(
                                [
                                    row["cpm_per_kb"]
                                    for row in rows_by_library[left]
                                ],
                                [
                                    row["cpm_per_kb"]
                                    for row in rows_by_library[right]
                                ],
                            )
                        )
    widths = [row["width_bp"] for row in first_rows]
    metrics = {
        "status": "ok",
        "schema_version": 1,
        "master_dhs_count": len(master_ids),
        "atlas_context_count": len(atlas_contexts),
        "library_count": len(rows_by_library),
        "width_bp_min": min(widths),
        "width_bp_median": statistics.median(widths),
        "width_bp_max": max(widths),
        "tie_group_counts": tie_counts,
        "replicate_correlations": correlations,
        "library_total_units": {
            library_id: rows[0]["total_units"]
            for library_id, rows in rows_by_library.items()
        },
        "library_zero_fractions": {
            library_id: sum(row["raw_count"] == 0 for row in rows) / len(rows)
            for library_id, rows in rows_by_library.items()
        },
        "reference_library_ids": reference_library_ids_by_assay,
    }
    write_json_if_changed(output_metrics, metrics)
    output_paths = {
        "library_signal": output_library_signal,
        "context_signal": output_context_signal,
        "qnorm_reference": output_reference,
        "activity_table": output_activity,
        "metrics": output_metrics,
        **{
            f"context:{context}": path
            for context, path in output_context_views.items()
        },
    }
    provenance_value = {
        **provenance,
        "schema_version": 1,
        "normalization": "cpm_per_kb_then_tie_aware_reference_qnorm_v1",
        "activity_formula": "sqrt_atac_times_h3k27ac_v1",
        "outputs": {
            name: {
                "path": str(path.resolve()),
                "sha256": sha256_file(path),
            }
            for name, path in sorted(output_paths.items())
        },
    }
    write_json_if_changed(output_provenance, provenance_value)
    return metrics
