"""Descriptive QC for master-DHS activity normalization outputs."""

from __future__ import annotations

from array import array
from collections import defaultdict
import csv
import gzip
import html
import json
import math
import os
from pathlib import Path
import tempfile
from typing import Any, Callable

from .activity import (
    CONTEXT_SIGNAL_FIELDS,
    FINAL_ACTIVITY_FIELDS,
    LIBRARY_SIGNAL_FIELDS,
    QNORM_REFERENCE_FIELDS,
    sha256_file,
)


QUANTILES = (0.0, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1.0)
QUANTILE_FIELDS = tuple(f"q{round(value * 100):03d}" for value in QUANTILES)
CORRELATION_FIELDS = (
    "cohort",
    "context",
    "assay",
    "library_a",
    "library_b",
    "pearson_cpm_per_kb",
    "pearson_log1p_cpm_per_kb",
    "spearman_cpm_per_kb",
)
DISTRIBUTION_FIELDS = (
    "scope",
    "context",
    "assay",
    "stage",
    "element_n",
    "zero_fraction",
    "tie_group_count",
    "largest_tie_group",
    *QUANTILE_FIELDS,
)
PALETTE = (
    "#0072b2",
    "#d55e00",
    "#009e73",
    "#cc79a7",
    "#e69f00",
    "#56b4e9",
    "#000000",
    "#7f3c8d",
    "#6f4e37",
)


def _open_tsv(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open(encoding="utf-8", newline="")


def _finite_float(value: str, *, path: Path, field: str) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{path}: {field} contains a non-finite value")
    return result


def _require_header(reader: csv.DictReader, fields: list[str], path: Path) -> None:
    if reader.fieldnames != fields:
        raise ValueError(f"{path}: unexpected columns")


def _read_library_profiles(
    path: Path,
) -> tuple[list[str], dict[str, dict[str, Any]]]:
    profiles: dict[str, dict[str, Any]] = {}
    master_ids: list[str] = []
    canonical_library: str | None = None
    with _open_tsv(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        _require_header(reader, LIBRARY_SIGNAL_FIELDS, path)
        for row in reader:
            library_id = row["library_id"]
            if canonical_library is None:
                canonical_library = library_id
            profile = profiles.setdefault(
                library_id,
                {
                    "library_id": library_id,
                    "assay": row["assay"],
                    "cohort": row["cohort"],
                    "context": row["context"],
                    "total_units": int(row["total_units"]),
                    "values": array("d"),
                    "zero_count": 0,
                },
            )
            if any(
                profile[key] != row[key]
                for key in ("assay", "cohort", "context")
            ) or profile["total_units"] != int(row["total_units"]):
                raise ValueError(f"{path}: inconsistent metadata for {library_id}")
            index = len(profile["values"])
            master_id = row["master_dhs_id"]
            if library_id == canonical_library:
                if index != len(master_ids):
                    raise ValueError(f"{path}: duplicate or disordered library rows")
                master_ids.append(master_id)
            elif index >= len(master_ids) or master_ids[index] != master_id:
                raise ValueError(f"{path}: master-DHS order differs for {library_id}")
            value = _finite_float(row["cpm_per_kb"], path=path, field="cpm_per_kb")
            if value < 0:
                raise ValueError(f"{path}: cpm_per_kb must be non-negative")
            profile["values"].append(value)
            profile["zero_count"] += int(value == 0)
    if not profiles or not master_ids:
        raise ValueError(f"{path}: empty library signal table")
    if any(len(profile["values"]) != len(master_ids) for profile in profiles.values()):
        raise ValueError(f"{path}: libraries do not cover the same master DHSs")
    return master_ids, profiles


def _read_context_profiles(
    path: Path,
    master_ids: list[str],
) -> dict[tuple[str, str], dict[str, Any]]:
    profiles: dict[tuple[str, str], dict[str, Any]] = {}
    with _open_tsv(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        _require_header(reader, CONTEXT_SIGNAL_FIELDS, path)
        for row in reader:
            key = (row["context"], row["assay"])
            profile = profiles.setdefault(
                key,
                {
                    "context": row["context"],
                    "assay": row["assay"],
                    "library_n": int(row["library_n"]),
                    "library_ids": row["library_ids"],
                    "values": array("d"),
                },
            )
            if (
                profile["library_n"] != int(row["library_n"])
                or profile["library_ids"] != row["library_ids"]
            ):
                raise ValueError(f"{path}: inconsistent aggregate metadata for {key}")
            index = len(profile["values"])
            if index >= len(master_ids) or row["master_dhs_id"] != master_ids[index]:
                raise ValueError(f"{path}: master-DHS order differs for {key}")
            value = _finite_float(row["cpm_per_kb"], path=path, field="cpm_per_kb")
            if value < 0:
                raise ValueError(f"{path}: cpm_per_kb must be non-negative")
            profile["values"].append(value)
    if any(len(profile["values"]) != len(master_ids) for profile in profiles.values()):
        raise ValueError(f"{path}: context/assay profiles are incomplete")
    return profiles


def _read_reference_profiles(
    path: Path,
    master_ids: list[str],
) -> tuple[str, dict[str, array]]:
    profiles: dict[str, array] = {}
    reference_context: str | None = None
    with _open_tsv(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        _require_header(reader, QNORM_REFERENCE_FIELDS, path)
        for row in reader:
            if reference_context is None:
                reference_context = row["reference_context"]
            elif reference_context != row["reference_context"]:
                raise ValueError(f"{path}: multiple reference contexts")
            values = profiles.setdefault(row["assay"], array("d"))
            index = len(values)
            if index >= len(master_ids) or row["master_dhs_id"] != master_ids[index]:
                raise ValueError(f"{path}: reference master-DHS order differs")
            value = _finite_float(
                row["reference_cpm_per_kb"],
                path=path,
                field="reference_cpm_per_kb",
            )
            if value < 0:
                raise ValueError(f"{path}: reference values must be non-negative")
            values.append(value)
    if reference_context is None or any(
        len(values) != len(master_ids) for values in profiles.values()
    ):
        raise ValueError(f"{path}: reference profiles are incomplete")
    return reference_context, profiles


def _read_final_profiles(
    path: Path,
    master_ids: list[str],
) -> dict[str, dict[str, array]]:
    profiles: dict[str, dict[str, array]] = {}
    with _open_tsv(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        _require_header(reader, FINAL_ACTIVITY_FIELDS, path)
        for row in reader:
            context = row["context"]
            profile = profiles.setdefault(
                context,
                {"atac": array("d"), "h3k27ac": array("d"), "activity": array("d")},
            )
            index = len(profile["activity"])
            if index >= len(master_ids) or row["master_dhs_id"] != master_ids[index]:
                raise ValueError(f"{path}: final master-DHS order differs for {context}")
            for field, key in (
                ("atac_qnorm", "atac"),
                ("h3k27ac_qnorm", "h3k27ac"),
                ("activity", "activity"),
            ):
                value = _finite_float(row[field], path=path, field=field)
                if value < 0:
                    raise ValueError(f"{path}: {field} must be non-negative")
                profile[key].append(value)
    if any(
        len(values) != len(master_ids)
        for profile in profiles.values()
        for values in profile.values()
    ):
        raise ValueError(f"{path}: final context profiles are incomplete")
    return profiles


def _quantile(sorted_values: list[float], probability: float) -> float:
    position = probability * (len(sorted_values) - 1)
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return sorted_values[lower]
    fraction = position - lower
    return sorted_values[lower] * (1 - fraction) + sorted_values[upper] * fraction


def _tie_stats(sorted_values: list[float]) -> tuple[int, int]:
    tie_groups = 0
    largest = 1
    start = 0
    while start < len(sorted_values):
        end = start + 1
        while end < len(sorted_values) and sorted_values[end] == sorted_values[start]:
            end += 1
        size = end - start
        if size > 1:
            tie_groups += 1
            largest = max(largest, size)
        start = end
    return tie_groups, largest


def _distribution_record(
    *,
    scope: str,
    context: str,
    assay: str,
    stage: str,
    values: array,
) -> dict[str, Any]:
    ordered = sorted(values)
    tie_groups, largest = _tie_stats(ordered)
    result: dict[str, Any] = {
        "scope": scope,
        "context": context,
        "assay": assay,
        "stage": stage,
        "element_n": len(ordered),
        "zero_fraction": sum(value == 0 for value in ordered) / len(ordered),
        "tie_group_count": tie_groups,
        "largest_tie_group": largest,
    }
    result.update(
        {
            field: _quantile(ordered, probability)
            for field, probability in zip(QUANTILE_FIELDS, QUANTILES)
        }
    )
    return result


def _pearson(
    left: array,
    right: array,
    transform: Callable[[float], float] = lambda value: value,
) -> float:
    if len(left) != len(right) or not left:
        raise ValueError("Correlation vectors must have equal nonzero length")
    left_mean = math.fsum(transform(value) for value in left) / len(left)
    right_mean = math.fsum(transform(value) for value in right) / len(right)
    numerator = math.fsum(
        (transform(x) - left_mean) * (transform(y) - right_mean)
        for x, y in zip(left, right)
    )
    left_ss = math.fsum((transform(value) - left_mean) ** 2 for value in left)
    right_ss = math.fsum((transform(value) - right_mean) ** 2 for value in right)
    if left_ss == 0 or right_ss == 0:
        return 0.0
    return numerator / math.sqrt(left_ss * right_ss)


def _rank(values: array) -> array:
    order = sorted(range(len(values)), key=lambda index: (values[index], index))
    ranks = array("d", [0.0]) * len(values)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        average_rank = (start + end - 1) / 2 + 1
        for position in range(start, end):
            ranks[order[position]] = average_rank
        start = end
    return ranks


def _replicate_qc(
    profiles: dict[str, dict[str, Any]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    groups: dict[tuple[str, str, str], list[str]] = defaultdict(list)
    for library_id, profile in profiles.items():
        groups[(profile["cohort"], profile["context"], profile["assay"])].append(
            library_id
        )
    correlations = []
    unreplicated = []
    for (cohort, context, assay), library_ids in sorted(groups.items()):
        library_ids.sort()
        if len(library_ids) == 1:
            unreplicated.append(
                {
                    "cohort": cohort,
                    "context": context,
                    "assay": assay,
                    "library_ids": library_ids,
                }
            )
            continue
        ranks = {
            library_id: _rank(profiles[library_id]["values"])
            for library_id in library_ids
        }
        for left_index, left_id in enumerate(library_ids):
            for right_id in library_ids[left_index + 1 :]:
                left = profiles[left_id]["values"]
                right = profiles[right_id]["values"]
                correlations.append(
                    {
                        "cohort": cohort,
                        "context": context,
                        "assay": assay,
                        "library_a": left_id,
                        "library_b": right_id,
                        "pearson_cpm_per_kb": _pearson(left, right),
                        "pearson_log1p_cpm_per_kb": _pearson(
                            left, right, math.log1p
                        ),
                        "spearman_cpm_per_kb": _pearson(
                            ranks[left_id], ranks[right_id]
                        ),
                    }
                )
    return correlations, unreplicated


def _sorted_distance(values: array, reference: array) -> dict[str, float | None]:
    left = sorted(values)
    right = sorted(reference)
    if len(left) != len(right) or not left:
        raise ValueError("Distribution vectors must have equal nonzero length")
    differences = [abs(x - y) for x, y in zip(left, right)]
    mean_absolute = math.fsum(differences) / len(differences)
    reference_mean = math.fsum(right) / len(right)
    return {
        "mean_absolute_sorted_difference": mean_absolute,
        "root_mean_square_sorted_difference": math.sqrt(
            math.fsum(value * value for value in differences) / len(differences)
        ),
        "maximum_absolute_sorted_difference": max(differences),
        "mean_absolute_fraction_reference_mean": (
            mean_absolute / reference_mean if reference_mean > 0 else None
        ),
    }


def _format_value(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return format(value, ".17g")
    return str(value)


def _tsv(fields: tuple[str, ...], rows: list[dict[str, Any]]) -> str:
    from io import StringIO

    buffer = StringIO()
    writer = csv.DictWriter(
        buffer,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    for row in rows:
        writer.writerow({field: _format_value(row.get(field)) for field in fields})
    return buffer.getvalue()


def _atomic_write(path: Path, content: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.is_file() and path.read_bytes() == content:
        return
    temporary_name: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb", prefix=f".{path.name}.", dir=path.parent, delete=False
        ) as handle:
            temporary_name = handle.name
            handle.write(content)
        os.replace(temporary_name, path)
    finally:
        if temporary_name and os.path.exists(temporary_name):
            os.unlink(temporary_name)


def _table(headers: list[str], rows: list[list[str]], classes: str = "") -> str:
    head = "".join(f"<th>{html.escape(value)}</th>" for value in headers)
    body = "".join(
        "<tr>" + "".join(f"<td>{value}</td>" for value in row) + "</tr>"
        for row in rows
    )
    return f'<div class="table-wrap"><table class="{classes}"><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>'


def _correlation_cell(value: float) -> str:
    bounded = min(1.0, max(-1.0, value))
    hue = 120 * ((bounded + 1) / 2)
    return (
        f'<span class="corr" style="background:hsl({hue:.1f} 70% 84%)">'
        f"{value:.3f}</span>"
    )


def _quantile_svg(
    distributions: list[dict[str, Any]],
    *,
    assay: str,
    stage: str,
) -> str:
    selected = [
        row
        for row in distributions
        if row["assay"] == assay and row["stage"] in {stage, "reference"}
    ]
    width, height = 760, 310
    left, right, top, bottom = 58, 18, 22, 48
    plot_width = width - left - right
    plot_height = height - top - bottom
    maximum = max(math.log1p(float(row["q100"])) for row in selected) or 1.0
    parts = [
        f'<svg viewBox="0 0 {width} {height}" role="img" aria-label="{assay} {stage} quantile profiles">',
        f'<rect x="{left}" y="{top}" width="{plot_width}" height="{plot_height}" fill="#fff" stroke="#ccd5df"/>',
    ]
    for probability in (0.0, 0.25, 0.5, 0.75, 1.0):
        x = left + probability * plot_width
        parts.append(
            f'<line x1="{x:.1f}" y1="{top}" x2="{x:.1f}" y2="{top + plot_height}" stroke="#edf0f4"/>'
        )
        parts.append(
            f'<text x="{x:.1f}" y="{height - 24}" text-anchor="middle">{probability:.2g}</text>'
        )
    for fraction in (0.0, 0.25, 0.5, 0.75, 1.0):
        y = top + plot_height * (1 - fraction)
        label = math.expm1(maximum * fraction)
        parts.append(
            f'<line x1="{left}" y1="{y:.1f}" x2="{left + plot_width}" y2="{y:.1f}" stroke="#edf0f4"/>'
        )
        parts.append(
            f'<text x="{left - 7}" y="{y + 4:.1f}" text-anchor="end">{label:.2g}</text>'
        )
    atlas_index = 0
    legend = []
    for row in selected:
        reference = row["stage"] == "reference"
        color = "#111827" if reference else PALETTE[atlas_index % len(PALETTE)]
        if not reference:
            atlas_index += 1
        points = []
        for probability, field in zip(QUANTILES, QUANTILE_FIELDS):
            x = left + probability * plot_width
            y = top + plot_height * (
                1 - math.log1p(float(row[field])) / maximum
            )
            points.append(f"{x:.1f},{y:.1f}")
        dash = ' stroke-dasharray="6 4"' if reference else ""
        stroke_width = 2.6 if reference else 1.7
        parts.append(
            f'<polyline points="{" ".join(points)}" fill="none" stroke="{color}" stroke-width="{stroke_width}"{dash}/>'
        )
        label = f'{row["context"]} ({row["stage"]})'
        legend.append((label, color, reference))
    parts.extend(
        [
            f'<text x="{left + plot_width / 2:.1f}" y="{height - 5}" text-anchor="middle">quantile</text>',
            f'<text transform="translate(14 {top + plot_height / 2:.1f}) rotate(-90)" text-anchor="middle">CPM per kb (log1p scale)</text>',
        ]
    )
    legend_x, legend_y = left + 8, top + 14
    for index, (label, color, reference) in enumerate(legend):
        column = index // 6
        row_index = index % 6
        x = legend_x + column * 180
        y = legend_y + row_index * 16
        dash = ' stroke-dasharray="5 3"' if reference else ""
        parts.append(
            f'<line x1="{x}" y1="{y}" x2="{x + 18}" y2="{y}" stroke="{color}" stroke-width="2"{dash}/>'
        )
        parts.append(
            f'<text x="{x + 23}" y="{y + 4}">{html.escape(label)}</text>'
        )
    parts.append("</svg>")
    return "".join(parts)


def _build_html(metrics: dict[str, Any], distributions: list[dict[str, Any]]) -> str:
    correlations = sorted(
        metrics["replicate_correlations"],
        key=lambda row: (row["spearman_cpm_per_kb"], row["library_a"], row["library_b"]),
    )
    correlation_rows = [
        [
            html.escape(row["cohort"]),
            html.escape(row["context"]),
            html.escape(row["assay"]),
            html.escape(row["library_a"]),
            html.escape(row["library_b"]),
            _correlation_cell(row["pearson_cpm_per_kb"]),
            _correlation_cell(row["pearson_log1p_cpm_per_kb"]),
            _correlation_cell(row["spearman_cpm_per_kb"]),
        ]
        for row in correlations
    ]
    unreplicated_rows = [
        [
            html.escape(row["cohort"]),
            html.escape(row["context"]),
            html.escape(row["assay"]),
            html.escape(", ".join(row["library_ids"])),
        ]
        for row in metrics["groups_without_replicates"]
    ]
    library_rows = [
        [
            html.escape(row["cohort"]),
            html.escape(row["context"]),
            html.escape(row["assay"]),
            html.escape(row["library_id"]),
            f'{row["total_units"]:,}',
            f'{row["zero_fraction"]:.3%}',
            f'{row["q050"]:.3g}',
            f'{row["q095"]:.3g}',
        ]
        for row in metrics["library_summaries"]
    ]
    matching_rows = []
    for row in metrics["reference_matching"]:
        pre = row["pre"]
        post = row["post"]
        matching_rows.append(
            [
                html.escape(row["context"]),
                html.escape(row["assay"]),
                f'{row["pre_zero_fraction"]:.3%}',
                f'{row["post_zero_fraction"]:.3%}',
                f'{pre["mean_absolute_sorted_difference"]:.4g}',
                f'{post["mean_absolute_sorted_difference"]:.4g}',
                (
                    f'{post["mean_absolute_fraction_reference_mean"]:.3%}'
                    if post["mean_absolute_fraction_reference_mean"] is not None
                    else "NA"
                ),
                f'{post["maximum_absolute_sorted_difference"]:.4g}',
            ]
        )
    cards = "".join(
        f'<div class="card"><strong>{value}</strong><span>{label}</span></div>'
        for value, label in (
            (f'{metrics["master_dhs_count"]:,}', "master DHSs"),
            (metrics["atlas_context_count"], "atlas contexts"),
            (metrics["library_count"], "libraries"),
            (len(correlations), "replicate pairs"),
        )
    )
    charts = "".join(
        f'<section class="chart"><h3>{assay.upper()} — {stage}</h3>{_quantile_svg(distributions, assay=assay, stage=stage)}</section>'
        for assay in ("atac", "h3k27ac")
        for stage in ("pre_qnorm", "post_qnorm")
    )
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Master-DHS activity QC</title>
<style>
body{{font:14px/1.45 system-ui,-apple-system,sans-serif;color:#172033;background:#f5f7fa;margin:0}}
main{{max-width:1280px;margin:auto;padding:28px}} h1{{margin-bottom:4px}} h2{{margin-top:34px;border-bottom:1px solid #d7dee7;padding-bottom:6px}}
.note{{background:#fff5cc;border-left:4px solid #d69e00;padding:12px 15px;margin:18px 0}} .cards{{display:flex;gap:12px;flex-wrap:wrap}}
.card{{background:white;border:1px solid #d7dee7;border-radius:7px;padding:14px 18px;min-width:150px}} .card strong{{display:block;font-size:23px}} .card span{{color:#58677c}}
.charts{{display:grid;grid-template-columns:repeat(auto-fit,minmax(560px,1fr));gap:16px}} .chart{{background:white;border:1px solid #d7dee7;border-radius:7px;padding:12px}}
svg{{width:100%;height:auto}} svg text{{font-size:10px;fill:#4a5568}} .table-wrap{{overflow:auto;background:white;border:1px solid #d7dee7;border-radius:7px}}
table{{border-collapse:collapse;width:100%;font-size:12px}} th,td{{padding:7px 9px;border-bottom:1px solid #e7ebf0;text-align:right;white-space:nowrap}} th{{position:sticky;top:0;background:#eef2f7}}
th:nth-child(-n+5),td:nth-child(-n+5){{text-align:left}} .corr{{display:inline-block;min-width:48px;padding:2px 4px;border-radius:3px;text-align:right}}
code{{background:#edf1f5;padding:1px 4px;border-radius:3px}} footer{{margin-top:32px;color:#66758a}}
</style></head><body><main>
<h1>Master-DHS activity QC</h1><p>Descriptive QC for <code>cpm_per_kb_then_tie_aware_reference_qnorm_v1</code>.</p>
<div class="note"><strong>Manual review required.</strong> Quantile normalization makes marginal distributions similar by construction and cannot rescue a discordant replicate. Review replicate agreement and pre-normalization sparsity before accepting activity values.</div>
<div class="cards">{cards}</div>
<h2>Reference matching</h2><p>Profiles show fixed quantiles. The black dashed line is the mean {html.escape(metrics["reference_context"])} reference for the matching assay.</p>
<div class="charts">{charts}</div>
{_table(["Context","Assay","Pre zero","Post zero","Pre mean |Δ rank|","Post mean |Δ rank|","Post |Δ| / ref mean","Post max |Δ rank|"], matching_rows)}
<h2>Replicate agreement</h2><p>Correlations use all master DHSs before context averaging or quantile normalization. Rows are ordered by Spearman correlation, weakest first.</p>
{_table(["Cohort","Context","Assay","Library A","Library B","Pearson","Pearson log1p","Spearman"], correlation_rows)}
<h3>Groups without biological replication</h3>
{_table(["Cohort","Context","Assay","Library"], unreplicated_rows) if unreplicated_rows else '<p>None.</p>'}
<h2>Library depth and sparsity</h2>
{_table(["Cohort","Context","Assay","Library","Assay units","Zero DHSs","Median CPM/kb","95th percentile"], library_rows)}
<h2>Interpretation</h2><ul>
<li>Libraries are counted independently and retain raw and CPM-per-kb values.</li>
<li>Biological libraries receive equal weight in each context/assay mean.</li>
<li>ATAC and H3K27ac are normalized separately to the corresponding mean reference profile.</li>
<li>Tied target values receive the mean of their corresponding reference ranks, so post-normalization distributions need not be exactly identical.</li>
<li>No automatic correlation, zero-fraction, or distribution threshold is applied by this report.</li>
</ul>
<footer>Schema version {metrics["schema_version"]}. Deterministic report; input SHA-256 values are recorded in the companion JSON.</footer>
</main></body></html>"""


def build_activity_qc_outputs(
    *,
    library_signal: Path,
    context_signal: Path,
    qnorm_reference: Path,
    activity_table: Path,
    activity_provenance: Path,
    atlas_contexts: list[str],
    reference_context: str,
    output_correlations: Path,
    output_distributions: Path,
    output_metrics: Path,
    output_report: Path,
) -> dict[str, Any]:
    """Validate activity tables and create deterministic descriptive QC outputs."""

    master_ids, libraries = _read_library_profiles(library_signal)
    contexts = _read_context_profiles(context_signal, master_ids)
    observed_reference_context, references = _read_reference_profiles(
        qnorm_reference, master_ids
    )
    final = _read_final_profiles(activity_table, master_ids)
    if observed_reference_context != reference_context:
        raise ValueError("Configured and observed reference contexts differ")
    expected_context_keys = {
        (context, assay)
        for context in atlas_contexts
        for assay in ("atac", "h3k27ac")
    }
    if set(contexts) != expected_context_keys:
        raise ValueError("Pre-normalization context/assay set differs from configuration")
    if set(final) != set(atlas_contexts):
        raise ValueError("Final activity contexts differ from configuration")
    if set(references) != {"atac", "h3k27ac"}:
        raise ValueError("Reference table must contain ATAC and H3K27ac")

    correlations, unreplicated = _replicate_qc(libraries)
    library_summaries = []
    for library_id, profile in sorted(
        libraries.items(),
        key=lambda item: (
            item[1]["cohort"],
            item[1]["context"],
            item[1]["assay"],
            item[0],
        ),
    ):
        ordered = sorted(profile["values"])
        library_summaries.append(
            {
                "library_id": library_id,
                "cohort": profile["cohort"],
                "context": profile["context"],
                "assay": profile["assay"],
                "total_units": profile["total_units"],
                "zero_fraction": profile["zero_count"] / len(master_ids),
                **{
                    field: _quantile(ordered, probability)
                    for field, probability in zip(QUANTILE_FIELDS, QUANTILES)
                },
            }
        )

    distributions = []
    for assay in ("atac", "h3k27ac"):
        distributions.append(
            _distribution_record(
                scope="reference",
                context=reference_context,
                assay=assay,
                stage="reference",
                values=references[assay],
            )
        )
    reference_matching = []
    for context in atlas_contexts:
        for assay in ("atac", "h3k27ac"):
            pre = contexts[(context, assay)]["values"]
            post = final[context][assay]
            pre_record = _distribution_record(
                scope="atlas",
                context=context,
                assay=assay,
                stage="pre_qnorm",
                values=pre,
            )
            post_record = _distribution_record(
                scope="atlas",
                context=context,
                assay=assay,
                stage="post_qnorm",
                values=post,
            )
            distributions.extend((pre_record, post_record))
            reference_matching.append(
                {
                    "context": context,
                    "assay": assay,
                    "pre_zero_fraction": pre_record["zero_fraction"],
                    "post_zero_fraction": post_record["zero_fraction"],
                    "pre_tie_group_count": pre_record["tie_group_count"],
                    "pre_largest_tie_group": pre_record["largest_tie_group"],
                    "post_tie_group_count": post_record["tie_group_count"],
                    "post_largest_tie_group": post_record["largest_tie_group"],
                    "pre": _sorted_distance(pre, references[assay]),
                    "post": _sorted_distance(post, references[assay]),
                }
            )
    combined_activity = [
        _distribution_record(
            scope="atlas",
            context=context,
            assay="combined",
            stage="activity",
            values=final[context]["activity"],
        )
        for context in atlas_contexts
    ]
    metrics = {
        "status": "descriptive_qc_complete",
        "schema_version": 1,
        "normalization": "cpm_per_kb_then_tie_aware_reference_qnorm_v1",
        "automatic_acceptance_thresholds": False,
        "master_dhs_count": len(master_ids),
        "atlas_context_count": len(atlas_contexts),
        "library_count": len(libraries),
        "reference_context": reference_context,
        "input_sha256": {
            "library_signal": sha256_file(library_signal),
            "context_signal": sha256_file(context_signal),
            "qnorm_reference": sha256_file(qnorm_reference),
            "activity_table": sha256_file(activity_table),
            "activity_provenance": sha256_file(activity_provenance),
        },
        "library_summaries": library_summaries,
        "replicate_correlations": correlations,
        "groups_without_replicates": unreplicated,
        "reference_matching": reference_matching,
        "combined_activity_distributions": combined_activity,
    }
    correlation_content = _tsv(CORRELATION_FIELDS, correlations).encode()
    distribution_content = _tsv(DISTRIBUTION_FIELDS, distributions).encode()
    metrics_content = (json.dumps(metrics, indent=2, sort_keys=True) + "\n").encode()
    report_content = _build_html(metrics, distributions).encode()
    _atomic_write(output_correlations, correlation_content)
    _atomic_write(output_distributions, distribution_content)
    _atomic_write(output_metrics, metrics_content)
    _atomic_write(output_report, report_content)
    return metrics
