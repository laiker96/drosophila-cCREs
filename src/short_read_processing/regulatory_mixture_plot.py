"""Plot guarded H3K27ac mixture fits for context-member master DHSs."""

from __future__ import annotations

import csv
import gzip
import html
import math
from pathlib import Path
from typing import Any

from .activity import _tsv_content, sha256_file, write_deterministic_gzip, write_json_if_changed
from .activity_tmm import _atomic_text_if_changed
from .regulatory_elements import CATALOG_FIELDS, MIXTURE_FIELDS


BIN_FIELDS = [
    "context",
    "bin_index",
    "bin_left_log10",
    "bin_right_log10",
    "bin_center_log10",
    "empirical_density",
    "low_component_density",
    "high_component_density",
    "mixture_density",
]


def _quantile(ordered: list[float], probability: float) -> float:
    position = (len(ordered) - 1) * probability
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def _normal_density(value: float, mean: float, sd: float) -> float:
    return math.exp(-0.5 * ((value - mean) / sd) ** 2) / (
        sd * math.sqrt(2.0 * math.pi)
    )


def _read_mixtures(path: Path) -> tuple[list[str], dict[str, dict[str, Any]]]:
    contexts = []
    mixtures = {}
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != MIXTURE_FIELDS:
            raise ValueError(f"{path}: unexpected mixture-model columns")
        for row in reader:
            context = row["context"]
            if not context or context in mixtures:
                raise ValueError(f"{path}: duplicate or empty context")
            parsed = {
                **row,
                "positive_member_n": int(row["positive_member_n"]),
                "mixture_supported": bool(int(row["mixture_supported"])),
            }
            for field in (
                "low_mean_log10",
                "high_mean_log10",
                "low_sd_log10",
                "high_sd_log10",
                "low_weight",
                "high_weight",
                "ashman_d",
                "delta_bic",
                "posterior_crossing_log10",
            ):
                parsed[field] = float(row[field]) if row[field] else None
            fitted_fields = (
                "low_mean_log10",
                "high_mean_log10",
                "low_sd_log10",
                "high_sd_log10",
                "low_weight",
                "high_weight",
            )
            fitted_n = sum(parsed[field] is not None for field in fitted_fields)
            if fitted_n not in (0, len(fitted_fields)):
                raise ValueError(f"{path}: context {context} has partial fitted parameters")
            if parsed["mixture_supported"] and fitted_n == 0:
                raise ValueError(
                    f"{path}: supported context {context} lacks fitted parameters"
                )
            contexts.append(context)
            mixtures[context] = parsed
    if not contexts:
        raise ValueError(f"{path}: no mixture models")
    return contexts, mixtures


def _read_member_values(
    path: Path,
    *,
    contexts: list[str],
) -> dict[str, list[float]]:
    values = {context: [] for context in contexts}
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != CATALOG_FIELDS:
            raise ValueError(f"{path}: unexpected regulatory-catalog columns")
        for row in reader:
            context = row["context"]
            if context not in values:
                raise ValueError(f"{path}: unexpected context {context!r}")
            if int(row["context_membership"]) != 1:
                continue
            signal = float(row["h3k27ac_max_500_normalized_cpm_per_kb"])
            if signal > 0:
                values[context].append(math.log10(signal))
    for context in contexts:
        values[context].sort()
    return values


def _distribution_rows(
    *,
    contexts: list[str],
    mixtures: dict[str, dict[str, Any]],
    values: dict[str, list[float]],
    bin_n: int,
) -> tuple[list[dict[str, Any]], dict[str, dict[str, float]]]:
    rows = []
    limits = {}
    for context in contexts:
        observed = values[context]
        mixture = mixtures[context]
        if len(observed) != mixture["positive_member_n"]:
            raise ValueError(
                f"Context {context}: catalog positive-member count differs from mixture fit"
            )
        has_fit = mixture["low_mean_log10"] is not None
        if not observed:
            lower, upper = 0.0, 1.0
        elif has_fit:
            lower = max(
                observed[0],
                min(
                    _quantile(observed, 0.0025),
                    mixture["low_mean_log10"] - 4 * mixture["low_sd_log10"],
                ),
            )
            upper = min(
                observed[-1],
                max(
                    _quantile(observed, 0.9975),
                    mixture["high_mean_log10"] + 4 * mixture["high_sd_log10"],
                ),
            )
        else:
            lower = _quantile(observed, 0.0025)
            upper = _quantile(observed, 0.9975)
        if not lower < upper:
            center = observed[0] if observed else 0.5
            lower, upper = center - 0.5, center + 0.5
        width = (upper - lower) / bin_n
        counts = [0] * bin_n
        clipped_low = clipped_high = 0
        for value in observed:
            if value < lower:
                clipped_low += 1
                continue
            if value > upper:
                clipped_high += 1
                continue
            index = min(bin_n - 1, int((value - lower) / width))
            counts[index] += 1
        for index, count in enumerate(counts):
            left = lower + index * width
            right = left + width
            center = (left + right) / 2
            low_density = (
                mixture["low_weight"]
                * _normal_density(
                    center,
                    mixture["low_mean_log10"],
                    mixture["low_sd_log10"],
                )
                if has_fit
                else 0.0
            )
            high_density = (
                mixture["high_weight"]
                * _normal_density(
                    center,
                    mixture["high_mean_log10"],
                    mixture["high_sd_log10"],
                )
                if has_fit
                else 0.0
            )
            rows.append(
                {
                    "context": context,
                    "bin_index": index + 1,
                    "bin_left_log10": left,
                    "bin_right_log10": right,
                    "bin_center_log10": center,
                    "empirical_density": (
                        count / (len(observed) * width) if observed else 0.0
                    ),
                    "low_component_density": low_density,
                    "high_component_density": high_density,
                    "mixture_density": low_density + high_density,
                }
            )
        limits[context] = {
            "lower": lower,
            "upper": upper,
            "clipped_low_n": clipped_low,
            "clipped_high_n": clipped_high,
        }
    return rows, limits


def _polyline(
    rows: list[dict[str, Any]],
    field: str,
    *,
    x0: float,
    y0: float,
    width: float,
    height: float,
    lower: float,
    upper: float,
    maximum: float,
) -> str:
    points = []
    for row in rows:
        x = x0 + (row["bin_center_log10"] - lower) / (upper - lower) * width
        y = y0 + height * (1.0 - row[field] / maximum)
        points.append(f"{x:.2f},{y:.2f}")
    return " ".join(points)


def _build_svg(
    *,
    contexts: list[str],
    mixtures: dict[str, dict[str, Any]],
    distribution_rows: list[dict[str, Any]],
    limits: dict[str, dict[str, float]],
) -> str:
    columns = min(3, len(contexts))
    rows_n = math.ceil(len(contexts) / columns)
    panel_width, panel_height = 410, 270
    top, bottom = 76, 48
    width = columns * panel_width
    height = top + rows_n * panel_height + bottom
    by_context = {
        context: [row for row in distribution_rows if row["context"] == context]
        for context in contexts
    }
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" role="img" aria-label="H3K27ac mixture distributions">',
        "<style>text{font-family:system-ui,-apple-system,sans-serif;fill:#253047;font-size:10px}.title{font-size:13px;font-weight:700}.subtitle{font-size:9px;fill:#526176}.axis{stroke:#6b778c;stroke-width:1}.grid{stroke:#e4e9ef;stroke-width:1}.hist{fill:#aab6c5;fill-opacity:.62}.low{stroke:#2878b5}.high{stroke:#d1495b}.sum{stroke:#202938}.cut{stroke:#16845b;stroke-dasharray:5 4}</style>",
        f'<text x="{width / 2:.1f}" y="24" text-anchor="middle" class="title">H3K27ac distributions in context-member DHSs</text>',
        f'<text x="{width / 2:.1f}" y="41" text-anchor="middle" class="subtitle">Positive log10 background-TMM max-window signal; bars are empirical density</text>',
        '<rect x="20" y="52" width="13" height="9" class="hist"/><text x="38" y="61">observed</text>',
        '<line x1="112" y1="57" x2="135" y2="57" class="low" stroke-width="2"/><text x="140" y="61">low component</text>',
        '<line x1="239" y1="57" x2="262" y2="57" class="high" stroke-width="2"/><text x="267" y="61">high component</text>',
        '<line x1="377" y1="57" x2="400" y2="57" class="sum" stroke-width="2.4"/><text x="405" y="61">mixture</text>',
        '<line x1="470" y1="57" x2="493" y2="57" class="cut" stroke-width="2"/><text x="498" y="61">posterior crossing for low/high assignment</text>',
    ]
    for context_index, context in enumerate(contexts):
        column = context_index % columns
        row_number = context_index // columns
        x0 = column * panel_width + 54
        y0 = top + row_number * panel_height + 36
        plot_width, plot_height = 330, 177
        rows = by_context[context]
        limit = limits[context]
        lower, upper = limit["lower"], limit["upper"]
        maximum = max(
            max(row[field] for row in rows)
            for field in ("empirical_density", "mixture_density")
        ) * 1.06
        if maximum <= 0:
            maximum = 1.0
        mixture = mixtures[context]
        status = (
            f"supported; D={mixture['ashman_d']:.2f}"
            if mixture["mixture_supported"]
            else f"WARNING: {mixture['support_reason']}"
        )
        parts.extend(
            [
                f'<text x="{x0}" y="{y0 - 19}" class="title">{html.escape(context)}</text>',
                f'<text x="{x0 + plot_width}" y="{y0 - 19}" text-anchor="end" class="subtitle">n={mixture["positive_member_n"]:,}; {html.escape(status)}</text>',
                f'<rect x="{x0}" y="{y0}" width="{plot_width}" height="{plot_height}" fill="#fff" stroke="#c9d2dd"/>',
            ]
        )
        for fraction in (0.0, 0.5, 1.0):
            y = y0 + plot_height * (1 - fraction)
            parts.append(
                f'<line x1="{x0}" y1="{y:.2f}" x2="{x0 + plot_width}" y2="{y:.2f}" class="grid"/>'
            )
            parts.append(
                f'<text x="{x0 - 5}" y="{y + 3:.2f}" text-anchor="end">{maximum * fraction:.2g}</text>'
            )
        bin_width_pixels = plot_width / len(rows)
        for index, record in enumerate(rows):
            bar_height = record["empirical_density"] / maximum * plot_height
            parts.append(
                f'<rect x="{x0 + index * bin_width_pixels:.2f}" y="{y0 + plot_height - bar_height:.2f}" width="{bin_width_pixels + 0.15:.2f}" height="{bar_height:.2f}" class="hist"/>'
            )
        if mixture["low_mean_log10"] is not None:
            for field, css_class, stroke_width in (
                ("low_component_density", "low", 2.0),
                ("high_component_density", "high", 2.0),
                ("mixture_density", "sum", 2.4),
            ):
                points = _polyline(
                    rows,
                    field,
                    x0=x0,
                    y0=y0,
                    width=plot_width,
                    height=plot_height,
                    lower=lower,
                    upper=upper,
                    maximum=maximum,
                )
                parts.append(
                    f'<polyline points="{points}" fill="none" class="{css_class}" stroke-width="{stroke_width}"/>'
                )
        crossing = mixture["posterior_crossing_log10"]
        if crossing is not None and lower <= crossing <= upper:
            x = x0 + (crossing - lower) / (upper - lower) * plot_width
            parts.append(
                f'<line x1="{x:.2f}" y1="{y0}" x2="{x:.2f}" y2="{y0 + plot_height}" class="cut" stroke-width="2"/>'
            )
        for fraction in (0.0, 0.5, 1.0):
            value = lower + fraction * (upper - lower)
            x = x0 + fraction * plot_width
            parts.append(
                f'<text x="{x:.2f}" y="{y0 + plot_height + 15}" text-anchor="middle">{value:.2f}</text>'
            )
        parts.append(
            f'<text x="{x0 + plot_width / 2}" y="{y0 + plot_height + 31}" text-anchor="middle">log10(max-window normalized H3K27ac)</text>'
        )
        parts.append(
            f'<text x="{x0 - 39}" y="{y0 + plot_height / 2}" text-anchor="middle" transform="rotate(-90 {x0 - 39} {y0 + plot_height / 2})">density</text>'
        )
    parts.append(
        f'<text x="{width / 2:.1f}" y="{height - 13}" text-anchor="middle" class="subtitle">Display range uses empirical q0.25%–q99.75% and fitted means ±4 SD when available; warnings name failed BIC, separation, weight, crossing, sample-size, or variance guards.</text>'
    )
    parts.append("</svg>")
    return "".join(parts)


def build_mixture_distribution_plot(
    *,
    catalog_path: Path,
    mixture_path: Path,
    output_svg: Path,
    output_bins: Path,
    output_metrics: Path,
    bin_n: int = 70,
) -> dict[str, Any]:
    """Create a deterministic faceted SVG and machine-readable density table."""

    if bin_n < 20:
        raise ValueError("At least 20 histogram bins are required")
    contexts, mixtures = _read_mixtures(mixture_path)
    values = _read_member_values(catalog_path, contexts=contexts)
    rows, limits = _distribution_rows(
        contexts=contexts,
        mixtures=mixtures,
        values=values,
        bin_n=bin_n,
    )
    write_deterministic_gzip(output_bins, _tsv_content(BIN_FIELDS, rows))
    _atomic_text_if_changed(
        output_svg,
        _build_svg(
            contexts=contexts,
            mixtures=mixtures,
            distribution_rows=rows,
            limits=limits,
        ),
    )
    metrics = {
        "method": "guarded_h3k27ac_mixture_distribution_plot_v1",
        "population": "positive_context_member_master_dhs",
        "signal": "log10_h3k27ac_max_500_normalized_cpm_per_kb",
        "histogram_bin_n": bin_n,
        "contexts": contexts,
        "supported_contexts": [
            context for context in contexts if mixtures[context]["mixture_supported"]
        ],
        "unsupported_contexts": [
            context for context in contexts if not mixtures[context]["mixture_supported"]
        ],
        "display_limits": limits,
        "inputs": {
            "catalog": {"path": str(catalog_path.resolve()), "sha256": sha256_file(catalog_path)},
            "mixtures": {"path": str(mixture_path.resolve()), "sha256": sha256_file(mixture_path)},
        },
        "outputs": {
            "svg": {"path": str(output_svg.resolve()), "sha256": sha256_file(output_svg)},
            "bins": {"path": str(output_bins.resolve()), "sha256": sha256_file(output_bins)},
        },
    }
    write_json_if_changed(output_metrics, metrics)
    return metrics
