"""Background-TMM normalization for master-DHS activity values."""

from __future__ import annotations

from contextlib import ExitStack
import csv
import gzip
import io
import math
import os
from pathlib import Path
import statistics
import tempfile
from typing import Any, Iterable

from .activity import (
    LIBRARY_SIGNAL_FIELDS,
    _aggregate_rows,
    _format_number,
    _read_library_signal,
    _tsv_content,
    sha256_file,
    write_deterministic_gzip,
    write_json_if_changed,
)


TMM_MASTER_METHOD = "tmm_master_dhs_v1"
TMM_BACKGROUND_METHOD = "tmm_background_10kb_v1"
TMM_METHODS = (TMM_MASTER_METHOD, TMM_BACKGROUND_METHOD)
BACKGROUND_BIN_WIDTH = 10_000

BACKGROUND_COUNT_FIELDS = [
    "background_bin_id",
    "chrom",
    "start",
    "end",
    "raw_count",
]
TMM_METADATA_FIELDS = [
    "library_id",
    "assay",
    "context",
    "total_units",
    "signal_sha256",
    "count_source",
    "count_source_sha256",
]
TMM_FACTOR_FIELDS = [
    "library_id",
    "assay",
    "context",
    "total_units",
    "feature_count",
    "tmm_normalization_factor",
    "effective_library_size",
    "normalization_method",
]
TMM_CONTEXT_SIGNAL_FIELDS = [
    "master_dhs_id",
    "chrom",
    "start",
    "end",
    "summit",
    "width_bp",
    "context",
    "assay",
    "normalization_method",
    "library_n",
    "library_ids",
    "raw_count_sum",
    "raw_count_mean",
    "total_units_sum",
    "cpm_per_kb",
    "cpm_per_kb_sd",
    "normalized_cpm_per_kb",
    "normalized_cpm_per_kb_sd",
]
TMM_ACTIVITY_FIELDS = [
    "master_dhs_id",
    "chrom",
    "start",
    "end",
    "summit",
    "width_bp",
    "context",
    "normalization_method",
    "atac_library_n",
    "atac_library_ids",
    "atac_raw_count_sum",
    "atac_raw_count_mean",
    "atac_total_units_sum",
    "atac_cpm_per_kb",
    "atac_cpm_per_kb_sd",
    "atac_normalized_cpm_per_kb",
    "atac_normalized_cpm_per_kb_sd",
    "h3k27ac_library_n",
    "h3k27ac_library_ids",
    "h3k27ac_raw_count_sum",
    "h3k27ac_raw_count_mean",
    "h3k27ac_total_units_sum",
    "h3k27ac_cpm_per_kb",
    "h3k27ac_cpm_per_kb_sd",
    "h3k27ac_normalized_cpm_per_kb",
    "h3k27ac_normalized_cpm_per_kb_sd",
    "activity",
]


def _atomic_text_if_changed(path: Path, content: str) -> None:
    encoded = content.encode("utf-8")
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.is_file() and path.read_bytes() == encoded:
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
            handle.write(encoded)
        os.replace(temporary_name, path)
    finally:
        if temporary_name and os.path.exists(temporary_name):
            os.unlink(temporary_name)


def write_background_bins(
    *,
    chrom_sizes_path: Path,
    autosomes_path: Path,
    output_path: Path,
    bin_width: int = BACKGROUND_BIN_WIDTH,
) -> int:
    """Write fixed genomic bins on configured autosomes in reference order."""

    if bin_width < 1:
        raise ValueError("Background bin width must be positive")
    autosomes = [
        line.strip()
        for line in autosomes_path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if not autosomes or len(autosomes) != len(set(autosomes)):
        raise ValueError("Autosome list must be non-empty and unique")
    wanted = set(autosomes)
    sizes: dict[str, int] = {}
    order = []
    with chrom_sizes_path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 2 or not fields[0]:
                raise ValueError(
                    f"{chrom_sizes_path}:{line_number}: expected chromosome and size"
                )
            chrom, size_text = fields
            size = int(size_text)
            if size < 1 or chrom in sizes:
                raise ValueError(
                    f"{chrom_sizes_path}:{line_number}: invalid chromosome size"
                )
            sizes[chrom] = size
            order.append(chrom)
    missing = sorted(wanted - set(sizes))
    if missing:
        raise ValueError("Autosomes absent from chromosome sizes: " + ", ".join(missing))

    lines = []
    bin_number = 0
    for chrom in order:
        if chrom not in wanted:
            continue
        for start in range(0, sizes[chrom], bin_width):
            end = min(start + bin_width, sizes[chrom])
            bin_number += 1
            lines.append(
                f"{chrom}\t{start}\t{end}\tBGBIN{bin_number:07d}\n"
            )
    _atomic_text_if_changed(output_path, "".join(lines))
    return bin_number


def _write_gzip_dict_rows(
    path: Path,
    fields: list[str],
    rows: Iterable[dict[str, Any]],
) -> None:
    """Atomically write deterministic gzip TSV without buffering every row."""

    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_name = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb",
            prefix=f".{path.name}.",
            dir=path.parent,
            delete=False,
        ) as raw:
            temporary_name = raw.name
            with gzip.GzipFile(
                filename="",
                fileobj=raw,
                mode="wb",
                compresslevel=6,
                mtime=0,
            ) as compressed:
                with io.TextIOWrapper(
                    compressed,
                    encoding="utf-8",
                    newline="",
                ) as text:
                    writer = csv.DictWriter(
                        text,
                        fieldnames=fields,
                        delimiter="\t",
                        lineterminator="\n",
                    )
                    writer.writeheader()
                    for raw_row in rows:
                        writer.writerow(
                            {
                                field: (
                                    _format_number(raw_row[field])
                                    if isinstance(raw_row.get(field), (int, float))
                                    else raw_row.get(field, "")
                                )
                                for field in fields
                            }
                        )
        if path.is_file() and path.read_bytes() == Path(temporary_name).read_bytes():
            os.unlink(temporary_name)
            temporary_name = None
            return
        os.replace(temporary_name, path)
        temporary_name = None
    finally:
        if temporary_name and os.path.exists(temporary_name):
            os.unlink(temporary_name)


def _signal_metadata(signal_paths: dict[str, Path]) -> dict[str, dict[str, Any]]:
    metadata = {}
    for library_id, path in signal_paths.items():
        with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != LIBRARY_SIGNAL_FIELDS:
                raise ValueError(f"{path}: unexpected library-signal columns")
            row = next(reader, None)
        if row is None:
            raise ValueError(f"{path}: library signal is empty")
        if row["library_id"] != library_id:
            raise ValueError(f"{path}: library ID does not match {library_id!r}")
        if row["cohort"] != "atlas" or row["assay"] not in {"atac", "h3k27ac"}:
            raise ValueError(f"{path}: TMM inputs must be atlas ATAC/H3K27ac libraries")
        total_units = int(row["total_units"])
        if total_units < 1:
            raise ValueError(f"{path}: total_units must be positive")
        metadata[library_id] = {
            "library_id": library_id,
            "assay": row["assay"],
            "context": row["context"],
            "total_units": total_units,
            "signal_sha256": sha256_file(path),
        }
    return metadata


def build_tmm_inputs(
    *,
    method: str,
    signal_paths: dict[str, Path],
    background_count_paths: dict[str, Path] | None,
    output_counts: Path,
    output_metadata: Path,
) -> dict[str, Any]:
    """Build the raw count matrix and metadata consumed by edgeR TMM."""

    if method not in TMM_METHODS:
        raise ValueError(f"Unsupported TMM method: {method}")
    if not signal_paths:
        raise ValueError("TMM requires atlas libraries")
    library_ids = sorted(signal_paths)
    metadata = _signal_metadata(signal_paths)
    for assay in ("atac", "h3k27ac"):
        if sum(metadata[item]["assay"] == assay for item in library_ids) < 2:
            raise ValueError(f"TMM requires at least two atlas {assay} libraries")

    if method == TMM_MASTER_METHOD:
        count_paths = signal_paths
        input_fields = LIBRARY_SIGNAL_FIELDS
        feature_field = "master_dhs_id"
    else:
        if background_count_paths is None or set(background_count_paths) != set(
            signal_paths
        ):
            raise ValueError("Background TMM requires one count file per atlas library")
        count_paths = background_count_paths
        input_fields = BACKGROUND_COUNT_FIELDS
        feature_field = "background_bin_id"

    matrix_fields = ["feature_id", "chrom", "start", "end", "width_bp", *library_ids]
    feature_count = 0

    def matrix_rows():
        nonlocal feature_count
        with ExitStack() as stack:
            readers = []
            for library_id in library_ids:
                handle = stack.enter_context(
                    gzip.open(
                        count_paths[library_id],
                        "rt",
                        encoding="utf-8",
                        newline="",
                    )
                )
                reader = csv.DictReader(handle, delimiter="\t")
                if reader.fieldnames != input_fields:
                    raise ValueError(
                        f"{count_paths[library_id]}: unexpected count columns"
                    )
                readers.append(reader)
            while True:
                rows = [next(reader, None) for reader in readers]
                if all(row is None for row in rows):
                    break
                if any(row is None for row in rows):
                    raise ValueError("TMM count inputs have different row counts")
                assert all(row is not None for row in rows)
                first = rows[0]
                identity = (
                    first[feature_field],
                    first["chrom"],
                    int(first["start"]),
                    int(first["end"]),
                )
                if identity[3] <= identity[2]:
                    raise ValueError("TMM input contains an invalid interval")
                counts = {}
                for library_id, row in zip(library_ids, rows):
                    observed = (
                        row[feature_field],
                        row["chrom"],
                        int(row["start"]),
                        int(row["end"]),
                    )
                    if observed != identity:
                        raise ValueError("TMM count inputs have different feature order")
                    if method == TMM_MASTER_METHOD and row["library_id"] != library_id:
                        raise ValueError("Master-DHS signal library order is inconsistent")
                    raw_count = int(row["raw_count"])
                    if raw_count < 0:
                        raise ValueError("TMM counts must be non-negative")
                    counts[library_id] = raw_count
                feature_count += 1
                yield {
                    "feature_id": identity[0],
                    "chrom": identity[1],
                    "start": identity[2],
                    "end": identity[3],
                    "width_bp": identity[3] - identity[2],
                    **counts,
                }

    _write_gzip_dict_rows(output_counts, matrix_fields, matrix_rows())
    if feature_count < 1:
        raise ValueError("TMM count matrix is empty")

    metadata_rows = []
    for library_id in library_ids:
        count_source = count_paths[library_id]
        metadata_rows.append(
            {
                **metadata[library_id],
                "count_source": str(count_source.resolve()),
                "count_source_sha256": sha256_file(count_source),
            }
        )
    _atomic_text_if_changed(
        output_metadata,
        _tsv_content(TMM_METADATA_FIELDS, metadata_rows),
    )
    return {
        "method": method,
        "feature_count": feature_count,
        "library_count": len(library_ids),
    }


def _read_tmm_factors(
    path: Path,
    *,
    expected_method: str,
) -> dict[str, dict[str, Any]]:
    factors = {}
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != TMM_FACTOR_FIELDS:
            raise ValueError(f"{path}: unexpected TMM-factor columns")
        for row in reader:
            library_id = row["library_id"]
            if library_id in factors:
                raise ValueError(f"{path}: duplicate factor for {library_id}")
            factor = float(row["tmm_normalization_factor"])
            effective_size = float(row["effective_library_size"])
            total_units = int(row["total_units"])
            if (
                row["normalization_method"] != expected_method
                or row["assay"] not in {"atac", "h3k27ac"}
                or factor <= 0
                or effective_size <= 0
                or not math.isfinite(factor)
                or not math.isfinite(effective_size)
                or not math.isclose(
                    effective_size,
                    total_units * factor,
                    rel_tol=1e-10,
                )
            ):
                raise ValueError(f"{path}: invalid factor for {library_id}")
            factors[library_id] = {
                **row,
                "total_units": total_units,
                "feature_count": int(row["feature_count"]),
                "tmm_normalization_factor": factor,
                "effective_library_size": effective_size,
            }
    if not factors:
        raise ValueError(f"{path}: TMM factor table is empty")
    return factors


def build_tmm_activity_outputs(
    *,
    method: str,
    signal_paths: dict[str, Path],
    factor_path: Path,
    contexts: list[str],
    output_context_signal: Path,
    output_activity: Path,
    output_metrics: Path,
    output_provenance: Path,
    provenance: dict[str, Any],
) -> dict[str, Any]:
    """Apply assay-specific TMM factors and build a comparative activity table."""

    if method not in TMM_METHODS:
        raise ValueError(f"Unsupported TMM method: {method}")
    factors = _read_tmm_factors(factor_path, expected_method=method)
    if set(factors) != set(signal_paths):
        raise ValueError("TMM factors and atlas signal libraries differ")
    rows_by_library = {
        library_id: _read_library_signal(path)
        for library_id, path in sorted(signal_paths.items())
    }
    first_rows = next(iter(rows_by_library.values()))
    master_ids = [row["master_dhs_id"] for row in first_rows]
    metadata = {}
    normalized = {}
    for library_id, rows in rows_by_library.items():
        if [row["master_dhs_id"] for row in rows] != master_ids:
            raise ValueError(f"Master IDs differ for library {library_id}")
        first = rows[0]
        factor = factors[library_id]
        if (
            first["cohort"] != "atlas"
            or first["assay"] != factor["assay"]
            or first["context"] != factor["context"]
            or first["total_units"] != factor["total_units"]
        ):
            raise ValueError(f"TMM metadata differ for library {library_id}")
        metadata[library_id] = {
            "assay": first["assay"],
            "context": first["context"],
        }
        effective_size = factor["effective_library_size"]
        normalized[library_id] = [
            row["raw_count"] * 1_000_000_000.0
            / (effective_size * row["width_bp"])
            for row in rows
        ]

    aggregates: dict[tuple[str, str], list[dict[str, Any]]] = {}
    context_rows = []
    for context in contexts:
        for assay in ("atac", "h3k27ac"):
            library_ids = sorted(
                library_id
                for library_id, values in metadata.items()
                if values == {"assay": assay, "context": context}
            )
            if not library_ids:
                raise ValueError(f"Atlas context {context!r} lacks {assay}")
            aggregate = _aggregate_rows(rows_by_library, library_ids)
            for element_index, row in enumerate(aggregate):
                values = [
                    normalized[library_id][element_index]
                    for library_id in library_ids
                ]
                row.update(
                    {
                        "context": context,
                        "assay": assay,
                        "normalization_method": method,
                        "normalized_cpm_per_kb": statistics.fmean(values),
                        "normalized_cpm_per_kb_sd": (
                            statistics.stdev(values) if len(values) > 1 else 0.0
                        ),
                    }
                )
            aggregates[(context, assay)] = aggregate
            context_rows.extend(aggregate)
    write_deterministic_gzip(
        output_context_signal,
        _tsv_content(TMM_CONTEXT_SIGNAL_FIELDS, context_rows),
    )

    final_rows = []
    for element_index, base in enumerate(first_rows):
        for context in contexts:
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
                "normalization_method": method,
            }
            for assay, prefix in (("atac", "atac"), ("h3k27ac", "h3k27ac")):
                values = aggregates[(context, assay)][element_index]
                for source, suffix in (
                    ("library_n", "library_n"),
                    ("library_ids", "library_ids"),
                    ("raw_count_sum", "raw_count_sum"),
                    ("raw_count_mean", "raw_count_mean"),
                    ("total_units_sum", "total_units_sum"),
                    ("cpm_per_kb", "cpm_per_kb"),
                    ("cpm_per_kb_sd", "cpm_per_kb_sd"),
                    ("normalized_cpm_per_kb", "normalized_cpm_per_kb"),
                    ("normalized_cpm_per_kb_sd", "normalized_cpm_per_kb_sd"),
                ):
                    row[f"{prefix}_{suffix}"] = values[source]
            row["activity"] = math.sqrt(
                row["atac_normalized_cpm_per_kb"]
                * row["h3k27ac_normalized_cpm_per_kb"]
            )
            final_rows.append(row)
    write_deterministic_gzip(
        output_activity,
        _tsv_content(TMM_ACTIVITY_FIELDS, final_rows),
    )

    assay_factor_products = {}
    for assay in ("atac", "h3k27ac"):
        assay_factors = [
            value["tmm_normalization_factor"]
            for value in factors.values()
            if value["assay"] == assay
        ]
        assay_factor_products[assay] = math.prod(assay_factors)
    metrics = {
        "status": "ok",
        "schema_version": 1,
        "normalization_method": method,
        "activity_formula": "sqrt_atac_times_h3k27ac_v1",
        "master_dhs_count": len(master_ids),
        "context_count": len(contexts),
        "library_count": len(rows_by_library),
        "factor_product_by_assay": assay_factor_products,
        "factor_range_by_assay": {
            assay: {
                "minimum": min(
                    value["tmm_normalization_factor"]
                    for value in factors.values()
                    if value["assay"] == assay
                ),
                "maximum": max(
                    value["tmm_normalization_factor"]
                    for value in factors.values()
                    if value["assay"] == assay
                ),
            }
            for assay in ("atac", "h3k27ac")
        },
    }
    write_json_if_changed(output_metrics, metrics)
    output_paths = {
        "context_signal": output_context_signal,
        "activity_table": output_activity,
        "metrics": output_metrics,
    }
    write_json_if_changed(
        output_provenance,
        {
            **provenance,
            "schema_version": 1,
            "normalization": method,
            "activity_formula": "sqrt_atac_times_h3k27ac_v1",
            "normalization_factor_table": {
                "path": str(factor_path.resolve()),
                "sha256": sha256_file(factor_path),
            },
            "input_signal_sha256": {
                library_id: sha256_file(path)
                for library_id, path in sorted(signal_paths.items())
            },
            "outputs": {
                name: {"path": str(path.resolve()), "sha256": sha256_file(path)}
                for name, path in sorted(output_paths.items())
            },
        },
    )
    return metrics
