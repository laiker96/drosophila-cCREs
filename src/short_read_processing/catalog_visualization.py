"""Build BED and normalized mean-signal tracks for catalog visualization."""

from __future__ import annotations

import csv
import gzip
import math
import os
from pathlib import Path
import subprocess
import tempfile
from typing import Any, Iterable, Iterator

from .activity import read_master_elements, sha256_file, write_json_if_changed
from .activity_tmm import (
    TMM_BACKGROUND_METHOD,
    TMM_FACTOR_FIELDS,
    _atomic_text_if_changed,
)


REGULATORY_CLASS_COLORS = {
    "promoter_associated": "215,48,39",
    "proximal_enhancer_like": "253,174,97",
    "distal_enhancer_like": "69,117,180",
    "unclassified_no_tss_on_contig": "117,112,179",
}


def _read_chrom_sizes(path: Path) -> list[tuple[str, int]]:
    rows: list[tuple[str, int]] = []
    seen = set()
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 2:
                raise ValueError(f"{path}:{line_number}: expected chromosome and size")
            chrom, size_text = fields
            size = int(size_text)
            if not chrom or chrom in seen or size < 1:
                raise ValueError(f"{path}:{line_number}: invalid chromosome size")
            rows.append((chrom, size))
            seen.add(chrom)
    if not rows:
        raise ValueError(f"Chromosome-size file is empty: {path}")
    return rows


def _read_context_matrix(
    path: Path,
    contexts: list[str],
    elements,
) -> dict[tuple[str, str], bool]:
    membership: dict[tuple[str, str], bool] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "master_dhs_id",
            "chrom",
            "start",
            "end",
            "summit",
            "context_n",
            *contexts,
        }
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path}: context matrix lacks required columns")
        rows = list(reader)
    if len(rows) != len(elements):
        raise ValueError("Context matrix and master registry have different row counts")
    for line_number, (row, element) in enumerate(zip(rows, elements), start=2):
        if (
            row["master_dhs_id"] != element.master_id
            or row["chrom"] != element.chrom
            or int(row["start"]) != element.start
            or int(row["end"]) != element.end
            or int(row["summit"]) != element.summit
        ):
            raise ValueError(f"{path}:{line_number}: master coordinates differ")
        values = []
        for context in contexts:
            if row[context] not in {"0", "1"}:
                raise ValueError(f"{path}:{line_number}: invalid {context} membership")
            value = row[context] == "1"
            values.append(value)
            membership[(element.master_id, context)] = value
        matrix_context_n = sum(
            int(row[field])
            for field in reader.fieldnames
            if field not in {"master_dhs_id", "chrom", "start", "end", "summit", "context_n"}
        )
        if int(row["context_n"]) != matrix_context_n:
            raise ValueError(f"{path}:{line_number}: context_n is inconsistent")
    return membership


def _active_bed_content(
    *,
    path: Path,
    context: str,
    element_by_id: dict[str, Any],
    membership: dict[tuple[str, str], bool],
) -> tuple[str, dict[str, int]]:
    lines = []
    class_counts: dict[str, int] = {}
    seen = set()
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "master_dhs_id",
            "chrom",
            "start",
            "end",
            "summit",
            "context",
            "context_membership",
            "regulatory_class",
            "mixture_component",
            "mixture_guardrail_warning",
            "mixture_high_posterior_probability",
        }
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path}: active-element table lacks required columns")
        for line_number, row in enumerate(reader, start=2):
            master_id = row["master_dhs_id"]
            element = element_by_id.get(master_id)
            if element is None or master_id in seen:
                raise ValueError(f"{path}:{line_number}: unknown or duplicate master DHS")
            if (
                row["context"] != context
                or row["context_membership"] != "1"
                or row["mixture_component"] != "high"
                or not membership[(master_id, context)]
                or row["chrom"] != element.chrom
                or int(row["start"]) != element.start
                or int(row["end"]) != element.end
                or int(row["summit"]) != element.summit
            ):
                raise ValueError(f"{path}:{line_number}: inconsistent active element")
            element_class = row["regulatory_class"]
            color = REGULATORY_CLASS_COLORS.get(element_class, "105,105,105")
            posterior = float(row["mixture_high_posterior_probability"])
            if not math.isfinite(posterior) or not 0 <= posterior <= 1:
                raise ValueError(f"{path}:{line_number}: invalid mixture posterior")
            warning = row["mixture_guardrail_warning"] == "1"
            state = "high_mixture_warning" if warning else "high_mixture_supported"
            name = f"{master_id}|{element_class}|{state}"
            score = round(1000 * posterior)
            lines.append(
                f"{element.chrom}\t{element.start}\t{element.end}\t{name}\t"
                f"{score}\t.\t{element.summit}\t{element.summit + 1}\t{color}\n"
            )
            class_counts[element_class] = class_counts.get(element_class, 0) + 1
            seen.add(master_id)
    return "".join(lines), class_counts


def build_catalog_beds(
    *,
    master_bed: Path,
    summit_bed: Path,
    context_matrix: Path,
    active_paths: dict[str, Path],
    output_master_bed: Path,
    output_context_dhs: dict[str, Path],
    output_active_beds: dict[str, Path],
    output_manifest: Path,
) -> dict[str, Any]:
    """Write portable master, context-DHS, and active-element BED tracks."""

    contexts = list(active_paths)
    if (
        not contexts
        or set(output_context_dhs) != set(contexts)
        or set(output_active_beds) != set(contexts)
    ):
        raise ValueError("Catalog BED context mappings must be complete and identical")
    elements = read_master_elements(master_bed, summit_bed)
    element_by_id = {element.master_id: element for element in elements}
    membership = _read_context_matrix(context_matrix, contexts, elements)

    _atomic_text_if_changed(
        output_master_bed,
        master_bed.read_text(encoding="utf-8"),
    )
    context_metrics = {}
    for context in contexts:
        dhs_lines = []
        for element in elements:
            if membership[(element.master_id, context)]:
                dhs_lines.append(
                    f"{element.chrom}\t{element.start}\t{element.end}\t"
                    f"{element.master_id}\t0\t.\t{element.summit}\t"
                    f"{element.summit + 1}\t0,145,130\n"
                )
        _atomic_text_if_changed(output_context_dhs[context], "".join(dhs_lines))
        active_content, class_counts = _active_bed_content(
            path=active_paths[context],
            context=context,
            element_by_id=element_by_id,
            membership=membership,
        )
        _atomic_text_if_changed(output_active_beds[context], active_content)
        active_n = sum(class_counts.values())
        if active_n > len(dhs_lines):
            raise ValueError(f"Active-element count exceeds context DHSs for {context}")
        context_metrics[context] = {
            "context_dhs_count": len(dhs_lines),
            "active_element_count": active_n,
            "active_by_regulatory_class": dict(sorted(class_counts.items())),
        }

    outputs = {
        "master_dhs_bed": {
            "path": str(output_master_bed.resolve()),
            "sha256": sha256_file(output_master_bed),
        },
        "contexts": {
            context: {
                "context_dhs_bed": {
                    "path": str(output_context_dhs[context].resolve()),
                    "sha256": sha256_file(output_context_dhs[context]),
                },
                "active_elements_bed": {
                    "path": str(output_active_beds[context].resolve()),
                    "sha256": sha256_file(output_active_beds[context]),
                },
            }
            for context in contexts
        },
    }
    manifest = {
        "status": "ok",
        "schema_version": 1,
        "method": "catalog_igv_bed_tracks_v1",
        "context_dhs_definition": "master DHS coordinates with context_matrix membership=1",
        "active_element_definition": "context-member master DHS assigned to the high H3K27ac mixture",
        "bed_format": "BED9; thick interval marks the representative summit",
        "regulatory_class_item_rgb": REGULATORY_CLASS_COLORS,
        "master_dhs_count": len(elements),
        "context_metrics": context_metrics,
        "inputs": {
            "master_bed": sha256_file(master_bed),
            "summit_bed": sha256_file(summit_bed),
            "context_matrix": sha256_file(context_matrix),
            "active_tables": {
                context: sha256_file(active_paths[context]) for context in contexts
            },
        },
        "outputs": outputs,
    }
    write_json_if_changed(output_manifest, manifest)
    return manifest


def read_track_factors(
    *,
    factor_path: Path,
    library_ids: list[str],
    assay: str,
    context: str,
) -> dict[str, dict[str, float | int | str]]:
    """Return validated TMM effective sizes for one context and assay."""

    selected: dict[str, dict[str, float | int | str]] = {}
    wanted = set(library_ids)
    if not library_ids or len(wanted) != len(library_ids):
        raise ValueError("Mean BigWig libraries must be non-empty and unique")
    with factor_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != TMM_FACTOR_FIELDS:
            raise ValueError(f"{factor_path}: unexpected TMM-factor columns")
        for row in reader:
            library_id = row["library_id"]
            if library_id not in wanted:
                continue
            effective_size = float(row["effective_library_size"])
            total_units = int(row["total_units"])
            normalization_factor = float(row["tmm_normalization_factor"])
            if (
                library_id in selected
                or row["assay"] != assay
                or row["context"] != context
                or row["normalization_method"] != TMM_BACKGROUND_METHOD
                or not math.isfinite(effective_size)
                or effective_size <= 0
                or not math.isfinite(normalization_factor)
                or normalization_factor <= 0
                or total_units < 1
            ):
                raise ValueError(f"Invalid TMM factor for {library_id}")
            selected[library_id] = {
                "total_units": total_units,
                "tmm_normalization_factor": normalization_factor,
                "effective_library_size": effective_size,
                "scale": 1_000_000.0 / effective_size,
            }
    if set(selected) != wanted:
        missing = sorted(wanted - set(selected))
        raise ValueError("TMM factors are missing libraries: " + ", ".join(missing))
    return {library_id: selected[library_id] for library_id in library_ids}


def mean_unionbedg_rows(
    lines: Iterable[str],
    *,
    chromosome_sizes: list[tuple[str, int]],
    library_n: int,
) -> Iterator[tuple[str, int, int, float]]:
    """Validate unionbedg rows and yield nonzero arithmetic mean segments."""

    if library_n < 1:
        raise ValueError("library_n must be positive")
    chrom_order = {chrom: index for index, (chrom, _) in enumerate(chromosome_sizes)}
    sizes = dict(chromosome_sizes)
    previous: tuple[int, int, int] | None = None
    for line_number, line in enumerate(lines, start=1):
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 3 + library_n:
            raise ValueError(f"unionbedg line {line_number}: unexpected column count")
        chrom = fields[0]
        if chrom not in chrom_order:
            raise ValueError(f"unionbedg line {line_number}: unknown chromosome {chrom}")
        start, end = int(fields[1]), int(fields[2])
        values = [float(value) for value in fields[3:]]
        if (
            start < 0
            or end <= start
            or end > sizes[chrom]
            or any(not math.isfinite(value) or value < 0 for value in values)
        ):
            raise ValueError(f"unionbedg line {line_number}: invalid interval or value")
        current = (chrom_order[chrom], start, end)
        if previous is not None and (
            current[0] < previous[0]
            or (current[0] == previous[0] and start < previous[2])
        ):
            raise ValueError("unionbedg rows are unsorted or overlapping")
        previous = current
        mean_value = sum(values) / library_n
        if mean_value > 0:
            yield chrom, start, end, mean_value


def _coverage_bedgraph(
    *,
    unit_path: Path,
    chrom_sizes_path: Path,
    scale: float,
    output_path: Path,
    threads: int,
) -> None:
    with output_path.open("wb") as output_handle:
        decompressor = subprocess.Popen(
            ["pigz", "-p", str(max(1, threads)), "-dc", str(unit_path)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        assert decompressor.stdout is not None
        coverage = subprocess.Popen(
            [
                "bedtools",
                "genomecov",
                "-bg",
                "-scale",
                f"{scale:.17g}",
                "-i",
                "stdin",
                "-g",
                str(chrom_sizes_path),
            ],
            stdin=decompressor.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        decompressor.stdout.close()
        assert coverage.stdout is not None
        sorter = subprocess.Popen(
            [
                "sort",
                "-k1,1",
                "-k2,2n",
                "--temporary-directory",
                str(output_path.parent),
            ],
            stdin=coverage.stdout,
            stdout=output_handle,
            stderr=subprocess.PIPE,
            env={**os.environ, "LC_ALL": "C"},
        )
        coverage.stdout.close()
        _, coverage_stderr = coverage.communicate()
        _, decompressor_stderr = decompressor.communicate()
        _, sorter_stderr = sorter.communicate()
    if decompressor.returncode or coverage.returncode or sorter.returncode:
        raise RuntimeError(
            f"Could not make coverage for {unit_path}: "
            + decompressor_stderr.decode(errors="replace")
            + coverage_stderr.decode(errors="replace")
            + sorter_stderr.decode(errors="replace")
        )


def _replace_if_changed(temporary: Path, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    if (
        output.is_file()
        and output.stat().st_size == temporary.stat().st_size
        and sha256_file(output) == sha256_file(temporary)
    ):
        temporary.unlink()
        return
    os.replace(temporary, output)


def build_context_mean_bigwig(
    *,
    unit_paths: dict[str, Path],
    factor_path: Path,
    chrom_sizes_path: Path,
    assay: str,
    context: str,
    output_bigwig: Path,
    output_metrics: Path,
    threads: int = 2,
) -> dict[str, Any]:
    """Write the mean of per-library background-TMM normalized coverage tracks."""

    if assay not in {"atac", "h3k27ac"}:
        raise ValueError(f"Unsupported track assay: {assay}")
    library_ids = list(unit_paths)
    factors = read_track_factors(
        factor_path=factor_path,
        library_ids=library_ids,
        assay=assay,
        context=context,
    )
    chromosome_sizes = sorted(_read_chrom_sizes(chrom_sizes_path))
    output_bigwig.parent.mkdir(parents=True, exist_ok=True)
    temporary_bigwig: Path | None = None
    try:
        with tempfile.TemporaryDirectory(
            prefix=f".{context}.{assay}.", dir=output_bigwig.parent
        ) as temporary_dir_text:
            temporary_dir = Path(temporary_dir_text)
            bedgraphs = []
            for index, library_id in enumerate(library_ids):
                bedgraph = temporary_dir / f"{index:04d}.bedgraph"
                _coverage_bedgraph(
                    unit_path=unit_paths[library_id],
                    chrom_sizes_path=chrom_sizes_path,
                    scale=float(factors[library_id]["scale"]),
                    output_path=bedgraph,
                    threads=threads,
                )
                bedgraphs.append(bedgraph)

            try:
                import pyBigWig
            except ImportError as error:  # pragma: no cover - rule environment owns this
                raise RuntimeError("pyBigWig is required to write catalog tracks") from error

            union = None
            if len(bedgraphs) == 1:
                mean_stream = bedgraphs[0].open(encoding="utf-8")
            else:
                union = subprocess.Popen(
                    [
                        "bedtools",
                        "unionbedg",
                        "-filler",
                        "0",
                        "-i",
                        *map(str, bedgraphs),
                    ],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                )
                assert union.stdout is not None
                mean_stream = union.stdout

            with tempfile.NamedTemporaryFile(
                prefix=f".{output_bigwig.name}.",
                dir=output_bigwig.parent,
                delete=False,
            ) as handle:
                temporary_bigwig = Path(handle.name)
            bigwig = pyBigWig.open(str(temporary_bigwig), "w")
            if bigwig is None:
                raise RuntimeError(f"Could not create BigWig: {temporary_bigwig}")
            try:
                bigwig.addHeader(chromosome_sizes)
                chrom_buffer: list[str] = []
                start_buffer: list[int] = []
                end_buffer: list[int] = []
                value_buffer: list[float] = []
                for chrom, start, end, value in mean_unionbedg_rows(
                    mean_stream,
                    chromosome_sizes=chromosome_sizes,
                    library_n=len(library_ids),
                ):
                    chrom_buffer.append(chrom)
                    start_buffer.append(start)
                    end_buffer.append(end)
                    value_buffer.append(value)
                    if len(chrom_buffer) == 100_000:
                        bigwig.addEntries(
                            chrom_buffer,
                            start_buffer,
                            ends=end_buffer,
                            values=value_buffer,
                        )
                        chrom_buffer.clear()
                        start_buffer.clear()
                        end_buffer.clear()
                        value_buffer.clear()
                if chrom_buffer:
                    bigwig.addEntries(
                        chrom_buffer,
                        start_buffer,
                        ends=end_buffer,
                        values=value_buffer,
                    )
            finally:
                bigwig.close()
                mean_stream.close()
            if union is not None:
                union_stderr = union.stderr.read() if union.stderr is not None else ""
                if union.wait() != 0:
                    raise RuntimeError("bedtools unionbedg failed: " + union_stderr)
            _replace_if_changed(temporary_bigwig, output_bigwig)
            temporary_bigwig = None
    finally:
        if temporary_bigwig is not None and temporary_bigwig.exists():
            temporary_bigwig.unlink()

    unit_semantics = "Tn5 insertion records" if assay == "atac" else "ChIP fragments"
    metrics = {
        "status": "ok",
        "schema_version": 1,
        "method": "mean_background_tmm_normalized_coverage_v1",
        "assay": assay,
        "context": context,
        "unit_semantics": unit_semantics,
        "value_definition": (
            "arithmetic mean across accepted libraries of basewise coverage "
            "multiplied by 1e6/effective_library_size"
        ),
        "normalization_method": TMM_BACKGROUND_METHOD,
        "library_n": len(library_ids),
        "libraries": {
            library_id: {
                **factors[library_id],
                "unit_bed": str(unit_paths[library_id].resolve()),
                "unit_bed_sha256": sha256_file(unit_paths[library_id]),
            }
            for library_id in library_ids
        },
        "inputs": {
            "normalization_factors_sha256": sha256_file(factor_path),
            "chrom_sizes_sha256": sha256_file(chrom_sizes_path),
        },
        "output": {
            "path": str(output_bigwig.resolve()),
            "sha256": sha256_file(output_bigwig),
        },
    }
    write_json_if_changed(output_metrics, metrics)
    return metrics
