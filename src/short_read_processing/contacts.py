"""Download-independent normalization and distance modelling for contact maps."""

from __future__ import annotations

from contextlib import contextmanager
import gzip
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
from typing import Any, Iterator

import numpy as np

from .artifacts import sha256_file
from .contact_metadata import read_contact_source_manifest, verify_reported_checksum


def _cooler():
    try:
        import cooler  # type: ignore
    except ImportError as error:
        raise RuntimeError("cooler is required for contact-map processing") from error
    return cooler


def canonical_chromosome(value: str) -> str:
    return value if value.startswith("chr") else f"chr{value}"


@contextmanager
def _atomic_output(path: Path, *, suffix: str = "") -> Iterator[Path]:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=suffix, dir=path.parent
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    temporary.unlink(missing_ok=True)
    try:
        yield temporary
        if not temporary.is_file() or temporary.stat().st_size == 0:
            raise RuntimeError(f"Atomic output was not created: {path}")
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def write_json_atomic(path: Path, value: dict[str, Any]) -> None:
    with _atomic_output(path, suffix=".json") as temporary:
        temporary.write_text(
            json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )


def mcool_resolutions(path: Path) -> list[int]:
    cooler = _cooler()
    resolutions = []
    for group in cooler.fileops.list_coolers(str(path)):
        token = group.rstrip("/").split("/")[-1]
        if token.isdigit():
            resolutions.append(int(token))
    return sorted(set(resolutions))


def select_source_resolution(resolutions: list[int], target: int) -> int:
    """Choose the closest stored resolution that exactly divides the target."""

    compatible = [
        value for value in resolutions if value <= target and target % value == 0
    ]
    if not compatible:
        raise ValueError(
            f"No source resolution in {resolutions} evenly divides target {target}"
        )
    return max(compatible)


def _decompress_cool(path: Path, destination: Path) -> Path:
    if destination.is_file() and destination.stat().st_size > 0:
        return destination
    with _atomic_output(destination, suffix=".cool") as temporary:
        with gzip.open(path, "rb") as source, temporary.open("wb") as target:
            shutil.copyfileobj(source, target, length=8 * 1024 * 1024)
    return destination


def _convert_h5(path: Path, destination: Path) -> Path:
    if destination.is_file() and destination.stat().st_size > 0:
        return destination
    converter = shutil.which("hicConvertFormat")
    if converter is None:
        raise RuntimeError("HiCExplorer hicConvertFormat is required for H5 contacts")
    with _atomic_output(destination, suffix=".cool") as temporary:
        subprocess.run(
            [
                converter,
                "--matrices",
                str(path),
                "--outFileName",
                str(temporary),
                "--inputFormat",
                "h5",
                "--outputFormat",
                "cool",
            ],
            check=True,
        )
    return destination


def _coarsen_if_needed(uri: str, target: int, destination: Path) -> str:
    cooler = _cooler()
    source = cooler.Cooler(uri)
    if source.binsize is None:
        raise ValueError(f"Variable-bin Cooler is unsupported: {uri}")
    source_resolution = int(source.binsize)
    if source_resolution == target:
        return uri
    if source_resolution > target or target % source_resolution:
        raise ValueError(
            f"Source resolution {source_resolution} does not divide target {target}: {uri}"
        )
    if destination.is_file() and destination.stat().st_size > 0:
        existing = cooler.Cooler(str(destination))
        if existing.binsize != target:
            raise ValueError(f"Cached contact has the wrong resolution: {destination}")
        return str(destination)
    with _atomic_output(destination, suffix=".cool") as temporary:
        cooler.coarsen_cooler(
            uri,
            str(temporary),
            factor=target // source_resolution,
            chunksize=20_000_000,
        )
    return str(destination)


def _prepare_source(
    source_id: str,
    path: Path,
    file_format: str,
    target_resolution: int,
    workdir: Path,
) -> tuple[str, int]:
    cooler = _cooler()
    if file_format == "mcool":
        resolutions = mcool_resolutions(path)
        selected = select_source_resolution(resolutions, target_resolution)
        uri = f"{path}::/resolutions/{selected}"
    elif file_format == "cool.gz":
        decompressed = _decompress_cool(path, workdir / f"{source_id}.cool")
        uri = str(decompressed)
        selected = int(cooler.Cooler(uri).binsize or 0)
    elif file_format == "h5":
        converted = _convert_h5(path, workdir / f"{source_id}.converted.cool")
        uri = str(converted)
        selected = int(cooler.Cooler(uri).binsize or 0)
    else:
        raise ValueError(f"Unsupported contact format: {file_format}")
    prepared = _coarsen_if_needed(
        uri,
        target_resolution,
        workdir / f"{source_id}.{target_resolution}.cool",
    )
    return prepared, selected


def standardize_context(
    *,
    context: str,
    source_manifest: Path,
    repository_root: Path,
    target_resolution: int,
    workdir: Path,
    output: Path,
) -> dict[str, Any]:
    """Coarsen exact divisors, sum replicate counts, and ICE-balance the merge."""

    cooler = _cooler()
    rows = [
        row for row in read_contact_source_manifest(source_manifest)
        if row["context"] == context
    ]
    if not rows:
        raise ValueError(f"No primary contact sources for {context}")
    if target_resolution < 1:
        raise ValueError("Contact resolution must be positive")

    workdir.mkdir(parents=True, exist_ok=True)
    prepared_uris: list[str] = []
    source_resolutions: list[int] = []
    source_paths: list[Path] = []
    intermediate_paths: list[Path] = []
    for row in rows:
        source = (repository_root / row["local_path"]).resolve()
        if not source.is_file() or source.stat().st_size == 0:
            raise FileNotFoundError(source)
        verify_reported_checksum(source, row["checksum"])
        uri, selected = _prepare_source(
            row["source_id"],
            source,
            row["format"],
            target_resolution,
            workdir,
        )
        prepared_uris.append(uri)
        source_resolutions.append(selected)
        source_paths.append(source)
        intermediate_paths.extend(
            [
                workdir / f"{row['source_id']}.cool",
                workdir / f"{row['source_id']}.converted.cool",
                workdir / f"{row['source_id']}.{target_resolution}.cool",
            ]
        )

    with _atomic_output(output, suffix=".cool") as temporary:
        cooler.merge_coolers(
            str(temporary),
            prepared_uris,
            mergebuf=20_000_000,
            columns=["count"],
        )
        matrix = cooler.Cooler(str(temporary))
        weights, balance_stats = cooler.balance_cooler(
            matrix,
            ignore_diags=2,
            mad_max=5,
            min_nnz=10,
            tol=1e-5,
            max_iters=200,
            chunksize=10_000_000,
            store=True,
        )
        if not np.isfinite(weights).any():
            raise ValueError(f"Balancing produced no finite weights for {context}")
        information = dict(matrix.info)

    metrics = {
        "schema_version": 1,
        "context": context,
        "contact_strategy": "observed",
        "contact_match": rows[0]["match_quality"],
        "assays": sorted({row["assay"] for row in rows}),
        "source_ids": [row["source_id"] for row in rows],
        "source_files": [
            {
                "path": str(path),
                "sha256": sha256_file(path),
                "url": row["url"],
                "reported_checksum": row["checksum"] or None,
            }
            for row, path in zip(rows, source_paths)
        ],
        "source_resolutions_bp": source_resolutions,
        "target_resolution_bp": target_resolution,
        "replicate_count": len(rows),
        "normalization": "select exact resolution; coarsen by count summation; "
        "sum replicate counts; ICE-balance merged matrix",
        "output": str(output),
        "output_sha256": sha256_file(output),
        "bins": int(information["nbins"]),
        "pixels": int(information["nnz"]),
        "finite_weight_fraction": float(np.isfinite(weights).mean()),
        "balance_converged": bool(balance_stats.get("converged", False)),
        "caveats": sorted({row["caveat"] for row in rows if row["caveat"]}),
        "intermediate_cleanup": "converted and coarsened work copies removed after success",
    }
    for path in intermediate_paths:
        path.unlink(missing_ok=True)
    try:
        workdir.rmdir()
    except OSError:
        pass
    return metrics


def fit_matrix_powerlaw(
    uri: str,
    canonical_chromosomes: list[str],
    maximum_distance: int,
) -> dict[str, float | int]:
    """Fit mean finite balanced cis contact, including zero pixels, by distance."""

    cooler = _cooler()
    matrix = cooler.Cooler(uri)
    if matrix.binsize is None:
        raise ValueError(f"Variable-bin Cooler is unsupported: {uri}")
    resolution = int(matrix.binsize)
    aliases = {
        canonical_chromosome(chromosome): chromosome
        for chromosome in matrix.chromnames
    }
    maximum_offset = max(2, maximum_distance // resolution)
    offsets = np.unique(np.geomspace(1, maximum_offset, 40).astype(int))
    contact_sums = np.zeros(len(offsets), dtype=float)
    contact_counts = np.zeros(len(offsets), dtype=np.int64)
    for chromosome in canonical_chromosomes:
        source_chromosome = aliases.get(canonical_chromosome(chromosome))
        if source_chromosome is not None:
            chromosome_matrix = (
                matrix.matrix(balance=True, sparse=True)
                .fetch(source_chromosome)
                .tocsr()
            )
            weights = np.asarray(
                matrix.bins().fetch(source_chromosome)["weight"], dtype=float
            )
            for index, offset in enumerate(offsets):
                valid = valid_balanced_diagonal(
                    chromosome_matrix, weights, int(offset)
                )
                contact_sums[index] += float(np.sum(valid))
                contact_counts[index] += len(valid)
    distances: list[float] = []
    contacts: list[float] = []
    for offset, contact_sum, contact_count in zip(
        offsets, contact_sums, contact_counts
    ):
        mean_contact = contact_sum / contact_count if contact_count else 0.0
        if mean_contact > 0:
            distances.append(float(offset * resolution))
            contacts.append(mean_contact)
    if len(distances) < 4:
        raise ValueError(f"Too few positive diagonals to fit contact decay: {uri}")
    gamma, intercept = np.polyfit(np.log(distances), np.log(contacts), 1)
    if not -3.0 < gamma < -0.05:
        raise ValueError(f"Implausible contact-decay exponent {gamma:.5g}: {uri}")
    return {
        "gamma": float(gamma),
        "scale": float(math.exp(intercept)),
        "resolution_bp": resolution,
        "fitted_diagonals": len(distances),
    }


def valid_balanced_diagonal(
    chromosome_matrix: Any,
    bin_weights: np.ndarray,
    offset: int,
) -> np.ndarray:
    """Return finite nonnegative pixels whose two bins both passed balancing."""

    if offset < 1:
        raise ValueError("Contact diagonal offset must be positive")
    if offset >= len(bin_weights):
        return np.asarray([], dtype=float)
    diagonal = np.asarray(chromosome_matrix.diagonal(offset), dtype=float).reshape(-1)
    valid_bins = np.isfinite(bin_weights[:-offset]) & np.isfinite(
        bin_weights[offset:]
    )
    if len(diagonal) != len(valid_bins):
        raise ValueError("Contact matrix and bin weights have inconsistent dimensions")
    return diagonal[valid_bins & np.isfinite(diagonal) & (diagonal >= 0)]


def fit_atlas_powerlaw(
    *,
    contacts: dict[str, Path],
    canonical_chromosomes: list[str],
    maximum_distance: int,
    output: Path,
) -> dict[str, Any]:
    """Fit each observed context and an atlas-wide model for missing contexts."""

    if not contacts:
        raise ValueError("At least one observed contact map is required")
    fits = {
        context: fit_matrix_powerlaw(
            str(path), canonical_chromosomes, maximum_distance
        )
        for context, path in sorted(contacts.items())
    }
    gamma = float(np.mean([float(fit["gamma"]) for fit in fits.values()]))
    scale = float(
        math.exp(np.mean([math.log(float(fit["scale"])) for fit in fits.values()]))
    )
    result = {
        "schema_version": 1,
        "contexts": fits,
        "atlas_powerlaw": {
            "gamma": gamma,
            "scale": scale,
            "source": (
                "mean gamma and geometric-mean scale fitted from observed atlas contacts"
            ),
        },
        "maximum_distance_bp": maximum_distance,
        "contact_sha256": {
            context: sha256_file(path) for context, path in sorted(contacts.items())
        },
    }
    write_json_atomic(output, result)
    return result
