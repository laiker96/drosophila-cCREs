"""Shared command-line plumbing for acquisition scripts."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
from typing import Sequence

from .accessions import AcquisitionError, RunPlan, normalize_accession, resolve_accessions
from .downloader import DownloadOptions, download_plans
from .manifest import write_manifest
from .sample_sheet import read_delimited_rows


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("value must be at least 1")
    return parsed


def add_download_arguments(parser: argparse.ArgumentParser) -> None:
    cpu_count = os.cpu_count() or 1
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/raw"),
        help="FASTQ root; each resolved run gets its own directory (default: data/raw)",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        help="Run manifest TSV (default: OUTPUT_DIR/download_manifest.tsv)",
    )
    parser.add_argument(
        "--backend",
        choices=("auto", "ena", "sra"),
        default="auto",
        help=(
            "Prefer ENA compressed FASTQs with per-run SRA Toolkit failover, "
            "or force ENA/SRA Toolkit (default: auto)"
        ),
    )
    parser.add_argument(
        "--resolve-jobs",
        type=positive_int,
        default=min(16, max(4, cpu_count)),
        help="Concurrent ENA metadata requests",
    )
    parser.add_argument(
        "--file-jobs",
        type=positive_int,
        default=min(8, cpu_count),
        help="FASTQ files downloaded concurrently by aria2c",
    )
    parser.add_argument(
        "--connections",
        type=positive_int,
        default=8,
        help="Segmented aria2c connections per FASTQ (maximum 16)",
    )
    parser.add_argument(
        "--checksum-retries",
        type=positive_int,
        default=3,
        help="Automatic retries after an aria2 checksum failure (default: 3)",
    )
    parser.add_argument(
        "--sra-jobs",
        type=positive_int,
        default=min(4, cpu_count),
        help="SRA fallback runs converted concurrently",
    )
    parser.add_argument(
        "--threads",
        type=positive_int,
        default=cpu_count,
        help="Total threads shared among concurrent SRA conversions",
    )
    parser.add_argument(
        "--keep-sra-cache",
        action="store_true",
        help="Keep prefetched .sra data after successful FASTQ conversion",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Resolve accessions and write a planned manifest without downloading",
    )


def read_accession_column(path: Path, column: str) -> list[str]:
    fieldnames, rows = read_delimited_rows(path)
    if column not in fieldnames:
        available = ", ".join(fieldnames) or "none"
        raise AcquisitionError(
            f"Accession column {column!r} not found in {path}; available columns: {available}"
        )
    accessions = [
        normalize_accession(row[column]) for row in rows if row.get(column) and row[column].strip()
    ]
    if not accessions:
        raise AcquisitionError(f"No accessions found in column {column!r} of {path}")
    return list(dict.fromkeys(accessions))


def _print_plan(plans: list[RunPlan]) -> None:
    ena_runs = sum(plan.backend == "ena" for plan in plans)
    sra_runs = sum(plan.backend == "sra" for plan in plans)
    file_count = sum(len(plan.files) for plan in plans)
    total_bytes = sum(item.size_bytes or 0 for plan in plans for item in plan.files)
    gib = total_bytes / 1024**3
    print(
        f"Resolved {len(plans)} run(s): {ena_runs} ENA, {sra_runs} SRA fallback; "
        f"{file_count} direct FASTQ file(s), {gib:.2f} GiB reported"
    )
    for plan in plans:
        print(
            f"  {plan.requested_accession} -> {plan.run_accession} "
            f"({plan.library_layout.lower()}, {plan.backend})"
        )


def execute_download(accessions: Sequence[str], args: argparse.Namespace) -> Path:
    output_dir = args.output_dir.resolve()
    manifest = (args.manifest or output_dir / "download_manifest.tsv").resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    plans = resolve_accessions(
        accessions,
        output_dir=output_dir,
        backend=args.backend,
        workers=args.resolve_jobs,
    )
    _print_plan(plans)
    if args.dry_run:
        write_manifest(manifest, plans)
        print(f"Dry run only; planned manifest written to {manifest}")
        return manifest

    options = DownloadOptions(
        file_jobs=args.file_jobs,
        connections=args.connections,
        checksum_retries=args.checksum_retries,
        sra_jobs=args.sra_jobs,
        threads=args.threads,
        keep_sra_cache=args.keep_sra_cache,
        ena_fallback=args.backend == "auto",
    )
    download_plans(plans, options)
    write_manifest(manifest, plans)
    print(f"Download complete; manifest written to {manifest}")
    return manifest


def cli_main(function) -> None:
    try:
        raise SystemExit(function())
    except (AcquisitionError, FileNotFoundError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(2) from exc
