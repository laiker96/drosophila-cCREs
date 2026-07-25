#!/usr/bin/env python3
"""Run accession sample sheet -> FASTQs -> resolved YAML -> Snakemake outputs."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

from short_read_processing.cli import add_download_arguments, cli_main, execute_download
from short_read_processing.configuration import generate_configs
from short_read_processing.sample_sheet import DEFAULT_SCHEMA, sample_sheet_accessions


REPO_ROOT = Path(__file__).resolve().parents[1]


def rule_threads(value: str) -> str:
    """Validate a Snakemake RULE=THREADS override."""
    rule, separator, threads = value.rpartition("=")
    if not separator or not rule or not threads.isdigit() or int(threads) < 1:
        raise argparse.ArgumentTypeError("use RULE=THREADS with a positive thread count")
    return value


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sample_sheet", type=Path, help="Canonical accession CSV/TSV")
    parser.add_argument("--schema", type=Path, default=DEFAULT_SCHEMA)
    parser.add_argument("--project", default="short-read")
    parser.add_argument("--run-id", default="baseline")
    parser.add_argument("--genome", choices=("dm6", "hg38"), default="dm6")
    parser.add_argument("--config-dir", type=Path, default=Path("configs"))
    parser.add_argument("--reference-root", type=Path, default=Path("references"))
    parser.add_argument("--snakefile", type=Path, default=REPO_ROOT / "workflow" / "Snakefile")
    parser.add_argument(
        "--workflow-profile",
        type=Path,
        default=REPO_ROOT / "profiles" / "local",
    )
    parser.add_argument(
        "--cores",
        type=int,
        help="Maximum aggregate cores (local or across submitted cluster jobs)",
    )
    parser.add_argument("--jobs", type=int, help="Maximum concurrent cluster jobs")
    parser.add_argument(
        "--max-threads",
        type=int,
        help="Maximum threads/CPUs requested by any individual rule",
    )
    parser.add_argument(
        "--set-threads",
        action="append",
        default=[],
        type=rule_threads,
        metavar="RULE=THREADS",
        help="Override one rule's thread count; repeat for additional rules",
    )
    parser.add_argument("--skip-download", action="store_true", help="Reuse --manifest")
    parser.add_argument("--download-only", action="store_true")
    parser.add_argument("--config-only", action="store_true")
    parser.add_argument("--snakemake-dry-run", action="store_true")
    parser.add_argument(
        "--from-stage",
        choices=("accessions", "final-bam", "master"),
        default="accessions",
        help=(
            "Explicit workflow starting artifact. final-bam and master are "
            "strict reuse-only modes with no upstream fallback."
        ),
    )
    parser.add_argument(
        "--final-bam-manifest",
        type=Path,
        help="Complete immutable final-BAM manifest for --from-stage final-bam",
    )
    parser.add_argument(
        "--master-manifest",
        type=Path,
        help="Immutable master-DHS bundle manifest for --from-stage master",
    )
    parser.add_argument("--atac-minimum-replicates", type=int, default=2)
    parser.add_argument("--atac-overlap-fraction", type=float, default=0.5)
    parser.add_argument(
        "--snakemake-arg",
        action="append",
        default=[],
        help="Additional single Snakemake argument; repeat as needed",
    )
    add_download_arguments(parser)
    args = parser.parse_args()
    if args.cores is not None and args.cores < 1:
        parser.error("--cores must be positive")
    if args.jobs is not None and args.jobs < 1:
        parser.error("--jobs must be positive")
    if args.max_threads is not None and args.max_threads < 1:
        parser.error("--max-threads must be positive")
    if args.download_only and args.config_only:
        parser.error("--download-only and --config-only are mutually exclusive")
    if args.from_stage == "accessions":
        if args.final_bam_manifest or args.master_manifest:
            parser.error(
                "artifact manifests require the matching --from-stage final-bam or master"
            )
    else:
        incompatible = (
            args.skip_download
            or args.download_only
            or args.dry_run
            or args.manifest is not None
        )
        if incompatible:
            parser.error(
                "--manifest, --skip-download, --download-only, and acquisition --dry-run "
                "are valid only with --from-stage accessions"
            )
        if args.from_stage == "final-bam":
            if not args.final_bam_manifest:
                parser.error("--from-stage final-bam requires --final-bam-manifest")
            if args.master_manifest:
                parser.error("--master-manifest requires --from-stage master")
        if args.from_stage == "master":
            if not args.master_manifest:
                parser.error("--from-stage master requires --master-manifest")
            if args.final_bam_manifest:
                parser.error(
                    "--from-stage master consumes the frozen master bundle, not BAMs"
                )

    sample_sheet = args.sample_sheet.resolve()
    manifest: Path | None = None
    if args.from_stage == "accessions":
        accessions = sample_sheet_accessions(
            sample_sheet, schema_path=args.schema.resolve()
        )
        manifest = (args.manifest or args.output_dir / "download_manifest.tsv").resolve()
        if args.skip_download:
            if not manifest.is_file():
                raise FileNotFoundError(f"Manifest does not exist: {manifest}")
            print(f"Reusing download manifest: {manifest}")
        else:
            manifest = execute_download(accessions, args)
            if args.dry_run:
                print("Download dry-run complete; processing was not started")
                return 0
        if args.download_only:
            return 0
    else:
        artifact = (
            args.final_bam_manifest
            if args.from_stage == "final-bam"
            else args.master_manifest
        )
        print(
            f"Strict reuse-only mode: starting from {args.from_stage}; "
            f"manifest={artifact.resolve()}"
        )
        print("FASTQ download, trimming, and alignment fallback are disabled")

    configs = generate_configs(
        manifest_path=manifest,
        sample_sheet_path=sample_sheet,
        output_dir=args.config_dir.resolve(),
        project=args.project,
        run_id=args.run_id,
        reference_root=args.reference_root,
        path_base=REPO_ROOT,
        require_fastq_files=True,
        schema_path=args.schema.resolve(),
        genome=args.genome,
        atac_minimum_replicates=args.atac_minimum_replicates,
        atac_overlap_fraction=args.atac_overlap_fraction,
        input_stage=args.from_stage,
        final_bam_manifest_path=(
            args.final_bam_manifest.resolve() if args.final_bam_manifest else None
        ),
        master_manifest_path=(
            args.master_manifest.resolve() if args.master_manifest else None
        ),
    )
    for config_path in configs:
        print(f"Resolved workflow config: {config_path}")
    if args.config_only:
        return 0

    snakemake = shutil.which("snakemake") or str(Path(sys.executable).with_name("snakemake"))
    if not Path(snakemake).is_file() and not shutil.which(snakemake):
        raise FileNotFoundError("snakemake is not available in PATH")
    for config_path in configs:
        command = [
            snakemake,
            "--snakefile",
            str(args.snakefile.resolve()),
            "--configfile",
            str(config_path),
            "--workflow-profile",
            str(args.workflow_profile.resolve()),
        ]
        if args.cores is not None:
            command.extend(["--cores", str(args.cores)])
        if args.jobs is not None:
            command.extend(["--jobs", str(args.jobs)])
        if args.max_threads is not None:
            command.extend(["--max-threads", str(args.max_threads)])
        if args.set_threads:
            command.extend(["--set-threads", *args.set_threads])
        if args.snakemake_dry_run:
            command.append("--dry-run")
        command.extend(args.snakemake_arg)
        print("Running: " + " ".join(command))
        subprocess.run(command, cwd=REPO_ROOT, check=True)
    return 0


if __name__ == "__main__":
    cli_main(main)
