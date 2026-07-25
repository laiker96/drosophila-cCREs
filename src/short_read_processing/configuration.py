"""Generate accession-first pipeline YAML from a validated sample sheet and manifest."""

from __future__ import annotations

from collections import defaultdict
from datetime import datetime, timezone
import os
from pathlib import Path
import re
import tempfile
from typing import Any

import yaml

from .accessions import AcquisitionError
from .artifacts import (
    read_final_bam_manifest,
    read_master_manifest,
    sha256_file,
)
from .manifest import read_manifest
from .sample_sheet import DEFAULT_SCHEMA, read_sample_sheet
from .workflow_config import workflow_semantic_sha256


GENOME_DEFAULTS = {
    "dm6": {"effective_genome_size": 142_573_017, "macs3_genome_size": "dm"},
    "hg38": {"effective_genome_size": 2_913_022_398, "macs3_genome_size": "hs"},
}
REFERENCE_SOURCES = {
    "dm6": {
        "fasta": {
            "url": "https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz",
            "checksum": "md5:4a5777324403eff92c3650c30af30120",
        },
        "annotation": {
            "url": (
                "https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/genes/"
                "dm6.ncbiRefSeq.gtf.gz"
            ),
            "checksum": "md5:17fa021bf9cf35a12703787cc18ea4c8",
        },
        "blacklist": {
            "url": (
                "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/v2.0/lists/"
                "dm6-blacklist.v2.bed.gz"
            ),
            "checksum": (
                "sha256:6174aa0304859a2ee836de7f1eaabc4e7ab7fc4c75ffd4b7587172b1e234026e"
            ),
        },
        "autosomes": ["chr2L", "chr2R", "chr3L", "chr3R", "chr4"],
    },
    "hg38": {
        "fasta": {
            "url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz",
            "checksum": "md5:efb4edf237dd3594e94610ed803c8a44",
        },
        "annotation": {
            "url": (
                "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/"
                "hg38.ncbiRefSeq.gtf.gz"
            ),
            "checksum": "md5:15e6cab3a16e87f2564839f7500700e2",
        },
        "blacklist": {
            "url": (
                "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/v2.0/lists/"
                "hg38-blacklist.v2.bed.gz"
            ),
            "checksum": (
                "sha256:c92e763af17271446194991e71917ac220593a5a3d40a06667be24178ef08cf2"
            ),
        },
        "autosomes": [f"chr{chromosome}" for chromosome in range(1, 23)],
    },
}
SAFE_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")

ATAC_QPOIS_DEFAULTS = {
    "fragment_maximum": 150,
    "minimum_exponent": 2,
    "maximum_exponent": 325,
    "merge_gap": 1,
    "minimum_length": 50,
    "maximum_length": 400,
}
ATAC_MASTER_DEFAULTS = {
    "summit_max_distance": 150,
    "minimum_summit_separation": 50,
}


def _safe_id(value: str, label: str) -> str:
    if not SAFE_ID_RE.fullmatch(value):
        raise AcquisitionError(
            f"Invalid {label} {value!r}; use letters, digits, periods, underscores, or hyphens"
        )
    return value


def _display_path(value: str | Path, path_base: Path) -> str:
    path = Path(value).resolve()
    try:
        return path.relative_to(path_base.resolve()).as_posix()
    except ValueError:
        return str(path)


def _split_extra_fastqs(value: str) -> list[str]:
    return [item for item in value.split(";") if item]


def _sample_fastqs(
    manifest_rows: list[dict[str, str]],
    *,
    path_base: Path,
    require_files: bool,
) -> tuple[str, list[str], list[str]]:
    layouts = {row["library_layout"].upper() for row in manifest_rows}
    if len(layouts) != 1 or layouts.pop() not in {"SINGLE", "PAIRED"}:
        raise AcquisitionError("All runs for one sample must have the same SINGLE or PAIRED layout")
    layout = manifest_rows[0]["library_layout"].upper()

    r1 = [row["fastq_1"] for row in manifest_rows if row["fastq_1"]]
    r2 = [row["fastq_2"] for row in manifest_rows if row["fastq_2"]]
    if not r1 or (layout == "PAIRED" and len(r1) != len(r2)):
        raise AcquisitionError("Manifest does not contain a complete FASTQ set for the sample")
    paths = r1 + r2 + [
        item for row in manifest_rows for item in _split_extra_fastqs(row["extra_fastqs"])
    ]
    if require_files:
        missing = [item for item in paths if not Path(item).is_file()]
        if missing:
            raise AcquisitionError("Manifest FASTQ files are missing:\n  " + "\n  ".join(missing))
    return (
        "paired" if layout == "PAIRED" else "single",
        [_display_path(item, path_base) for item in r1],
        [_display_path(item, path_base) for item in r2],
    )


def _reference_config(genome: str, reference_root: Path, path_base: Path) -> dict[str, object]:
    root = reference_root / genome
    defaults = GENOME_DEFAULTS[genome]
    return {
        "name": genome,
        "fasta": _display_path(root / f"{genome}.fa", path_base),
        "bowtie2_index": _display_path(root / "bowtie2" / genome, path_base),
        "chrom_sizes": _display_path(root / f"{genome}.chrom.sizes", path_base),
        "blacklist_bed": _display_path(root / f"{genome}.blacklist.bed", path_base),
        "tss_bed": _display_path(root / f"{genome}.tss.bed", path_base),
        "autosomes_file": _display_path(root / f"{genome}.autosomes.txt", path_base),
        "mitochondrial_contig": "chrM",
        "effective_genome_size": defaults["effective_genome_size"],
        "macs3_genome_size": defaults["macs3_genome_size"],
        "preparation": {
            "mode": "download",
            **REFERENCE_SOURCES[genome],
        },
    }


def _peak_config(row: dict[str, Any], layout: str) -> dict[str, object]:
    assay = str(row["assay"])
    caller = str(row["peak_caller"])
    if caller == "hmmratac":
        if assay != "atac":
            raise AcquisitionError("HMMRATAC is only valid for ATAC-seq")
        if layout != "paired":
            raise AcquisitionError(
                f"Library {row['library_id']}: HMMRATAC requires paired-end ATAC-seq; "
                "set peak_caller=callpeak for single-end data"
            )
        return {
            "command": "hmmratac",
            "lower": row["hmmratac_lower"],
            "upper": row["hmmratac_upper"],
            "prescan_cutoff": row["hmmratac_prescan_cutoff"],
        }

    if assay == "atac":
        if row["macs3_broad"]:
            raise AcquisitionError("ATAC two-ended qpois calling does not support broad peaks")
        return {
            "command": "callpeak",
            "mode": "tn5_qpois",
            "format": "BED",
            "qvalue": row["macs3_qvalue"],
            "broad": False,
            "nomodel": True,
            "shift": row["macs3_shift"] if row["macs3_shift"] is not None else -75,
            "extsize": row["macs3_extsize"] if row["macs3_extsize"] is not None else 150,
            "write_bedgraph": True,
            "spmr": False,
        }

    macs_format = str(row["macs3_format"])
    if macs_format == "AUTO":
        macs_format = "BAMPE" if layout == "paired" else "BAM"
    config: dict[str, object] = {
        "command": "callpeak",
        "format": macs_format,
        "qvalue": row["macs3_qvalue"],
        "broad": row["macs3_broad"],
        "nomodel": row["macs3_nomodel"],
        "shift": row["macs3_shift"],
        "extsize": row["macs3_extsize"],
        "write_bedgraph": True,
        "spmr": True,
        "bedgraph_outputs": {
            "treatment_suffix": "_treat_pileup.bdg",
            "control_suffix": "_control_lambda.bdg",
        },
    }
    if row["macs3_broad"]:
        config["broad_cutoff"] = row["macs3_broad_cutoff"]
    return config


def _processing_config(row: dict[str, Any], path_base: Path) -> dict[str, object]:
    adapter_fasta = row.get("adapter_fasta")
    return {
        "trimming": {
            "adapter_preset": row["adapter_preset"],
            "adapter_fasta": _display_path(adapter_fasta, path_base) if adapter_fasta else None,
            "quality_cutoff": row["trim_quality_cutoff"],
            "minimum_length": row["trim_minimum_length"],
            "error_rate": 0.1,
            "minimum_overlap": 3,
        },
        "alignment": {
            "preset": row["bowtie2_preset"],
            "maximum_fragment_length": row["maximum_fragment_length"],
            "mapq_minimum": row["mapq_minimum"],
        },
        "filtering": {
            "remove_duplicates": row["remove_duplicates"],
            "remove_mitochondrial": row["remove_mitochondrial"],
        },
    }


def _manifest_runs_for_accessions(
    accessions: list[str], manifest_by_request: dict[str, list[dict[str, str]]]
) -> list[dict[str, str]]:
    by_run: dict[str, dict[str, str]] = {}
    for accession in accessions:
        runs = manifest_by_request.get(accession, [])
        if not runs:
            raise AcquisitionError(f"No manifest rows found for sample-sheet accession {accession}")
        for run in runs:
            by_run[run["run_accession"]] = run
    return [by_run[key] for key in sorted(by_run)]


def _write_config_if_changed(path: Path, config: dict[str, object]) -> None:
    """Atomically write a config without changing an identical file's mtime."""

    if path.is_file():
        existing = yaml.safe_load(path.read_text(encoding="utf-8"))
        if isinstance(existing, dict):
            existing_provenance = existing.get("provenance")
            provenance = config.get("provenance")
            if isinstance(existing_provenance, dict) and isinstance(provenance, dict):
                new_generated_at = provenance.get("generated_at_utc")
                generated_at = existing_provenance.get("generated_at_utc")
                if generated_at:
                    provenance["generated_at_utc"] = generated_at
                if existing == config:
                    return
                provenance["generated_at_utc"] = new_generated_at
            elif existing == config:
                return

    content = yaml.safe_dump(config, sort_keys=False, default_flow_style=False)
    temporary_name: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
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


def generate_configs(
    *,
    manifest_path: Path | None,
    sample_sheet_path: Path,
    output_dir: Path,
    project: str,
    run_id: str,
    reference_root: Path,
    path_base: Path,
    require_fastq_files: bool,
    schema_path: Path = DEFAULT_SCHEMA,
    genome: str = "dm6",
    atac_minimum_replicates: int = 2,
    atac_overlap_fraction: float = 0.5,
    input_stage: str = "accessions",
    final_bam_manifest_path: Path | None = None,
    master_manifest_path: Path | None = None,
) -> list[Path]:
    project = _safe_id(project, "project ID")
    run_id = _safe_id(run_id, "run ID")
    if genome not in GENOME_DEFAULTS:
        raise AcquisitionError(f"Unsupported genome: {genome!r}")
    if input_stage not in {"accessions", "final-bam", "master"}:
        raise AcquisitionError(f"Unsupported input stage: {input_stage!r}")
    if input_stage == "accessions" and manifest_path is None:
        raise AcquisitionError("Accession mode requires a download manifest")
    if input_stage == "final-bam" and final_bam_manifest_path is None:
        raise AcquisitionError("Final-BAM mode requires a final-BAM manifest")
    if input_stage == "master" and master_manifest_path is None:
        raise AcquisitionError("Master mode requires a master manifest")

    sheet_rows = read_sample_sheet(sample_sheet_path, schema_path=schema_path)
    sample_sheet_sha256 = sha256_file(sample_sheet_path)
    schema_sha256 = sha256_file(schema_path)

    if input_stage == "master":
        if not any(row["assay"] == "atac" for row in sheet_rows):
            raise AcquisitionError("Master reuse requires ATAC metadata in the sample sheet")
        master = read_master_manifest(
            master_manifest_path,
            require_files=require_fastq_files,
        )
        if master["genome"] != genome:
            raise AcquisitionError(
                f"Master manifest genome {master['genome']!r} does not match {genome!r}"
            )
        assays = {str(row["assay"]) for row in sheet_rows}
        group_project = project if assays == {"atac"} else f"{project}.atac.{genome}"
        external_master = {
            key: (
                _display_path(value, path_base)
                if key
                in {
                    "master_bed",
                    "summits_bed",
                    "membership_tsv",
                    "context_matrix_tsv",
                    "stats_json",
                }
                else value
            )
            for key, value in master.items()
        }
        config: dict[str, Any] = {
            "project": group_project,
            "run_id": run_id,
            "output_dir": "results",
            "assay": "atac",
            "input_stage": "master",
            "reference": _reference_config(genome, reference_root, path_base),
            "samples": [],
            "external_master": external_master,
            "provenance": {
                "generated_at_utc": datetime.now(timezone.utc).isoformat(),
                "input_mode": "external_master_manifest",
                "sample_sheet": _display_path(sample_sheet_path, path_base),
                "sample_sheet_sha256": sample_sheet_sha256,
                "sample_sheet_schema": _display_path(schema_path, path_base),
                "sample_sheet_schema_sha256": schema_sha256,
                "master_manifest": _display_path(master_manifest_path, path_base),
                "master_manifest_sha256": sha256_file(master_manifest_path),
            },
        }
        config["provenance"]["semantic_sha256"] = workflow_semantic_sha256(config)
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / f"{group_project}.yaml"
        _write_config_if_changed(output_path, config)
        return [output_path.resolve()]

    manifest_rows = read_manifest(manifest_path) if input_stage == "accessions" else []
    manifest_by_request: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in manifest_rows:
        manifest_by_request[row["requested_accession"]].append(row)

    final_bams = (
        read_final_bam_manifest(
            final_bam_manifest_path,
            require_files=require_fastq_files,
        )
        if input_stage == "final-bam"
        else {}
    )
    if final_bams:
        sheet_library_assay = {
            str(row["library_id"]): str(row["assay"]) for row in sheet_rows
        }
        unknown = sorted(set(final_bams) - set(sheet_library_assay))
        if unknown:
            raise AcquisitionError(
                "Final-BAM manifest contains libraries absent from the sample sheet: "
                + ", ".join(unknown)
            )
        selected_assays = {
            sheet_library_assay[library_id] for library_id in final_bams
        }
        expected = {
            library_id
            for library_id, assay in sheet_library_assay.items()
            if assay in selected_assays
        }
        missing = sorted(expected - set(final_bams))
        if missing:
            raise AcquisitionError(
                "Final-BAM manifest is incomplete for its selected assay(s): "
                + ", ".join(missing)
            )
        sheet_rows = [
            row for row in sheet_rows if str(row["library_id"]) in final_bams
        ]

    sheet_by_library: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in sheet_rows:
        sheet_by_library[str(row["library_id"])].append(row)

    prepared: list[dict[str, Any]] = []
    for library_id, rows in sheet_by_library.items():
        first = rows[0]
        accessions = list(dict.fromkeys(str(row["accession"]) for row in rows))
        artifact: dict[str, str] | None = None
        if input_stage == "accessions":
            runs = _manifest_runs_for_accessions(accessions, manifest_by_request)
            layout, r1, r2 = _sample_fastqs(
                runs,
                path_base=path_base,
                require_files=require_fastq_files,
            )
        else:
            artifact = final_bams[library_id]
            for field in ("assay", "context", "role"):
                if str(artifact[field]) != str(first[field]):
                    raise AcquisitionError(
                        f"Final-BAM manifest {field} for {library_id!r} "
                        "does not match the sample sheet"
                    )
            if artifact["genome"] != genome:
                raise AcquisitionError(
                    f"Final-BAM manifest genome for {library_id!r} "
                    f"is {artifact['genome']!r}, expected {genome!r}"
                )
            layout = artifact["layout"]
            r1, r2 = [], []
        prepared.append(
            {
                "id": library_id,
                "accessions": accessions,
                "assay": first["assay"],
                "genome": genome,
                "context": first["context"],
                "role": first["role"],
                "control": first["control_library"],
                "layout": layout,
                "r1": r1,
                "r2": r2,
                "peak_caller": _peak_config(first, layout),
                "parameters": _processing_config(first, path_base),
                "artifact": artifact,
            }
        )

    groups: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for item in prepared:
        groups[(str(item["assay"]), str(item["genome"]))].append(item)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_paths: list[Path] = []
    for (assay, genome), items in sorted(groups.items()):
        sample_ids = {str(item["id"]) for item in items}
        role_by_id = {str(item["id"]): str(item["role"]) for item in items}
        item_by_id = {str(item["id"]): item for item in items}
        samples: list[dict[str, object]] = []
        for item in items:
            sample: dict[str, object] = {
                "id": item["id"],
                "accessions": item["accessions"],
                "context": item["context"],
                "role": item["role"],
                "layout": item["layout"],
                "parameters": item["parameters"],
            }
            if input_stage == "accessions":
                sample["r1"] = item["r1"]
                if item["layout"] == "paired":
                    sample["r2"] = item["r2"]
            else:
                artifact = dict(item["artifact"])
                artifact["bam"] = _display_path(artifact["bam"], path_base)
                artifact["bai"] = _display_path(artifact["bai"], path_base)
                sample["final_bam"] = artifact
            if item["role"] == "treatment":
                sample["peak_caller"] = item["peak_caller"]
            if assay.startswith("chip") and item["role"] == "treatment":
                control_library = str(item["control"] or "")
                if control_library:
                    if (
                        control_library not in sample_ids
                        or role_by_id.get(control_library) != "control"
                    ):
                        raise AcquisitionError(
                            f"ChIP treatment {item['id']} has an invalid matched control"
                        )
                    if item_by_id[control_library]["layout"] != item["layout"]:
                        raise AcquisitionError(
                            f"ChIP treatment {item['id']} and control {control_library} "
                            "must have the same read layout"
                        )
                    sample["control"] = control_library
            samples.append(sample)

        group_project = project if len(groups) == 1 else f"{project}.{assay}.{genome}"
        config = {
            "project": group_project,
            "run_id": run_id,
            "output_dir": "results",
            "assay": assay,
            "input_stage": input_stage,
            "reference": _reference_config(genome, reference_root, path_base),
            "samples": samples,
            "provenance": {
                "generated_at_utc": datetime.now(timezone.utc).isoformat(),
                "input_mode": (
                    "accession_sample_sheet"
                    if input_stage == "accessions"
                    else "external_final_bam_manifest"
                ),
                "sample_sheet": _display_path(sample_sheet_path, path_base),
                "sample_sheet_sha256": sample_sheet_sha256,
                "sample_sheet_schema": _display_path(schema_path, path_base),
                "sample_sheet_schema_sha256": schema_sha256,
            },
        }
        if input_stage == "accessions":
            config["provenance"]["download_manifest"] = _display_path(
                manifest_path, path_base
            )
            config["provenance"]["download_manifest_sha256"] = sha256_file(
                manifest_path
            )
        else:
            config["provenance"]["final_bam_manifest"] = _display_path(
                final_bam_manifest_path, path_base
            )
            config["provenance"]["final_bam_manifest_sha256"] = sha256_file(
                final_bam_manifest_path
            )
        if assay == "atac":
            from .consensus import condition_specs

            config["atac_qpois"] = dict(ATAC_QPOIS_DEFAULTS)
            samples_by_context: dict[str, list[str]] = defaultdict(list)
            for item in items:
                samples_by_context[str(item["context"])].append(str(item["id"]))
            condition_values = [
                {
                    "id": context,
                    "label": context,
                    "samples": sorted(context_samples),
                }
                for context, context_samples in sorted(samples_by_context.items())
            ]
            conditions = condition_specs(
                condition_values,
                sample_ids=sample_ids,
                minimum_replicates=atac_minimum_replicates,
            )
            for condition in conditions:
                methods = {
                    str(item_by_id[sample]["peak_caller"]["command"])
                    for sample in condition.samples
                }
                layouts = {
                    str(item_by_id[sample]["layout"])
                    for sample in condition.samples
                }
                if len(methods) != 1 or len(layouts) != 1:
                    raise AcquisitionError(
                        f"ATAC condition {condition.condition_id} mixes peak callers or layouts"
                    )
            if not 0 < atac_overlap_fraction <= 1:
                raise AcquisitionError("ATAC overlap fraction must be in (0, 1]")
            config["atac_consensus"] = {
                "enabled": True,
                "conditions": condition_values,
                "minimum_replicates": atac_minimum_replicates,
                "replicate_overlap_fraction": atac_overlap_fraction,
            }
            config["atac_master"] = dict(ATAC_MASTER_DEFAULTS)
        config["provenance"]["semantic_sha256"] = workflow_semantic_sha256(config)
        output_path = output_dir / f"{group_project}.yaml"
        _write_config_if_changed(output_path, config)
        output_paths.append(output_path.resolve())
    return output_paths
