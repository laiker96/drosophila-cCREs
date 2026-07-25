"""Semantic validation for resolved workflow YAML files."""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

from .accessions import AcquisitionError
from .artifacts import semantic_sha256


SAFE_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
CHECKSUM_RE = re.compile(r"^(?:md5:[0-9a-f]{32}|sha256:[0-9a-f]{64})$")
ATAC_QPOIS_FIELDS = {
    "fragment_maximum",
    "minimum_exponent",
    "maximum_exponent",
    "merge_gap",
    "minimum_length",
    "maximum_length",
}
ATAC_CONSENSUS_FIELDS = {
    "enabled",
    "conditions",
    "minimum_replicates",
    "replicate_overlap_fraction",
}
ATAC_MASTER_FIELDS = {"summit_max_distance", "minimum_summit_separation"}
FINAL_BAM_FIELDS = {
    "bam",
    "bai",
    "bam_sha256",
    "bai_sha256",
    "filtering_contract",
    "qc_status",
    "genome",
}
EXTERNAL_MASTER_FIELDS = {
    "genome",
    "method",
    "source_project",
    "source_run_id",
    "master_bed",
    "master_bed_sha256",
    "summits_bed",
    "summits_bed_sha256",
    "membership_tsv",
    "membership_tsv_sha256",
    "context_matrix_tsv",
    "context_matrix_tsv_sha256",
    "stats_json",
    "stats_json_sha256",
}
REFERENCE_FIELDS = {
    "name",
    "fasta",
    "bowtie2_index",
    "chrom_sizes",
    "blacklist_bed",
    "tss_bed",
    "autosomes_file",
    "mitochondrial_contig",
    "effective_genome_size",
    "macs3_genome_size",
}


def wildcard_regex(values: list[str]) -> str:
    return "(?:" + "|".join(re.escape(value) for value in values) + ")" if values else r"(?!)"


def aria2_checksum(source: dict[str, Any]) -> str:
    """Translate a catalog checksum into aria2's ALGORITHM=DIGEST form."""

    algorithm, digest = str(source["checksum"]).split(":", 1)
    return f"{'sha-256' if algorithm == 'sha256' else algorithm}={digest}"


def workflow_semantic_sha256(config: dict[str, Any]) -> str:
    """Return the timestamp-independent digest used to identify a scientific run."""

    semantic_input = {
        key: value for key, value in config.items() if key != "provenance"
    }
    provenance = config.get("provenance", {})
    semantic_input["provenance_inputs"] = {
        key: value
        for key, value in provenance.items()
        if key not in {"generated_at_utc", "semantic_sha256"}
    }
    return semantic_sha256(semantic_input)


def guard_result_namespace(config: dict[str, Any], resolved_config: Path) -> None:
    """Refuse to reuse a result namespace for a different scientific config."""

    if not resolved_config.is_file():
        return
    requested = config.get("provenance", {}).get("semantic_sha256")
    # Hand-authored legacy configs are outside the generated-config contract.
    # Canonical configs from run_pipeline.py always carry this digest.
    if not requested:
        return
    try:
        existing = json.loads(resolved_config.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise AcquisitionError(
            f"Cannot verify existing result namespace {resolved_config}: {error}. "
            "Use a new run_id."
        ) from error
    recorded = existing.get("provenance", {}).get("semantic_sha256")
    if not recorded:
        raise AcquisitionError(
            f"Cannot prove that existing result namespace {resolved_config.parent.parent} "
            "has the same scientific configuration. Use a new run_id."
        )
    if requested != recorded:
        raise AcquisitionError(
            f"Result namespace {resolved_config.parent.parent} belongs to semantic "
            f"configuration {recorded}, not {requested}. Use a new run_id."
        )


def _required(mapping: dict[str, Any], fields: set[str], label: str) -> None:
    missing = sorted(field for field in fields if field not in mapping or mapping[field] is None)
    if missing:
        raise AcquisitionError(f"{label} is missing: " + ", ".join(missing))


def validate_workflow_config(config: dict[str, Any]) -> None:
    """Reject inconsistent resolved configurations before Snakemake builds the DAG."""

    _required(config, {"project", "run_id", "output_dir", "assay", "reference", "samples"}, "Config")
    for field in ("project", "run_id"):
        if not SAFE_ID_RE.fullmatch(str(config[field])):
            raise AcquisitionError(f"Invalid {field}: {config[field]!r}")
    if config["assay"] not in {"atac", "chip_tf", "chip_histone"}:
        raise AcquisitionError(f"Unsupported assay: {config['assay']!r}")
    if not isinstance(config["reference"], dict):
        raise AcquisitionError("reference must be a mapping")
    _required(config["reference"], REFERENCE_FIELDS, "Reference")
    input_stage = str(config.get("input_stage", "accessions"))
    if input_stage not in {"accessions", "final-bam", "master"}:
        raise AcquisitionError(f"Unsupported input_stage: {input_stage!r}")
    provenance = config.get("provenance")
    if provenance is not None:
        if not isinstance(provenance, dict):
            raise AcquisitionError("provenance must be a mapping")
        recorded_semantic = provenance.get("semantic_sha256")
        if recorded_semantic is not None:
            if not re.fullmatch(r"[0-9a-f]{64}", str(recorded_semantic)):
                raise AcquisitionError(
                    "provenance semantic_sha256 is not a SHA-256 digest"
                )
            calculated_semantic = workflow_semantic_sha256(config)
            if recorded_semantic != calculated_semantic:
                raise AcquisitionError(
                    "provenance semantic_sha256 does not match the resolved configuration"
                )
    external_master = config.get("external_master")
    if input_stage == "master":
        if config["assay"] != "atac":
            raise AcquisitionError("Master reuse is valid only for ATAC")
        if not isinstance(external_master, dict):
            raise AcquisitionError("Master reuse requires external_master")
        _required(external_master, EXTERNAL_MASTER_FIELDS, "External master")
        if external_master["genome"] != config["reference"].get("name"):
            raise AcquisitionError("External master and reference genomes differ")
        for field in EXTERNAL_MASTER_FIELDS:
            if field.endswith("_sha256") and not re.fullmatch(
                r"[0-9a-f]{64}", str(external_master[field])
            ):
                raise AcquisitionError(f"External master {field} is not a SHA-256 digest")
    elif external_master is not None:
        raise AcquisitionError("external_master requires input_stage=master")

    qpois = config.get("atac_qpois")
    if (
        config["assay"] == "atac"
        and input_stage != "master"
        and not isinstance(qpois, dict)
    ):
        raise AcquisitionError("ATAC configurations require atac_qpois parameters")
    if qpois is not None:
        if config["assay"] != "atac" or not isinstance(qpois, dict):
            raise AcquisitionError("atac_qpois is only valid for ATAC configurations")
        _required(qpois, ATAC_QPOIS_FIELDS, "ATAC qpois")
        if (
            int(qpois["fragment_maximum"]) < 2
            or int(qpois["minimum_exponent"]) < 0
            or int(qpois["maximum_exponent"]) < int(qpois["minimum_exponent"])
            or int(qpois["merge_gap"]) < 0
            or int(qpois["minimum_length"]) < 1
            or int(qpois["maximum_length"]) < int(qpois["minimum_length"])
        ):
            raise AcquisitionError("ATAC qpois parameters are invalid")
    consensus = config.get("atac_consensus")
    if consensus is not None:
        if config["assay"] != "atac" or not isinstance(consensus, dict):
            raise AcquisitionError("atac_consensus is only valid for ATAC configurations")
        _required(consensus, ATAC_CONSENSUS_FIELDS, "ATAC consensus")
        if not isinstance(consensus["enabled"], bool):
            raise AcquisitionError("ATAC consensus enabled must be true or false")
        if (
            int(consensus["minimum_replicates"]) < 2
            or not 0 < float(consensus["replicate_overlap_fraction"]) <= 1
            or not isinstance(consensus["conditions"], list)
            or not consensus["conditions"]
        ):
            raise AcquisitionError("ATAC consensus parameters are invalid")
    master = config.get("atac_master")
    if master is not None:
        if config["assay"] != "atac" or not isinstance(master, dict):
            raise AcquisitionError("atac_master is only valid for ATAC configurations")
        _required(master, ATAC_MASTER_FIELDS, "ATAC master DHS")
        if (
            int(master["summit_max_distance"]) < 0
            or int(master["minimum_summit_separation"]) < 0
            or int(master["minimum_summit_separation"])
            > int(master["summit_max_distance"])
        ):
            raise AcquisitionError("ATAC master DHS parameters are invalid")
        if not consensus or not consensus.get("enabled"):
            raise AcquisitionError("ATAC master DHS construction requires ATAC consensus")
    if int(config["reference"]["effective_genome_size"]) < 1:
        raise AcquisitionError("effective_genome_size must be positive")
    preparation = config["reference"].get("preparation")
    if preparation is not None:
        if not isinstance(preparation, dict) or preparation.get("mode") != "download":
            raise AcquisitionError("Reference preparation mode must be 'download'")
        _required(
            preparation,
            {"mode", "fasta", "annotation", "blacklist", "autosomes"},
            "Reference preparation",
        )
        for source_name in ("fasta", "annotation", "blacklist"):
            source = preparation[source_name]
            if not isinstance(source, dict):
                raise AcquisitionError(f"Reference {source_name} source must be a mapping")
            _required(source, {"url", "checksum"}, f"Reference {source_name} source")
            if not str(source["url"]).startswith("https://"):
                raise AcquisitionError(f"Reference {source_name} URL must use HTTPS")
            if not CHECKSUM_RE.fullmatch(str(source["checksum"])):
                raise AcquisitionError(
                    f"Reference {source_name} checksum must be md5:<hex> or sha256:<hex>"
                )
        autosomes = preparation["autosomes"]
        if not isinstance(autosomes, list) or not autosomes or any(not item for item in autosomes):
            raise AcquisitionError("Reference autosomes must be a non-empty list")

    samples = config["samples"]
    if not isinstance(samples, list) or (not samples and input_stage != "master"):
        raise AcquisitionError("samples must be a non-empty list")
    if input_stage == "master" and samples:
        raise AcquisitionError("Master reuse configuration must not schedule sample processing")
    sample_ids = [str(sample.get("id", "")) for sample in samples]
    if any(not SAFE_ID_RE.fullmatch(sample_id) for sample_id in sample_ids):
        raise AcquisitionError("Every sample must have a safe non-empty id")
    if len(sample_ids) != len(set(sample_ids)):
        raise AcquisitionError("Sample IDs must be unique")
    role_by_id = {str(sample["id"]): sample.get("role") for sample in samples}
    context_by_id = {str(sample["id"]): str(sample.get("context", "")) for sample in samples}

    for sample in samples:
        sample_id = str(sample["id"])
        sample_fields = {"accessions", "context", "role", "layout", "parameters"}
        if input_stage == "accessions":
            sample_fields.add("r1")
        else:
            sample_fields.add("final_bam")
        _required(sample, sample_fields, f"Sample {sample_id}")
        if not SAFE_ID_RE.fullmatch(str(sample["context"])):
            raise AcquisitionError(f"Sample {sample_id}: invalid context")
        if sample["role"] not in {"treatment", "control"}:
            raise AcquisitionError(f"Sample {sample_id}: invalid role")
        if sample["layout"] not in {"single", "paired"}:
            raise AcquisitionError(f"Sample {sample_id}: invalid layout")
        if input_stage == "accessions":
            if not isinstance(sample["r1"], list) or not sample["r1"]:
                raise AcquisitionError(
                    f"Sample {sample_id}: r1 must contain at least one FASTQ"
                )
            if sample["layout"] == "paired":
                if not isinstance(sample.get("r2"), list) or len(sample["r1"]) != len(
                    sample["r2"]
                ):
                    raise AcquisitionError(
                        f"Sample {sample_id}: paired r1/r2 lane counts differ"
                    )
            elif sample.get("r2"):
                raise AcquisitionError(
                    f"Sample {sample_id}: single-end input must not contain r2"
                )
        else:
            if sample.get("r1") or sample.get("r2"):
                raise AcquisitionError(
                    f"Sample {sample_id}: final-BAM mode must not contain FASTQs"
                )
            final_bam = sample["final_bam"]
            if not isinstance(final_bam, dict):
                raise AcquisitionError(
                    f"Sample {sample_id}: final_bam must be a mapping"
                )
            _required(final_bam, FINAL_BAM_FIELDS, f"Sample {sample_id} final BAM")
            if final_bam["genome"] != config["reference"]["name"]:
                raise AcquisitionError(
                    f"Sample {sample_id}: final BAM and reference genomes differ"
                )
            if final_bam["qc_status"] != "accepted":
                raise AcquisitionError(
                    f"Sample {sample_id}: final BAM QC is not accepted"
                )
            if final_bam["filtering_contract"] != "short-read-processing-final-v1":
                raise AcquisitionError(
                    f"Sample {sample_id}: unsupported final BAM filtering contract"
                )
            for field in ("bam_sha256", "bai_sha256"):
                if not re.fullmatch(r"[0-9a-f]{64}", str(final_bam[field])):
                    raise AcquisitionError(
                        f"Sample {sample_id}: {field} is not a SHA-256 digest"
                    )

        if sample["role"] == "treatment":
            peak = sample.get("peak_caller")
            if not isinstance(peak, dict) or peak.get("command") not in {"callpeak", "hmmratac"}:
                raise AcquisitionError(f"Sample {sample_id}: treatment requires a peak caller")
            if peak["command"] == "hmmratac":
                if config["assay"] != "atac" or sample["layout"] != "paired":
                    raise AcquisitionError(
                        f"Sample {sample_id}: HMMRATAC requires paired-end ATAC-seq"
                    )
            else:
                _required(
                    peak,
                    {"format", "qvalue", "broad", "nomodel", "write_bedgraph", "spmr"},
                    f"Sample {sample_id} callpeak",
                )
                if config["assay"] == "atac":
                    if (
                        peak.get("mode") != "tn5_qpois"
                        or peak["format"] != "BED"
                        or peak["broad"]
                        or not peak["nomodel"]
                        or not peak["write_bedgraph"]
                        or peak["spmr"]
                        or peak.get("shift") is None
                        or peak.get("extsize") is None
                    ):
                        raise AcquisitionError(
                            f"Sample {sample_id}: invalid two-ended Tn5 qpois configuration"
                        )
                elif not peak["write_bedgraph"] or not peak["spmr"]:
                    raise AcquisitionError(
                        f"Sample {sample_id}: ChIP callpeak must write -B --SPMR bedGraphs"
                    )

        if config["assay"].startswith("chip") and sample["role"] == "treatment":
            control = str(sample.get("control") or "")
            if control and role_by_id.get(control) != "control":
                raise AcquisitionError(f"Sample {sample_id}: invalid matched ChIP control")
            if control and context_by_id[control] != str(sample["context"]):
                raise AcquisitionError(
                    f"Sample {sample_id}: treatment and control contexts differ"
                )
        if config["assay"] == "atac" and sample["role"] != "treatment":
            raise AcquisitionError("ATAC configurations cannot contain control samples")

    if consensus and consensus["enabled"] and config["assay"] != "atac":
        raise AcquisitionError("ATAC consensus requires an ATAC configuration")
    if consensus and consensus["enabled"]:
        from .consensus import condition_specs

        specifications = condition_specs(
            consensus["conditions"],
            sample_ids=[
                sample_id
                for sample_id in sample_ids
                if role_by_id[sample_id] == "treatment"
            ],
            minimum_replicates=int(consensus["minimum_replicates"]),
        )
        for specification in specifications:
            contexts = {
                context_by_id[sample_id] for sample_id in specification.samples
            }
            if contexts != {specification.condition_id}:
                raise AcquisitionError(
                    f"ATAC condition {specification.condition_id!r} does not match "
                    "its sample contexts"
                )


def resolve_input_paths(config: dict[str, Any], base: Path) -> None:
    """Resolve relative FASTQ/reference paths against a config's launch directory in place."""

    reference = config["reference"]
    for key in ("fasta", "bowtie2_index", "chrom_sizes", "blacklist_bed", "tss_bed", "autosomes_file"):
        path = Path(reference[key])
        reference[key] = str(path if path.is_absolute() else (base / path).resolve())
    for sample in config["samples"]:
        for key in ("r1", "r2"):
            if key in sample:
                sample[key] = [
                    str(path if path.is_absolute() else (base / path).resolve())
                    for value in sample[key]
                    for path in [Path(value)]
                ]
        adapter = sample["parameters"]["trimming"].get("adapter_fasta")
        if adapter:
            path = Path(adapter)
            sample["parameters"]["trimming"]["adapter_fasta"] = str(
                path if path.is_absolute() else (base / path).resolve()
            )
        final_bam = sample.get("final_bam")
        if final_bam:
            for key in ("bam", "bai"):
                path = Path(final_bam[key])
                final_bam[key] = str(
                    path if path.is_absolute() else (base / path).resolve()
                )
    external_master = config.get("external_master")
    if external_master:
        for key in (
            "master_bed",
            "summits_bed",
            "membership_tsv",
            "context_matrix_tsv",
            "stats_json",
        ):
            path = Path(external_master[key])
            external_master[key] = str(
                path if path.is_absolute() else (base / path).resolve()
            )
