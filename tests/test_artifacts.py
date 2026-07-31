import json
from pathlib import Path

import pytest

from short_read_processing.accessions import AcquisitionError
from short_read_processing.artifacts import (
    CATALOG_FILE_FIELDS,
    FINAL_BAM_FILTERING_CONTRACT,
    read_catalog_manifest,
    read_final_bam_manifest,
    read_master_manifest,
    semantic_sha256,
    sha256_file,
)


def _digest(character: str) -> str:
    return character * 64


def test_final_bam_manifest_resolves_relative_paths(tmp_path):
    bam = tmp_path / "inputs" / "rep1.bam"
    bai = tmp_path / "inputs" / "rep1.bam.bai"
    bam.parent.mkdir()
    bam.write_bytes(b"bam")
    bai.write_bytes(b"bai")
    manifest = tmp_path / "final-bams.tsv"
    manifest.write_text(
        "library_id\tassay\tcontext\trole\tlayout\tbam\tbai\tgenome\t"
        "filtering_contract\tbam_sha256\tbai_sha256\tqc_status\n"
        f"rep1\tatac\teye\ttreatment\tpaired\tinputs/rep1.bam\t"
        f"inputs/rep1.bam.bai\tdm6\t{FINAL_BAM_FILTERING_CONTRACT}\t"
        f"{_digest('a')}\t{_digest('b')}\taccepted\n"
    )

    rows = read_final_bam_manifest(manifest)

    assert rows["rep1"]["bam"] == str(bam.resolve())
    assert rows["rep1"]["bai"] == str(bai.resolve())
    assert rows["rep1"]["assay"] == "atac"


@pytest.mark.parametrize(
    ("replacement", "message"),
    [
        ("\trejected\n", "qc_status"),
        ("\twrong-contract\t", "filtering_contract"),
        ("\tshort\t", "SHA-256"),
    ],
)
def test_final_bam_manifest_rejects_invalid_contract_fields(
    tmp_path, replacement, message
):
    bam = tmp_path / "rep1.bam"
    bai = tmp_path / "rep1.bam.bai"
    bam.write_bytes(b"bam")
    bai.write_bytes(b"bai")
    row = (
        f"rep1\tatac\teye\ttreatment\tpaired\t{bam}\t{bai}\tdm6\t"
        f"{FINAL_BAM_FILTERING_CONTRACT}\t{_digest('a')}\t{_digest('b')}\taccepted\n"
    )
    if replacement == "\trejected\n":
        row = row.replace("\taccepted\n", replacement)
    elif replacement == "\twrong-contract\t":
        row = row.replace(f"\t{FINAL_BAM_FILTERING_CONTRACT}\t", replacement)
    else:
        row = row.replace(f"\t{_digest('a')}\t", replacement)
    manifest = tmp_path / "final-bams.tsv"
    manifest.write_text(
        "library_id\tassay\tcontext\trole\tlayout\tbam\tbai\tgenome\t"
        "filtering_contract\tbam_sha256\tbai_sha256\tqc_status\n" + row
    )

    with pytest.raises(AcquisitionError, match=message):
        read_final_bam_manifest(manifest)


def test_master_manifest_requires_one_complete_row(tmp_path):
    files = {}
    for name in (
        "master_bed",
        "summits_bed",
        "membership_tsv",
        "context_matrix_tsv",
        "stats_json",
    ):
        path = tmp_path / name
        path.write_text(name)
        files[name] = path
    columns = [
        "genome",
        "method",
        "input_filtering_contract",
        "source_project",
        "source_run_id",
    ]
    columns.extend(
        item
        for name in files
        for item in (name, f"{name}_sha256")
    )
    values = [
        "dm6",
        "method-v1",
        FINAL_BAM_FILTERING_CONTRACT,
        "atlas",
        "run-v1",
    ]
    values.extend(
        item
        for name, path in files.items()
        for item in (path.name, sha256_file(path))
    )
    manifest = tmp_path / "master.tsv"
    manifest.write_text("\t".join(columns) + "\n" + "\t".join(values) + "\n")

    parsed = read_master_manifest(manifest)

    assert parsed["master_bed"] == str(files["master_bed"].resolve())
    assert parsed["method"] == "method-v1"

    manifest.write_text(
        manifest.read_text().replace(
            FINAL_BAM_FILTERING_CONTRACT,
            "short-read-processing-final-v1",
        )
    )
    with pytest.raises(AcquisitionError, match="input_filtering_contract"):
        read_master_manifest(manifest)


def _write_catalog_manifest(tmp_path: Path) -> tuple[Path, dict[str, Path]]:
    catalog = tmp_path / "master_elements_long.tsv.gz"
    catalog.write_bytes(b"catalog")
    catalog_digest = sha256_file(catalog)
    metrics = tmp_path / "metrics.json"
    metrics.write_text(
        json.dumps(
            {
                "catalog_sha256": catalog_digest,
                "context_count": 2,
                "method": "catalog-v1",
            }
        )
    )
    provenance = tmp_path / "provenance.json"
    provenance.write_text(
        json.dumps({"outputs": {"catalog": {"sha256": catalog_digest}}})
    )
    resolved = tmp_path / "resolved.json"
    resolved.write_text(
        json.dumps(
            {
                "project": "atlas",
                "run_id": "catalog-v1",
                "reference": {"name": "dm6"},
                "activity": {"contexts": ["eye", "wing"]},
                "provenance": {
                    "semantic_sha256": _digest("c"),
                    "sample_sheet_sha256": _digest("d"),
                },
            }
        )
    )
    files = {
        "catalog": catalog,
        "metrics": metrics,
        "provenance": provenance,
        "resolved_config": resolved,
    }
    columns = [
        "genome",
        "method",
        "contexts",
        "source_project",
        "source_run_id",
    ]
    columns.extend(
        item for field in CATALOG_FILE_FIELDS for item in (field, f"{field}_sha256")
    )
    values = ["dm6", "catalog-v1", "eye,wing", "atlas", "catalog-v1"]
    values.extend(
        item
        for field, path in files.items()
        for item in (path.name, sha256_file(path))
    )
    manifest = tmp_path / "catalog.tsv"
    manifest.write_text("\t".join(columns) + "\n" + "\t".join(values) + "\n")
    return manifest, files


def test_catalog_manifest_verifies_bundle_and_source_metadata(tmp_path):
    manifest, files = _write_catalog_manifest(tmp_path)

    parsed = read_catalog_manifest(manifest)

    assert parsed["catalog"] == str(files["catalog"].resolve())
    assert parsed["contexts"] == ["eye", "wing"]
    assert parsed["source_semantic_sha256"] == _digest("c")


def test_catalog_manifest_rejects_changed_artifact(tmp_path):
    manifest, files = _write_catalog_manifest(tmp_path)
    files["catalog"].write_bytes(b"changed")

    with pytest.raises(AcquisitionError, match="SHA-256 mismatch"):
        read_catalog_manifest(manifest)


def test_semantic_sha256_is_order_independent():
    assert semantic_sha256({"a": 1, "b": 2}) == semantic_sha256({"b": 2, "a": 1})


def test_final_bam_manifest_allows_pending_review_for_downstream_qc(tmp_path):
    bam = tmp_path / "rep1.bam"
    bai = tmp_path / "rep1.bam.bai"
    bam.write_bytes(b"bam")
    bai.write_bytes(b"bai")
    manifest = tmp_path / "final-bams.tsv"
    manifest.write_text(
        "library_id\tassay\tcontext\trole\tlayout\tbam\tbai\tgenome\t"
        "filtering_contract\tbam_sha256\tbai_sha256\tqc_status\n"
        f"rep1\tatac\teye\ttreatment\tpaired\t{bam}\t{bai}\tdm6\t"
        f"{FINAL_BAM_FILTERING_CONTRACT}\t{_digest('a')}\t{_digest('b')}\t"
        "pending_review\n"
    )

    assert read_final_bam_manifest(manifest)["rep1"]["qc_status"] == "pending_review"


def test_final_bam_manifest_rejects_rejected_qc(tmp_path):
    bam = tmp_path / "rep1.bam"
    bai = tmp_path / "rep1.bam.bai"
    bam.write_bytes(b"bam")
    bai.write_bytes(b"bai")
    manifest = tmp_path / "final-bams.tsv"
    manifest.write_text(
        "library_id\tassay\tcontext\trole\tlayout\tbam\tbai\tgenome\t"
        "filtering_contract\tbam_sha256\tbai_sha256\tqc_status\n"
        f"rep1\tatac\teye\ttreatment\tpaired\t{bam}\t{bai}\tdm6\t"
        f"{FINAL_BAM_FILTERING_CONTRACT}\t{_digest('a')}\t{_digest('b')}\t"
        "rejected\n"
    )

    with pytest.raises(AcquisitionError, match="rejected"):
        read_final_bam_manifest(manifest)


def test_final_bam_manifest_allows_documented_rejection_for_activity_review(
    tmp_path,
):
    bam = tmp_path / "rep1.bam"
    bai = tmp_path / "rep1.bam.bai"
    bam.write_bytes(b"bam")
    bai.write_bytes(b"bai")
    manifest = tmp_path / "reviewed.tsv"
    manifest.write_text(
        "library_id\tassay\tcontext\trole\tlayout\tbam\tbai\tgenome\t"
        "filtering_contract\tbam_sha256\tbai_sha256\tqc_status\t"
        "estimated_fragment_length_bp\tnotes\n"
        f"rep1\th3k27ac\teye\ttreatment\tsingle\t{bam}\t{bai}\tdm6\t"
        f"{FINAL_BAM_FILTERING_CONTRACT}\t{_digest('a')}\t{_digest('b')}\t"
        "rejected\t380\tinsufficient depth\n"
    )

    rows = read_final_bam_manifest(manifest, allow_rejected=True)

    assert rows["rep1"]["qc_status"] == "rejected"
    assert rows["rep1"]["estimated_fragment_length_bp"] == "380"
    assert rows["rep1"]["notes"] == "insufficient depth"


@pytest.mark.parametrize("value", ["0", "not-an-integer"])
def test_final_bam_manifest_rejects_invalid_estimated_fragment_length(
    tmp_path, value
):
    bam = tmp_path / "rep1.bam"
    bai = tmp_path / "rep1.bam.bai"
    bam.write_bytes(b"bam")
    bai.write_bytes(b"bai")
    manifest = tmp_path / "reviewed.tsv"
    manifest.write_text(
        "library_id\tassay\tcontext\trole\tlayout\tbam\tbai\tgenome\t"
        "filtering_contract\tbam_sha256\tbai_sha256\tqc_status\t"
        "estimated_fragment_length_bp\n"
        f"rep1\th3k27ac\teye\ttreatment\tsingle\t{bam}\t{bai}\tdm6\t"
        f"{FINAL_BAM_FILTERING_CONTRACT}\t{_digest('a')}\t{_digest('b')}\t"
        f"accepted\t{value}\n"
    )

    with pytest.raises(AcquisitionError, match="positive integer"):
        read_final_bam_manifest(manifest)
