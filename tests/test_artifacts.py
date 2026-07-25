from pathlib import Path

import pytest

from short_read_processing.accessions import AcquisitionError
from short_read_processing.artifacts import (
    FINAL_BAM_FILTERING_CONTRACT,
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
    columns = ["genome", "method", "source_project", "source_run_id"]
    columns.extend(
        item
        for name in files
        for item in (name, f"{name}_sha256")
    )
    values = ["dm6", "method-v1", "atlas", "run-v1"]
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
