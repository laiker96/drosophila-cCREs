import sys
from pathlib import Path

import pytest

from review_final_bam_manifest import main
from short_read_processing.accessions import AcquisitionError
from short_read_processing.artifacts import (
    FINAL_BAM_FILTERING_CONTRACT,
    read_final_bam_manifest,
    sha256_file,
)


def _source_manifest(tmp_path: Path) -> Path:
    bam = tmp_path / "input" / "h3.final.bam"
    bai = tmp_path / "input" / "h3.final.bam.bai"
    bam.parent.mkdir()
    bam.write_bytes(b"bam")
    bai.write_bytes(b"bai")
    manifest = tmp_path / "pending.tsv"
    manifest.write_text(
        "library_id\tassay\tcontext\trole\tlayout\tbam\tbai\tgenome\t"
        "filtering_contract\tbam_sha256\tbai_sha256\tqc_status\n"
        f"h3\tchip_histone\teye\ttreatment\tsingle\t{bam}\t{bai}\tdm6\t"
        f"{FINAL_BAM_FILTERING_CONTRACT}\t{sha256_file(bam)}\t"
        f"{sha256_file(bai)}\tpending_review\n"
    )
    return manifest


def test_review_manifest_records_single_end_fragment_length(tmp_path, monkeypatch):
    source = _source_manifest(tmp_path)
    decisions = tmp_path / "decisions.tsv"
    decisions.write_text(
        "library_id\tqc_status\testimated_fragment_length_bp\tnotes\n"
        "h3\taccepted\t165\taccepted after QC\n"
    )
    output = tmp_path / "reviewed.tsv"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "review_final_bam_manifest.py",
            str(source),
            "--decisions",
            str(decisions),
            "--output",
            str(output),
        ],
    )

    assert main() == 0
    row = read_final_bam_manifest(output)["h3"]
    assert row["qc_status"] == "accepted"
    assert row["estimated_fragment_length_bp"] == "165"


def test_review_manifest_requires_single_end_fragment_length(tmp_path, monkeypatch):
    source = _source_manifest(tmp_path)
    decisions = tmp_path / "decisions.tsv"
    decisions.write_text(
        "library_id\tqc_status\testimated_fragment_length_bp\tnotes\n"
        "h3\taccepted\t\taccepted after QC\n"
    )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "review_final_bam_manifest.py",
            str(source),
            "--decisions",
            str(decisions),
            "--output",
            str(tmp_path / "reviewed.tsv"),
        ],
    )

    with pytest.raises(AcquisitionError, match="requires a positive"):
        main()
