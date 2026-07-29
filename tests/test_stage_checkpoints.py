import json
from pathlib import Path

import pytest

from short_read_processing.accessions import AcquisitionError
from short_read_processing.artifacts import sha256_file
from short_read_processing.stage_checkpoints import read_stage_checkpoint


def _checkpoint(tmp_path: Path, *, stage: str = "alignment") -> Path:
    artifact = tmp_path / "artifact.txt"
    artifact.write_text("complete\n")
    manifest = tmp_path / f"{stage}.checkpoint.json"
    manifest.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "stage": stage,
                "source_project": "atlas",
                "source_run_id": "qc-v1",
                "semantic_sha256": "a" * 64,
                "parameters": {},
                "artifacts": {
                    "artifact": {
                        "path": artifact.name,
                        "sha256": sha256_file(artifact),
                    }
                },
            }
        )
        + "\n"
    )
    return manifest


def test_stage_checkpoint_resolves_and_verifies_artifacts(tmp_path):
    checkpoint = _checkpoint(tmp_path)

    parsed = read_stage_checkpoint(checkpoint, expected_stage="alignment")

    assert parsed["stage"] == "alignment"
    assert parsed["artifacts"]["artifact"]["path"] == str(
        (tmp_path / "artifact.txt").resolve()
    )


def test_stage_checkpoint_rejects_wrong_boundary(tmp_path):
    with pytest.raises(AcquisitionError, match="expected 'qc'"):
        read_stage_checkpoint(_checkpoint(tmp_path), expected_stage="qc")


def test_stage_checkpoint_rejects_changed_artifact(tmp_path):
    checkpoint = _checkpoint(tmp_path)
    (tmp_path / "artifact.txt").write_text("changed\n")

    with pytest.raises(AcquisitionError, match="SHA-256 mismatch"):
        read_stage_checkpoint(checkpoint)
