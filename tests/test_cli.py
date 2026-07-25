import argparse
from pathlib import Path
import sys

import pytest

from short_read_processing.accessions import AcquisitionError
from short_read_processing.cli import read_accession_column
from run_pipeline import main as run_pipeline_main, rule_threads


def test_read_configurable_accession_column(tmp_path):
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text("run_id\tlabel\nSRR123\ta\nERX456\tb\nSRR123\tduplicate\n")
    assert read_accession_column(metadata, "run_id") == ["SRR123", "ERX456"]


def test_read_accession_column_from_csv(tmp_path):
    metadata = tmp_path / "metadata.csv"
    metadata.write_text("run_id,label\nSRR123,a\nERR456,b\n")
    assert read_accession_column(metadata, "run_id") == ["SRR123", "ERR456"]


def test_missing_accession_column_is_clear(tmp_path):
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text("wrong\nSRR123\n")
    with pytest.raises(AcquisitionError, match="available columns: wrong"):
        read_accession_column(metadata, "run_id")


@pytest.mark.parametrize("value", ["align_lane", "align_lane=0", "=4", "align_lane=x"])
def test_rule_thread_override_rejects_invalid_values(value):
    with pytest.raises(argparse.ArgumentTypeError, match="RULE=THREADS"):
        rule_threads(value)


def test_rule_thread_override_accepts_positive_count():
    assert rule_threads("align_lane=16") == "align_lane=16"


def test_final_bam_stage_requires_explicit_manifest(monkeypatch, tmp_path):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "accession\tlibrary_id\tassay\tcontext\n"
        "SRR100001\tatac_rep1\tatac\teye\n"
    )
    monkeypatch.setattr(
        sys,
        "argv",
        ["run_pipeline.py", str(sheet), "--from-stage", "final-bam"],
    )

    with pytest.raises(SystemExit) as error:
        run_pipeline_main()

    assert error.value.code == 2


def test_reuse_stage_rejects_download_fallback_flags(monkeypatch, tmp_path):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "accession\tlibrary_id\tassay\tcontext\n"
        "SRR100001\tatac_rep1\tatac\teye\n"
    )
    manifest = tmp_path / "final-bams.tsv"
    manifest.write_text("")
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "run_pipeline.py",
            str(sheet),
            "--from-stage",
            "final-bam",
            "--final-bam-manifest",
            str(manifest),
            "--skip-download",
        ],
    )

    with pytest.raises(SystemExit) as error:
        run_pipeline_main()

    assert error.value.code == 2


def test_activity_stage_requires_all_explicit_artifact_inputs(monkeypatch, tmp_path):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "accession\tlibrary_id\tassay\tcontext\n"
        "SRR100001\tatac_rep1\tatac\teye\n"
    )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "run_pipeline.py",
            str(sheet),
            "--from-stage",
            "activity",
            "--master-manifest",
            str(tmp_path / "master.tsv"),
            "--activity-atlas-bam-manifest",
            str(tmp_path / "atlas.tsv"),
        ],
    )

    with pytest.raises(SystemExit) as error:
        run_pipeline_main()

    assert error.value.code == 2


def test_accession_stage_rejects_activity_artifacts(monkeypatch, tmp_path):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "accession\tlibrary_id\tassay\tcontext\n"
        "SRR100001\tatac_rep1\tatac\teye\n"
    )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "run_pipeline.py",
            str(sheet),
            "--activity-reference-sheet",
            str(tmp_path / "reference.tsv"),
        ],
    )

    with pytest.raises(SystemExit) as error:
        run_pipeline_main()

    assert error.value.code == 2


def test_nonactivity_reuse_stage_rejects_activity_options(monkeypatch, tmp_path):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "accession\tlibrary_id\tassay\tcontext\n"
        "SRR100001\tatac_rep1\tatac\teye\n"
    )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "run_pipeline.py",
            str(sheet),
            "--from-stage",
            "master",
            "--master-manifest",
            str(tmp_path / "master.tsv"),
            "--activity-reference-context",
            "s2_t0",
        ],
    )

    with pytest.raises(SystemExit) as error:
        run_pipeline_main()

    assert error.value.code == 2
