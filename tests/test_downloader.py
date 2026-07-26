from pathlib import Path
import hashlib
import subprocess

import pytest

from short_read_processing.accessions import AcquisitionError, FilePlan, RunPlan
from short_read_processing.downloader import (
    _discard_untracked_partial_files,
    _download_one_sra,
    _run_aria2_with_checksum_retries,
)


def _file_plan(path: Path, size: int) -> FilePlan:
    return FilePlan(url="https://example.org/reads.fastq.gz", md5="", size_bytes=size, path=path)


def test_discards_size_mismatched_file_without_aria2_state(tmp_path):
    partial = tmp_path / "reads.fastq.gz"
    partial.write_bytes(b"partial")

    discarded = _discard_untracked_partial_files([_file_plan(partial, 100)])

    assert discarded == [partial]
    assert not partial.exists()


def test_keeps_aria2_managed_partial_file(tmp_path):
    partial = tmp_path / "reads.fastq.gz"
    partial.write_bytes(b"partial")
    partial.with_name(partial.name + ".aria2").touch()

    discarded = _discard_untracked_partial_files([_file_plan(partial, 100)])

    assert discarded == []
    assert partial.exists()


def test_checksum_failure_retries_only_incomplete_files(tmp_path, monkeypatch):
    complete = _file_plan(tmp_path / "complete.fastq.gz", 8)
    failed = _file_plan(tmp_path / "failed.fastq.gz", 7)
    complete.path.write_bytes(b"complete")
    failed.path.write_bytes(b"corrupt")
    control = failed.path.with_name(failed.path.name + ".aria2")
    control.touch()
    input_path = tmp_path / "aria2-input.txt"
    input_path.write_text("initial\n")
    attempts = []

    def fake_run(command):
        attempts.append(command)
        if len(attempts) == 1:
            return subprocess.CompletedProcess(command, 32)
        failed.path.write_bytes(b"correct")
        control.unlink()
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr("short_read_processing.downloader.subprocess.run", fake_run)

    _run_aria2_with_checksum_retries(
        ["aria2c", f"--input-file={input_path}"],
        input_path=input_path,
        files=[complete, failed],
        checksum_retries=3,
    )

    assert len(attempts) == 2
    assert failed.path.name in input_path.read_text()
    assert complete.path.name not in input_path.read_text()


def test_checksum_failure_stops_after_configured_retries(tmp_path, monkeypatch):
    failed = _file_plan(tmp_path / "failed.fastq.gz", 7)
    failed.path.write_bytes(b"corrupt")
    failed.path.with_name(failed.path.name + ".aria2").touch()
    input_path = tmp_path / "aria2-input.txt"
    input_path.write_text("initial\n")
    attempts = []

    def fake_run(command):
        attempts.append(command)
        return subprocess.CompletedProcess(command, 32)

    monkeypatch.setattr("short_read_processing.downloader.subprocess.run", fake_run)

    with pytest.raises(AcquisitionError, match="after 3 attempt"):
        _run_aria2_with_checksum_retries(
            ["aria2c", f"--input-file={input_path}"],
            input_path=input_path,
            files=[failed],
            checksum_retries=2,
        )

    assert len(attempts) == 3


def test_checksum_retry_restarts_complete_size_corrupt_file(tmp_path, monkeypatch):
    failed = FilePlan(
        url="https://example.org/failed.fastq.gz",
        md5=hashlib.md5(b"correct").hexdigest(),
        size_bytes=7,
        path=tmp_path / "failed.fastq.gz",
    )
    failed.path.write_bytes(b"corrupt")
    input_path = tmp_path / "aria2-input.txt"
    input_path.write_text("initial\n")
    attempts = []

    def fake_run(command):
        attempts.append(command)
        if len(attempts) == 1:
            return subprocess.CompletedProcess(command, 32)
        assert not failed.path.exists()
        failed.path.write_bytes(b"correct")
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr("short_read_processing.downloader.subprocess.run", fake_run)

    _run_aria2_with_checksum_retries(
        ["aria2c", f"--input-file={input_path}"],
        input_path=input_path,
        files=[failed],
        checksum_retries=1,
    )

    assert len(attempts) == 2
    assert failed.path.read_bytes() == b"correct"


def test_success_exit_still_rechecks_and_retries_corrupt_file(tmp_path, monkeypatch):
    failed = FilePlan(
        url="https://example.org/failed.fastq.gz",
        md5=hashlib.md5(b"correct").hexdigest(),
        size_bytes=7,
        path=tmp_path / "failed.fastq.gz",
    )
    failed.path.write_bytes(b"corrupt")
    input_path = tmp_path / "aria2-input.txt"
    input_path.write_text("initial\n")
    attempts = []

    def fake_run(command):
        attempts.append(command)
        if len(attempts) == 1:
            return subprocess.CompletedProcess(command, 0)
        assert not failed.path.exists()
        failed.path.write_bytes(b"correct")
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr("short_read_processing.downloader.subprocess.run", fake_run)

    _run_aria2_with_checksum_retries(
        ["aria2c", f"--input-file={input_path}"],
        input_path=input_path,
        files=[failed],
        checksum_retries=1,
    )

    assert len(attempts) == 2
    assert failed.path.read_bytes() == b"correct"


def test_sra_download_stages_fastqs_before_completion(tmp_path, monkeypatch):
    run_dir = tmp_path / "raw" / "SRR123456"
    run_dir.mkdir(parents=True)
    (run_dir / "SRR123456_1.fastq.gz").write_bytes(b"partial")
    plan = RunPlan(
        requested_accession="SRR123456",
        experiment_accession="SRX999999",
        run_accession="SRR123456",
        library_layout="PAIRED",
        backend="sra",
        run_dir=run_dir,
    )

    monkeypatch.setattr(
        "short_read_processing.downloader._require_executable", lambda name: name
    )

    def fake_run(command, *, label):
        if label.startswith("fasterq-dump"):
            output = Path(command[command.index("--outdir") + 1])
            (output / "SRR123456_1.fastq").write_bytes(b"r1")
            (output / "SRR123456_2.fastq").write_bytes(b"r2")

    def fake_gzip(path, *, threads):
        compressed = path.with_suffix(path.suffix + ".gz")
        compressed.write_bytes(path.read_bytes())
        path.unlink()
        return compressed

    monkeypatch.setattr("short_read_processing.downloader._run", fake_run)
    monkeypatch.setattr("short_read_processing.downloader._gzip_fastq", fake_gzip)

    _download_one_sra(plan, threads=2, keep_cache=False)

    assert (run_dir / ".download-complete").read_text() == "complete\n"
    assert (run_dir / "SRR123456_1.fastq.gz").read_bytes() == b"r1"
    assert (run_dir / "SRR123456_2.fastq.gz").read_bytes() == b"r2"
    assert [item.mate for item in plan.files] == ["r1", "r2"]
