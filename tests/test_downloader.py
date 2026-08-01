from pathlib import Path
import hashlib
import os
import subprocess

import pytest

from short_read_processing.accessions import AcquisitionError, FilePlan, RunPlan
from short_read_processing.downloader import (
    DownloadOptions,
    _discard_untracked_partial_files,
    _download_one_sra,
    _record_ena_verification,
    _run_aria2_with_checksum_retries,
    _verified_ena_file,
    download_ena,
    download_plans,
)


def _file_plan(path: Path, size: int) -> FilePlan:
    return FilePlan(
        url=f"https://example.org/{path.name}",
        md5="",
        size_bytes=size,
        path=path,
    )


def _checksummed_file_plan(path: Path, content: bytes) -> FilePlan:
    return FilePlan(
        url=f"https://example.org/{path.name}",
        md5=hashlib.md5(content).hexdigest(),
        size_bytes=len(content),
        path=path,
    )


def test_checksum_record_reuses_only_unchanged_fastq(tmp_path):
    item = _checksummed_file_plan(tmp_path / "reads.fastq.gz", b"correct")
    item.path.write_bytes(b"correct")

    _record_ena_verification(item)

    assert _verified_ena_file(item)
    stat = item.path.stat()
    os.utime(
        item.path,
        ns=(stat.st_atime_ns, stat.st_mtime_ns + 1),
    )
    assert not _verified_ena_file(item)


def test_download_ena_queues_only_unverified_fastqs(tmp_path, monkeypatch):
    existing_item = _checksummed_file_plan(
        tmp_path / "existing" / "existing.fastq.gz", b"existing"
    )
    pending_item = _checksummed_file_plan(
        tmp_path / "pending" / "pending.fastq.gz", b"pending"
    )
    existing_item.path.parent.mkdir(parents=True)
    pending_item.path.parent.mkdir(parents=True)
    existing_item.path.write_bytes(b"existing")
    pending_item.path.write_bytes(b"pending")
    _record_ena_verification(existing_item)
    plans = [
        RunPlan(
            requested_accession="SRR111111",
            experiment_accession="SRX999999",
            run_accession="SRR111111",
            library_layout="SINGLE",
            backend="ena",
            run_dir=existing_item.path.parent,
            files=[existing_item],
        ),
        RunPlan(
            requested_accession="SRR222222",
            experiment_accession="SRX999999",
            run_accession="SRR222222",
            library_layout="SINGLE",
            backend="ena",
            run_dir=pending_item.path.parent,
            files=[pending_item],
        ),
    ]

    def fake_run(command):
        input_option = next(
            item for item in command if item.startswith("--input-file=")
        )
        input_text = Path(input_option.split("=", 1)[1]).read_text()
        assert pending_item.path.name in input_text
        assert existing_item.path.name not in input_text
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr(
        "short_read_processing.downloader._require_executable", lambda name: name
    )
    monkeypatch.setattr("short_read_processing.downloader.subprocess.run", fake_run)

    download_ena(
        plans,
        DownloadOptions(file_jobs=1, connections=1, sra_jobs=1, threads=1),
    )

    assert [plan.status for plan in plans] == ["existing", "downloaded"]
    assert _verified_ena_file(existing_item)
    assert _verified_ena_file(pending_item)


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


def test_unexplained_aria_failure_does_not_record_fastq(tmp_path, monkeypatch):
    item = _checksummed_file_plan(tmp_path / "complete.fastq.gz", b"complete")
    item.path.write_bytes(b"complete")
    input_path = tmp_path / "aria2-input.txt"
    input_path.write_text("initial\n")
    session_path = tmp_path / "aria2.session"

    monkeypatch.setattr(
        "short_read_processing.downloader.subprocess.run",
        lambda command: subprocess.CompletedProcess(command, 22),
    )

    with pytest.raises(AcquisitionError, match="exit code 22"):
        _run_aria2_with_checksum_retries(
            ["aria2c", f"--input-file={input_path}"],
            input_path=input_path,
            files=[item],
            checksum_retries=0,
            session_path=session_path,
        )

    assert not _verified_ena_file(item)


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


def test_success_exit_trusts_aria2_checksum_verification(tmp_path, monkeypatch):
    verified = FilePlan(
        url="https://example.org/verified.fastq.gz",
        md5=hashlib.md5(b"correct").hexdigest(),
        size_bytes=7,
        path=tmp_path / "verified.fastq.gz",
    )
    verified.path.write_bytes(b"correct")
    input_path = tmp_path / "aria2-input.txt"
    input_path.write_text("initial\n")
    attempts = []

    def fake_run(command):
        attempts.append(command)
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr("short_read_processing.downloader.subprocess.run", fake_run)

    def unexpected_rehash(*args, **kwargs):
        raise AssertionError("successful aria2 downloads must not be rehashed")

    monkeypatch.setattr("short_read_processing.downloader.hashlib.md5", unexpected_rehash)

    _run_aria2_with_checksum_retries(
        ["aria2c", f"--input-file={input_path}"],
        input_path=input_path,
        files=[verified],
        checksum_retries=1,
    )

    assert len(attempts) == 1
    assert verified.path.read_bytes() == b"correct"


def _run_plan(tmp_path, run, *, backend="ena"):
    run_dir = tmp_path / run
    files = [
        _file_plan(run_dir / f"{run}_1.fastq.gz", 2),
        _file_plan(run_dir / f"{run}_2.fastq.gz", 2),
    ]
    return RunPlan(
        requested_accession=run,
        experiment_accession="SRX999999",
        run_accession=run,
        library_layout="PAIRED",
        backend=backend,
        run_dir=run_dir,
        files=files if backend == "ena" else [],
    )


def test_auto_backend_falls_back_only_failed_ena_runs(tmp_path, monkeypatch):
    complete = _run_plan(tmp_path, "SRR111111")
    failed = _run_plan(tmp_path, "SRR222222")
    direct_sra = _run_plan(tmp_path, "SRR333333", backend="sra")
    received_sra = []
    aria2_commands = []

    for plan in (complete, failed):
        for item in plan.files:
            item.path.parent.mkdir(parents=True, exist_ok=True)
            item.path.write_bytes(b"ok")
            item.md5 = hashlib.md5(b"ok").hexdigest()

    failed_file = failed.files[0]
    failed_url = failed.files[0].url

    def fake_run(command):
        aria2_commands.append(command)
        session_option = next(
            item for item in command if item.startswith("--save-session=")
        )
        session_path = Path(session_option.split("=", 1)[1])
        session_path.write_text(
            f"{failed_url}\n"
            f"  dir={failed.run_dir}\n"
            f"  out={failed.files[0].path.name}\n",
            encoding="utf-8",
        )
        return subprocess.CompletedProcess(command, 22)

    def fake_download_sra(plans, options):
        received_sra.extend(plans)
        for plan in plans:
            plan.status = "downloaded"

    monkeypatch.setattr(
        "short_read_processing.downloader._require_executable", lambda name: name
    )
    monkeypatch.setattr("short_read_processing.downloader.subprocess.run", fake_run)
    monkeypatch.setattr("short_read_processing.downloader.download_sra", fake_download_sra)

    download_plans(
        [complete, failed, direct_sra],
        DownloadOptions(
            file_jobs=1,
            connections=1,
            sra_jobs=1,
            threads=1,
            ena_fallback=True,
        ),
    )

    assert complete.backend == "ena"
    assert complete.status == "downloaded"
    assert failed.backend == "sra"
    assert failed.files == []
    assert received_sra == [direct_sra, failed]
    assert len(aria2_commands) == 1
    assert any(item.startswith("--save-session=") for item in aria2_commands[0])
    assert all(_verified_ena_file(item) for item in complete.files)
    assert not _verified_ena_file(failed_file)


def test_explicit_ena_backend_does_not_fall_back(tmp_path, monkeypatch):
    failed = _run_plan(tmp_path, "SRR222222")

    def fake_download_ena(plans, options):
        raise AcquisitionError("aria2c failed with exit code 22")

    def unexpected_download_sra(plans, options):
        raise AssertionError("explicit ENA mode must not use SRA Toolkit")

    monkeypatch.setattr("short_read_processing.downloader.download_ena", fake_download_ena)
    monkeypatch.setattr(
        "short_read_processing.downloader.download_sra", unexpected_download_sra
    )

    with pytest.raises(AcquisitionError, match="exit code 22"):
        download_plans(
            [failed],
            DownloadOptions(
                file_jobs=1,
                connections=1,
                sra_jobs=1,
                threads=1,
            ),
        )


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
        if label.startswith("prefetch"):
            output = Path(command[command.index("--output-directory") + 1])
            (output / "SRR123456").mkdir(parents=True)
        if label.startswith("fasterq-dump"):
            expected = (
                tmp_path / "raw" / ".sra-cache" / "SRR123456" / "SRR123456"
            )
            assert Path(command[1]) == expected
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
