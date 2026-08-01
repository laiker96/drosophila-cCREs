import hashlib
import shutil
import threading
from functools import partial
from http.server import BaseHTTPRequestHandler, SimpleHTTPRequestHandler, ThreadingHTTPServer

import pytest

from short_read_processing.accessions import FilePlan, RunPlan
from short_read_processing.downloader import (
    DownloadOptions,
    _verified_ena_file,
    download_ena,
    download_plans,
)


class _QuietHandler(SimpleHTTPRequestHandler):
    def log_message(self, format, *args):
        return


class _ForbiddenHandler(BaseHTTPRequestHandler):
    def do_GET(self):
        self.send_error(403)

    def do_HEAD(self):
        self.send_error(403)

    def log_message(self, format, *args):
        return


@pytest.mark.skipif(shutil.which("aria2c") is None, reason="aria2c is not installed")
def test_real_aria2_download_and_checksum(tmp_path):
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    source = source_dir / "SRR123.fastq.gz"
    payload = b"deterministic-fastq-block\n" * 100_000
    source.write_bytes(payload)
    checksum = hashlib.md5(payload).hexdigest()

    handler = partial(_QuietHandler, directory=str(source_dir))
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        destination = tmp_path / "output" / "SRR123" / source.name
        plan = RunPlan(
            requested_accession="SRR123",
            experiment_accession="SRX123",
            run_accession="SRR123",
            library_layout="SINGLE",
            backend="ena",
            run_dir=destination.parent,
            files=[
                FilePlan(
                    url=f"http://127.0.0.1:{server.server_port}/{source.name}",
                    md5=checksum,
                    size_bytes=len(payload),
                    path=destination,
                    mate="r1",
                )
            ],
        )
        download_ena(
            [plan],
            DownloadOptions(file_jobs=2, connections=4, sra_jobs=1, threads=2),
        )
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)

    assert destination.read_bytes() == payload
    assert plan.status == "downloaded"
    assert _verified_ena_file(plan.files[0])

    download_ena(
        [plan],
        DownloadOptions(file_jobs=2, connections=4, sra_jobs=1, threads=2),
    )

    assert plan.status == "existing"


@pytest.mark.skipif(shutil.which("aria2c") is None, reason="aria2c is not installed")
def test_real_aria2_http_failure_selects_expected_size_file_for_fallback(
    tmp_path, monkeypatch
):
    server = ThreadingHTTPServer(("127.0.0.1", 0), _ForbiddenHandler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        destination = tmp_path / "output" / "SRR403" / "SRR403.fastq.gz"
        destination.parent.mkdir(parents=True)
        expected = b"expected"
        destination.write_bytes(b"corrupt!")
        plan = RunPlan(
            requested_accession="SRR403",
            experiment_accession="SRX403",
            run_accession="SRR403",
            library_layout="SINGLE",
            backend="ena",
            run_dir=destination.parent,
            files=[
                FilePlan(
                    url=f"http://127.0.0.1:{server.server_port}/SRR403.fastq.gz",
                    md5=hashlib.md5(expected).hexdigest(),
                    size_bytes=len(expected),
                    path=destination,
                    mate="r1",
                )
            ],
        )
        received_sra = []

        def fake_download_sra(plans, options):
            received_sra.extend(plans)
            for item in plans:
                item.status = "downloaded"

        monkeypatch.setattr(
            "short_read_processing.downloader.download_sra", fake_download_sra
        )
        download_plans(
            [plan],
            DownloadOptions(
                file_jobs=1,
                connections=1,
                sra_jobs=1,
                threads=1,
                checksum_retries=0,
                ena_fallback=True,
            ),
        )
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)

    assert received_sra == [plan]
    assert plan.backend == "sra"
    assert plan.status == "downloaded"
