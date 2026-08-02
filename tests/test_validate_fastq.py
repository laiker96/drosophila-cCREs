import gzip
import subprocess
import sys
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "workflow" / "scripts" / "validate_fastq.py"


def _write_fastq(path: Path, records: list[tuple[bytes, bytes, bytes]]) -> None:
    with gzip.open(path, "wb") as handle:
        for header, sequence, quality in records:
            handle.write(header + b"\n" + sequence + b"\n+\n" + quality + b"\n")


def _run(*paths: Path, paired: bool = False) -> subprocess.CompletedProcess[str]:
    command = [sys.executable, str(SCRIPT)]
    if paired:
        command.append("--paired")
    command.extend(str(path) for path in paths)
    return subprocess.run(command, text=True, capture_output=True, check=False)


def test_validate_fastq_accepts_valid_paired_files(tmp_path):
    r1 = tmp_path / "R1.fastq.gz"
    r2 = tmp_path / "R2.fastq.gz"
    records = [(b"@read1", b"ACGTN", b"IIIII")]
    _write_fastq(r1, records)
    _write_fastq(r2, records)

    result = _run(r1, r2, paired=True)

    assert result.returncode == 0
    assert result.stdout.count("records=1") == 2


def test_validate_fastq_rejects_invalid_sequence_byte(tmp_path):
    path = tmp_path / "bad-sequence.fastq.gz"
    _write_fastq(path, [(b"@read1", b"AC@T", b"IIII")])

    result = _run(path)

    assert result.returncode != 0
    assert "invalid sequence byte(s): 64" in result.stderr


def test_validate_fastq_rejects_invalid_quality_byte(tmp_path):
    path = tmp_path / "bad-quality.fastq.gz"
    _write_fastq(path, [(b"@read1", b"ACGT", b"II\xffI")])

    result = _run(path)

    assert result.returncode != 0
    assert "invalid quality byte(s): 255" in result.stderr


def test_validate_fastq_rejects_unequal_paired_record_counts(tmp_path):
    r1 = tmp_path / "R1.fastq.gz"
    r2 = tmp_path / "R2.fastq.gz"
    record = (b"@read1", b"ACGT", b"IIII")
    _write_fastq(r1, [record, record])
    _write_fastq(r2, [record])

    result = _run(r1, r2, paired=True)

    assert result.returncode != 0
    assert "paired FASTQs contain different record counts" in result.stderr
