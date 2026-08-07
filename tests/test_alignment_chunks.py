import gzip
import importlib.util
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_script(name: str):
    path = REPO_ROOT / "workflow" / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


split_fastq_chunks = _load_script("split_fastq_chunks")
aggregate_bowtie2_logs = _load_script("aggregate_bowtie2_logs")


def _write_fastq(path: Path, mate: int, names: list[str]) -> bytes:
    content = b"".join(
        f"@{name}/{mate}\nACGT\n+\nIIII\n".encode() for name in names
    )
    with gzip.open(path, "wb") as handle:
        handle.write(content)
    return content


def test_split_fastqs_preserves_pairs_and_is_deterministic(tmp_path):
    r1 = tmp_path / "lane_R1.fastq.gz"
    r2 = tmp_path / "lane_R2.fastq.gz"
    expected_r1 = _write_fastq(r1, 1, ["read1", "read2", "read3"])
    expected_r2 = _write_fastq(r2, 2, ["read1", "read2", "read3"])

    first = tmp_path / "first"
    second = tmp_path / "second"
    manifest = split_fastq_chunks.split_fastqs(
        [r1, r2],
        first,
        unit="sample.L001",
        layout="paired",
        records_per_chunk=2,
    )
    split_fastq_chunks.split_fastqs(
        [r1, r2],
        second,
        unit="sample.L001",
        layout="paired",
        records_per_chunk=2,
    )

    assert manifest["total_records"] == 3
    assert [chunk["records"] for chunk in manifest["chunks"]] == [2, 1]
    for mate, expected in (("r1", expected_r1), ("r2", expected_r2)):
        names = [chunk[mate] for chunk in manifest["chunks"]]
        restored = b"".join(gzip.open(first / name, "rb").read() for name in names)
        assert restored == expected
        assert [(first / name).read_bytes() for name in names] == [
            (second / name).read_bytes() for name in names
        ]


def test_split_fastqs_rejects_mismatched_mates(tmp_path):
    r1 = tmp_path / "lane_R1.fastq.gz"
    r2 = tmp_path / "lane_R2.fastq.gz"
    _write_fastq(r1, 1, ["read1"])
    _write_fastq(r2, 2, ["different"])

    with pytest.raises(ValueError, match="identifiers differ"):
        split_fastq_chunks.split_fastqs(
            [r1, r2],
            tmp_path / "chunks",
            unit="sample.L001",
            layout="paired",
            records_per_chunk=2,
        )

    assert not (tmp_path / "chunks").exists()


def test_aggregate_bowtie2_chunk_logs(tmp_path):
    paths = []
    for index, counts in enumerate(((10, 2, 7, 1), (5, 1, 3, 1)), start=1):
        total, zero, once, multiple = counts
        path = tmp_path / f"chunk{index}.log"
        path.write_text(
            f"{total} reads; of these:\n"
            f"  {total} (100.00%) were paired; of these:\n"
            f"    {zero} (0.00%) aligned concordantly 0 times\n"
            f"    {once} (0.00%) aligned concordantly exactly 1 time\n"
            f"    {multiple} (0.00%) aligned concordantly >1 times\n"
            "0.00% overall alignment rate\n"
        )
        paths.append(path)

    combined = aggregate_bowtie2_logs.aggregate(paths, "paired")

    assert "15 reads; of these:" in combined
    assert "3 (20.00%) aligned concordantly 0 times" in combined
    assert "80.00% overall alignment rate" in combined
    assert "restartable chunks: 2" in combined
