import gzip
from pathlib import Path

import pytest

from short_read_processing.activity import (
    count_units,
    cpm_per_kb,
    read_master_elements,
    sha256_file,
    write_library_signal,
)


def _master(tmp_path: Path):
    master = tmp_path / "master.bed"
    summits = tmp_path / "summits.bed"
    master.write_text(
        "chr1\t0\t31\tDHS0000001\t0\t.\n"
        "chr1\t31\t100\tDHS0000002\t0\t.\n"
    )
    summits.write_text(
        "chr1\t15\t16\tDHS0000001\t0\t.\n"
        "chr1\t60\t61\tDHS0000002\t0\t.\n"
    )
    return master, summits


def test_cpm_per_kb_uses_original_31bp_width():
    assert cpm_per_kb(31, 1_000_000, 31) == pytest.approx(1000.0)
    with pytest.raises(ValueError, match="width"):
        cpm_per_kb(1, 10, 0)
    with pytest.raises(ValueError, match="units"):
        cpm_per_kb(1, 0, 31)


def test_insertion_boundary_counting(tmp_path):
    master, summits = _master(tmp_path)
    units = tmp_path / "insertions.bed.gz"
    with gzip.open(units, "wt") as handle:
        handle.write(
            "chr1\t0\t1\n"
            "chr1\t30\t31\n"
            "chr1\t31\t32\n"
            "chr1\t99\t100\n"
            "chr1\t100\t101\n"
        )

    counts = count_units(
        read_master_elements(master, summits),
        units,
        expected_total=5,
    )

    assert counts == [2, 2]


def test_fragment_can_overlap_two_nonoverlapping_elements(tmp_path):
    master, summits = _master(tmp_path)
    units = tmp_path / "fragments.bed"
    units.write_text("chr1\t20\t40\nchr1\t100\t120\n")

    counts = count_units(
        read_master_elements(master, summits),
        units,
        expected_total=2,
    )

    assert counts == [1, 1]


def test_units_on_valid_chromosome_without_master_elements_stay_in_denominator(
    tmp_path,
):
    master, summits = _master(tmp_path)
    units = tmp_path / "fragments.bed"
    units.write_text("chr1\t20\t40\nchr2\t0\t10\n")

    counts = count_units(
        read_master_elements(master, summits),
        units,
        expected_total=2,
        chromosome_order=["chr1", "chr2"],
    )

    assert counts == [1, 1]


def test_library_signal_is_deterministic_and_keeps_variable_widths(tmp_path):
    master, summits = _master(tmp_path)
    elements = read_master_elements(master, summits)
    signal = tmp_path / "signals" / "atac.tsv.gz"
    summary = tmp_path / "signals" / "atac.json"
    arguments = {
        "elements": elements,
        "counts": [1, 3],
        "total_units": 100,
        "library_id": "atac_1",
        "assay": "atac",
        "cohort": "atlas",
        "context": "ctx",
        "output": signal,
        "summary": summary,
    }
    write_library_signal(**arguments)
    first_hashes = (sha256_file(signal), sha256_file(summary))
    write_library_signal(**arguments)
    assert (sha256_file(signal), sha256_file(summary)) == first_hashes
    with gzip.open(signal, "rt") as handle:
        widths = [int(line.split("\t")[5]) for line in list(handle)[1:]]
    assert widths == [31, 69]
