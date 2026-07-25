import csv
import gzip
from pathlib import Path

import pytest

from short_read_processing.activity import (
    build_activity_outputs,
    count_units,
    cpm_per_kb,
    read_master_elements,
    sha256_file,
    tie_aware_quantile_normalize,
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


def test_tie_aware_quantile_normalization():
    normalized, tie_groups = tie_aware_quantile_normalize(
        [0.0, 0.0, 2.0],
        [1.0, 3.0, 9.0],
    )
    assert normalized == [2.0, 2.0, 9.0]
    assert tie_groups == 1


def test_activity_outputs_are_deterministic_and_keep_raw_values(tmp_path):
    master, summits = _master(tmp_path)
    elements = read_master_elements(master, summits)
    specifications = [
        ("atlas_atac_1", "atac", "atlas", "ctx", [1, 3], 100),
        ("atlas_atac_2", "atac", "atlas", "ctx", [4, 2], 1000),
        ("atlas_h3_1", "h3k27ac", "atlas", "ctx", [2, 5], 200),
        ("atlas_h3_2", "h3k27ac", "atlas", "ctx", [3, 1], 400),
        ("ref_atac_1", "atac", "reference", "ref", [2, 6], 100),
        ("ref_atac_2", "atac", "reference", "ref", [4, 2], 200),
        ("ref_h3_1", "h3k27ac", "reference", "ref", [1, 8], 100),
        ("ref_h3_2", "h3k27ac", "reference", "ref", [2, 4], 200),
    ]
    signal_paths = {}
    for library_id, assay, cohort, context, counts, total in specifications:
        signal = tmp_path / "signals" / f"{library_id}.tsv.gz"
        summary = tmp_path / "signals" / f"{library_id}.json"
        write_library_signal(
            elements=elements,
            counts=counts,
            total_units=total,
            library_id=library_id,
            assay=assay,
            cohort=cohort,
            context=context,
            output=signal,
            summary=summary,
        )
        signal_paths[library_id] = signal

    outputs = {
        "output_library_signal": tmp_path / "library.tsv.gz",
        "output_context_signal": tmp_path / "context.tsv.gz",
        "output_reference": tmp_path / "reference.tsv.gz",
        "output_activity": tmp_path / "activity.tsv.gz",
        "output_context_views": {"ctx": tmp_path / "ctx.tsv.gz"},
        "output_metrics": tmp_path / "metrics.json",
        "output_provenance": tmp_path / "provenance.json",
    }
    build_activity_outputs(
        signal_paths=signal_paths,
        atlas_contexts=["ctx"],
        reference_context="ref",
        provenance={"test": True},
        **outputs,
    )
    hashes = {
        name: sha256_file(path)
        for name, path in outputs.items()
        if isinstance(path, Path)
    }

    build_activity_outputs(
        signal_paths=dict(reversed(list(signal_paths.items()))),
        atlas_contexts=["ctx"],
        reference_context="ref",
        provenance={"test": True},
        **outputs,
    )

    assert {
        name: sha256_file(path)
        for name, path in outputs.items()
        if isinstance(path, Path)
    } == hashes
    with gzip.open(outputs["output_activity"], "rt", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [int(row["width_bp"]) for row in rows] == [31, 69]
    assert [int(row["atac_raw_count_sum"]) for row in rows] == [5, 5]
    assert all(float(row["activity"]) >= 0 for row in rows)
