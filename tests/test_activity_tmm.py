import csv
import gzip
import math
from pathlib import Path
import subprocess

import pytest

from short_read_processing.activity import read_master_elements, write_library_signal
from short_read_processing.activity_tmm import (
    BACKGROUND_COUNT_FIELDS,
    TMM_ACTIVITY_FIELDS,
    TMM_BACKGROUND_METHOD,
    TMM_FACTOR_FIELDS,
    TMM_MASTER_METHOD,
    _write_gzip_dict_rows,
    build_tmm_activity_outputs,
    build_tmm_inputs,
    write_background_bins,
)


def test_streaming_gzip_has_no_temporary_filename_and_is_deterministic(tmp_path):
    first = tmp_path / "first.tsv.gz"
    second = tmp_path / "second.tsv.gz"
    rows = [{"element": "DHS0000001", "value": 1}]

    _write_gzip_dict_rows(first, ["element", "value"], rows)
    _write_gzip_dict_rows(second, ["element", "value"], rows)

    first_bytes = first.read_bytes()
    assert first_bytes == second.read_bytes()
    assert first_bytes[:3] == b"\x1f\x8b\x08"
    assert first_bytes[3] & 0x08 == 0  # gzip FNAME flag is absent


def _signals(tmp_path: Path) -> dict[str, Path]:
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
    elements = read_master_elements(master, summits)
    specifications = [
        ("atac_1", "atac", [1, 3], 100),
        ("atac_2", "atac", [4, 2], 1000),
        ("h3_1", "h3k27ac", [2, 5], 200),
        ("h3_2", "h3k27ac", [3, 1], 400),
    ]
    paths = {}
    for library_id, assay, counts, total in specifications:
        path = tmp_path / "signals" / f"{library_id}.tsv.gz"
        write_library_signal(
            elements=elements,
            counts=counts,
            total_units=total,
            library_id=library_id,
            assay=assay,
            cohort="atlas",
            context="ctx",
            output=path,
            summary=tmp_path / "signals" / f"{library_id}.json",
        )
        paths[library_id] = path
    return paths


def test_background_bins_use_only_configured_autosomes(tmp_path):
    chrom_sizes = tmp_path / "chrom.sizes"
    chrom_sizes.write_text("chr1\t25\nchrX\t20\nchr2\t10\n")
    autosomes = tmp_path / "autosomes.txt"
    autosomes.write_text("chr2\nchr1\n")
    output = tmp_path / "bins.bed"

    count = write_background_bins(
        chrom_sizes_path=chrom_sizes,
        autosomes_path=autosomes,
        output_path=output,
        bin_width=10,
    )

    assert count == 4
    assert output.read_text().splitlines() == [
        "chr1\t0\t10\tBGBIN0000001",
        "chr1\t10\t20\tBGBIN0000002",
        "chr1\t20\t25\tBGBIN0000003",
        "chr2\t0\t10\tBGBIN0000004",
    ]


def test_background_counter_does_not_require_sorted_fragment_midpoints(tmp_path):
    bins = tmp_path / "bins.bed"
    bins.write_text(
        "chr1\t0\t10000\tBGBIN0000001\n"
        "chr1\t10000\t20000\tBGBIN0000002\n"
    )
    units = tmp_path / "units.bed.gz"
    with gzip.open(units, "wt") as handle:
        handle.write(
            "chr1\t0\t18000\n"  # midpoint 9000
            "chr1\t5000\t7000\n"  # midpoint 6000: decreases despite sorted starts
            "chr1\t8000\t22000\n"  # midpoint 15000
            "chrX\t0\t100\n"
        )
    output = tmp_path / "counts.tsv.gz"
    subprocess.run(
        [
            "bash",
            str(
                Path(__file__).resolve().parents[1]
                / "src"
                / "count_activity_background.sh"
            ),
            str(units),
            str(bins),
            str(output),
            "1",
        ],
        check=True,
    )
    with gzip.open(output, "rt", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [int(row["raw_count"]) for row in rows] == [2, 1]


def test_tmm_inputs_preserve_raw_master_and_background_counts(tmp_path):
    signals = _signals(tmp_path)
    master_counts = tmp_path / "master-counts.tsv.gz"
    master_metadata = tmp_path / "master-metadata.tsv"
    metrics = build_tmm_inputs(
        method=TMM_MASTER_METHOD,
        signal_paths=signals,
        background_count_paths=None,
        output_counts=master_counts,
        output_metadata=master_metadata,
    )
    assert metrics == {
        "method": TMM_MASTER_METHOD,
        "feature_count": 2,
        "library_count": 4,
    }
    with gzip.open(master_counts, "rt", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows[0]["feature_id"] == "DHS0000001"
    assert [int(rows[0][item]) for item in sorted(signals)] == [1, 4, 2, 3]

    background_paths = {}
    for index, library_id in enumerate(sorted(signals), start=1):
        path = tmp_path / "background" / f"{library_id}.tsv.gz"
        path.parent.mkdir(exist_ok=True)
        with gzip.open(path, "wt", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=BACKGROUND_COUNT_FIELDS,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerow(
                {
                    "background_bin_id": "BGBIN0000001",
                    "chrom": "chr1",
                    "start": 0,
                    "end": 10,
                    "raw_count": index,
                }
            )
        background_paths[library_id] = path
    background_counts = tmp_path / "background-counts.tsv.gz"
    build_tmm_inputs(
        method=TMM_BACKGROUND_METHOD,
        signal_paths=signals,
        background_count_paths=background_paths,
        output_counts=background_counts,
        output_metadata=tmp_path / "background-metadata.tsv",
    )
    with gzip.open(background_counts, "rt", newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["feature_id"] == "BGBIN0000001"
    assert [int(row[item]) for item in sorted(signals)] == [1, 2, 3, 4]


def test_tmm_activity_uses_effective_library_sizes_and_keeps_plain_cpm(tmp_path):
    signals = _signals(tmp_path)
    factors = tmp_path / "factors.tsv"
    factor_values = {"atac_1": 2.0, "atac_2": 0.5, "h3_1": 2.0, "h3_2": 0.5}
    with factors.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=TMM_FACTOR_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for library_id in sorted(signals):
            assay = "atac" if library_id.startswith("atac") else "h3k27ac"
            total = {"atac_1": 100, "atac_2": 1000, "h3_1": 200, "h3_2": 400}[
                library_id
            ]
            factor = factor_values[library_id]
            writer.writerow(
                {
                    "library_id": library_id,
                    "assay": assay,
                    "context": "ctx",
                    "total_units": total,
                    "feature_count": 2,
                    "tmm_normalization_factor": factor,
                    "effective_library_size": total * factor,
                    "normalization_method": TMM_MASTER_METHOD,
                }
            )

    activity = tmp_path / "activity.tsv.gz"
    metrics = build_tmm_activity_outputs(
        method=TMM_MASTER_METHOD,
        signal_paths=signals,
        factor_path=factors,
        contexts=["ctx"],
        output_context_signal=tmp_path / "context.tsv.gz",
        output_activity=activity,
        output_metrics=tmp_path / "metrics.json",
        output_provenance=tmp_path / "provenance.json",
        provenance={"test": True},
    )
    assert metrics["factor_product_by_assay"] == {"atac": 1.0, "h3k27ac": 1.0}
    with gzip.open(activity, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assert reader.fieldnames == TMM_ACTIVITY_FIELDS
        first = next(reader)
    expected_atac = ((1e9 / (200 * 31)) + (4e9 / (500 * 31))) / 2
    expected_h3 = ((2e9 / (400 * 31)) + (3e9 / (200 * 31))) / 2
    assert float(first["atac_normalized_cpm_per_kb"]) == pytest.approx(expected_atac)
    assert float(first["h3k27ac_normalized_cpm_per_kb"]) == pytest.approx(expected_h3)
    assert float(first["activity"]) == pytest.approx(math.sqrt(expected_atac * expected_h3))
    assert float(first["atac_cpm_per_kb"]) != pytest.approx(expected_atac)
