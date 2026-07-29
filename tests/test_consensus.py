import json
from pathlib import Path

import pytest

from short_read_processing.accessions import AcquisitionError
from short_read_processing.consensus import build_condition_consensus, condition_specs


def peak(start: int, end: int, name: str) -> str:
    return f"chr1\t{start}\t{end}\t{name}\t100\t.\n"


def test_conditions_require_every_library_once():
    values = [
        {
            "id": "eye",
            "label": "Eye disc",
            "samples": ["eye_rep1", "eye_rep2"],
        }
    ]
    conditions = condition_specs(
        values,
        sample_ids=["eye_rep1", "eye_rep2"],
        minimum_replicates=2,
    )
    assert conditions[0].samples == ("eye_rep1", "eye_rep2")

    with pytest.raises(AcquisitionError, match="do not assign"):
        condition_specs(
            values,
            sample_ids=["eye_rep1", "eye_rep2", "eye_rep3"],
            minimum_replicates=2,
        )


def test_condition_consensus_requires_fractional_coverage_in_two_replicates(
    tmp_path: Path,
):
    pooled = tmp_path / "pooled.bed"
    pooled.write_text(peak(0, 100, "pooled_1") + peak(200, 300, "pooled_2"))
    replicate_1 = tmp_path / "rep1.bed"
    replicate_1.write_text(peak(0, 100, "rep1_1") + peak(200, 240, "rep1_2"))
    replicate_2 = tmp_path / "rep2.bed"
    replicate_2.write_text(peak(0, 60, "rep2_1") + peak(200, 300, "rep2_2"))
    output = tmp_path / "consensus.bed"
    support = tmp_path / "support.tsv"
    stats = tmp_path / "stats.json"

    metrics = build_condition_consensus(
        condition_id="eye",
        pooled_peaks=pooled,
        replicate_peaks={"eye_rep1": replicate_1, "eye_rep2": replicate_2},
        output_bed=output,
        support_tsv=support,
        stats_json=stats,
    )

    rows = [line.split("\t") for line in output.read_text().splitlines()]
    assert [row[3] for row in rows] == ["pooled_1"]
    assert rows[0][6:12] == ["eye", "2", "2", "1", "eye_rep1,eye_rep2", "qpois"]
    assert support.read_text().splitlines()[2].split("\t")[2:4] == ["0", "1"]
    assert metrics["retained_replicate_supported_peaks"] == 1
    assert json.loads(stats.read_text()) == metrics
