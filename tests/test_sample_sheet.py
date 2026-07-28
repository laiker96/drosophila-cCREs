from collections import defaultdict
from pathlib import Path

import pytest

from short_read_processing.accessions import AcquisitionError
from short_read_processing.sample_sheet import read_sample_sheet


def test_four_column_csv_is_accepted_and_defaults_are_resolved(tmp_path: Path):
    sheet = tmp_path / "samples.csv"
    sheet.write_text(
        "accession,library_id,assay,context\n"
        "SRR123,atac_rep1,atac,eye\n"
    )
    row = read_sample_sheet(sheet)[0]

    assert row["accession"] == "SRR123"
    assert row["library_id"] == "atac_rep1"
    assert row["role"] == "treatment"
    assert row["peak_caller"] == "callpeak"
    assert row["macs3_qvalue"] == 0.1
    assert row["macs3_shift"] == -75
    assert row["macs3_extsize"] == 150
    assert row["adapter_preset"] == "nextera"
    assert row["mapq_minimum"] == 30


def test_h3k27ac_alias_resolves_histone_defaults(tmp_path: Path):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "accession\tlibrary_id\tassay\tcontext\n"
        "SRR123\th3_rep1\th3k27ac\teye\n"
    )
    row = read_sample_sheet(sheet)[0]

    assert row["assay"] == "chip_histone"
    assert row["macs3_broad"] is True
    assert row["macs3_broad_cutoff"] == 0.1
    assert row["adapter_preset"] == "truseq"


def test_missing_required_column_fails(tmp_path: Path):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "accession\tlibrary_id\tassay\n"
        "SRR123\tatac_rep1\tatac\n"
    )
    with pytest.raises(AcquisitionError, match="missing required columns: context"):
        read_sample_sheet(sheet)


def test_technical_runs_must_have_identical_contexts(tmp_path: Path):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "accession\tlibrary_id\tassay\tcontext\n"
        "SRR123\tatac_rep1\tatac\teye\n"
        "SRR124\tatac_rep1\tatac\twing\n"
    )
    with pytest.raises(AcquisitionError, match="disagree on: context"):
        read_sample_sheet(sheet)


def test_matched_input_must_name_control_in_same_context(tmp_path: Path):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "accession\tlibrary_id\tassay\tcontext\trole\tcontrol_library\n"
        "SRR123\th3_rep1\th3k27ac\teye\ttreatment\tinput_rep1\n"
        "SRR124\tinput_rep1\th3k27ac\twing\tcontrol\t\n"
    )
    with pytest.raises(AcquisitionError, match="same assay and context"):
        read_sample_sheet(sheet)


def test_chip_callpeak_shift_requires_bam_format(tmp_path: Path):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "accession\tlibrary_id\tassay\tcontext\tmacs3_shift\tmacs3_extsize\n"
        "SRR123\ttf_rep1\tchip_tf\teye\t-75\t150\n"
    )
    with pytest.raises(AcquisitionError, match="requires macs3_format=BAM"):
        read_sample_sheet(sheet)


def test_curated_ip_only_table_is_minimal_and_covers_atlas_libraries():
    resources = Path(__file__).parents[1] / "resources"
    rows = read_sample_sheet(resources / "atlas_samples_ip_only.tsv")

    assert len(rows) == 38
    assert len({row["accession"] for row in rows}) == 38
    assert len({row["library_id"] for row in rows}) == 37
    assert [row["library_id"] for row in rows].count("e5_atac_rep1") == 2
    assert {row["context"] for row in rows if row["assay"] == "atac"} == {
        "ab", "e11", "e13", "e5", "ead", "hid", "lb", "o", "wid"
    }
    atac_libraries: dict[str, set[str]] = defaultdict(set)
    for row in rows:
        if row["assay"] == "atac":
            atac_libraries[str(row["context"])].add(str(row["library_id"]))
    assert all(len(libraries) >= 2 for libraries in atac_libraries.values())
    assert all(row["control_library"] is None for row in rows)


def test_curated_controlled_table_links_three_matched_inputs():
    resources = Path(__file__).parents[1] / "resources"
    rows = read_sample_sheet(resources / "atlas_samples_with_inputs.tsv")
    controls = {row["library_id"] for row in rows if row["role"] == "control"}
    treatments = {
        row["library_id"]: row["control_library"]
        for row in rows
        if row["control_library"]
    }

    assert len(rows) == 41
    assert controls == {"ab_input_rep1", "e5_input_rep1", "e5_input_rep2"}
    assert treatments == {
        "ab_h3k27ac_rep1": "ab_input_rep1",
        "e5_h3k27ac_rep1": "e5_input_rep1",
        "e5_h3k27ac_rep2": "e5_input_rep2",
    }
