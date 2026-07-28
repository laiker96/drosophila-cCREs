import xml.etree.ElementTree as ET
from pathlib import Path

from build_igv_session import build_catalog_session, build_session


def touch_all(paths):
    for path in paths:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()


def test_session_contains_final_qpois_consensus_and_chip_tracks(tmp_path: Path):
    atac = tmp_path / "atac"
    condition = atac / "conditions" / "e5"
    touch_all(
        [
            condition / "tracks/e5.MACS3-pileup.unscaled.bw",
            condition / "tracks/e5.qpois.bw",
            condition / "peaks/e5.candidates.narrowPeak",
            condition / "peaks/e5.qpois-refined.bed",
            condition / "peaks/e5.replicate-supported.bed",
        ]
    )
    chip = tmp_path / "chip"
    touch_all(
        [
            chip / "tracks/e5_h3k27ac_rep1.CPM.bw",
            chip / "peaks/e5_h3k27ac_rep1/e5_h3k27ac_rep1_peaks.broadPeak",
        ]
    )
    output = tmp_path / "session.xml"

    counts = build_session(atac, output, "dm6", "All", chip)

    root = ET.parse(output).getroot()
    resources = root.findall("./Resources/Resource")
    tracks = root.findall("./Panel/Track")
    assert counts == (1, 1, 7)
    assert len(resources) == len(tracks) == 7
    assert all((output.parent / item.attrib["path"]).is_file() for item in resources)
    assert [track.attrib["name"] for track in tracks][-2:] == [
        "e5_h3k27ac_rep1 | ChIP CPM signal",
        "e5_h3k27ac_rep1 | MACS3 peaks",
    ]


def test_final_session_selects_one_chip_replicate_per_context(tmp_path: Path):
    atac = tmp_path / "atac"
    condition = atac / "conditions" / "e5"
    touch_all(
        [
            condition / "tracks/e5.MACS3-pileup.unscaled.bw",
            condition / "tracks/e5.qpois.bw",
            condition / "peaks/e5.replicate-supported.bed",
        ]
    )
    chip = tmp_path / "chip"
    touch_all(
        [
            chip / "tracks/e5_h3k27ac_rep1.CPM.bw",
            chip / "peaks/e5_h3k27ac_rep1/e5_h3k27ac_rep1_peaks.broadPeak",
            chip / "tracks/e5_h3k27ac_rep2.CPM.bw",
            chip / "peaks/e5_h3k27ac_rep2/e5_h3k27ac_rep2_peaks.broadPeak",
        ]
    )
    output = tmp_path / "session.xml"

    counts = build_session(
        atac,
        output,
        "dm6",
        "All",
        chip,
        final_atac_only=True,
        chip_one_per_context=True,
    )

    tracks = ET.parse(output).getroot().findall("./Panel/Track")
    assert counts == (1, 1, 5)
    assert [track.attrib["name"] for track in tracks] == [
        "E5 | MACS3 insertion pileup",
        "E5 | qpois signal",
        "E5 | replicate-supported peaks",
        "e5_h3k27ac_rep1 | ChIP CPM signal",
        "e5_h3k27ac_rep1 | MACS3 peaks",
    ]


def test_session_auto_includes_master_dhs_track(tmp_path: Path):
    atac = tmp_path / "atac"
    condition = atac / "conditions" / "e5"
    touch_all(
        [
            atac / "master/master_dhs.bed",
            condition / "tracks/e5.MACS3-pileup.unscaled.bw",
            condition / "tracks/e5.qpois.bw",
            condition / "peaks/e5.replicate-supported.bed",
        ]
    )
    output = tmp_path / "session.xml"

    counts = build_session(atac, output, "dm6", "All", final_atac_only=True)

    tracks = ET.parse(output).getroot().findall("./Panel/Track")
    assert counts == (1, 0, 4)
    assert tracks[0].attrib["name"] == "Master DHS registry"
    assert tracks[0].attrib["color"] == "106,27,154"


def test_catalog_session_has_relative_mean_signal_and_element_tracks(tmp_path: Path):
    inputs = {
        "atac_bigwig": tmp_path / "catalog/tracks/ctx.atac.bw",
        "h3k27ac_bigwig": tmp_path / "catalog/tracks/ctx.h3k27ac.bw",
        "context_dhs_bed": tmp_path / "catalog/bed/ctx.dhs.bed",
        "master_dhs_bed": tmp_path / "catalog/bed/master_dhs.bed",
        "active_elements_bed": tmp_path / "catalog/bed/ctx.active_elements.bed",
    }
    touch_all(list(inputs.values()))
    output = tmp_path / "catalog/igv/ctx.xml"

    count = build_catalog_session(
        context="ctx",
        genome="dm6",
        output=output,
        **inputs,
    )

    root = ET.parse(output).getroot()
    resources = root.findall("./Resources/Resource")
    tracks = root.findall("./Panel/Track")
    assert count == len(resources) == len(tracks) == 5
    assert all((output.parent / item.attrib["path"]).is_file() for item in resources)
    assert [track.attrib["name"] for track in tracks] == [
        "CTX | mean ATAC Tn5 signal (background-TMM)",
        "CTX | mean H3K27ac fragment coverage (background-TMM)",
        "CTX | context DHSs (master coordinates)",
        "CTX | active regulatory elements",
        "Master DHS registry",
    ]
