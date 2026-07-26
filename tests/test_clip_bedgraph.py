from pathlib import Path

import pytest

from clip_bedgraph import clip_bedgraph


def test_clip_bedgraph_clips_reference_boundaries_and_discards_outside(tmp_path):
    sizes = tmp_path / "chrom.sizes"
    sizes.write_text("chr1\t10\n", encoding="utf-8")
    source = tmp_path / "input.bdg"
    source.write_text(
        "track type=bedGraph\n"
        "chr1\t-5\t3\t1.0\n"
        "chr1\t3\t7\t2.0\n"
        "chr1\t8\t12\t3.0\n"
        "chr1\t12\t15\t4.0\n",
        encoding="utf-8",
    )
    output = tmp_path / "output.bdg"

    stats = clip_bedgraph(source, output, sizes)

    assert output.read_text(encoding="utf-8") == (
        "chr1\t0\t3\t1.0\n"
        "chr1\t3\t7\t2.0\n"
        "chr1\t8\t10\t3.0\n"
    )
    assert stats == {
        "records": 4,
        "written": 3,
        "left_clipped": 1,
        "right_clipped": 2,
        "discarded_outside": 1,
        "maximum_right_clip": 5,
    }


@pytest.mark.parametrize(
    ("line", "message"),
    [
        ("chr2\t0\t1\t1\n", "unknown chromosome"),
        ("chr1\t4\t4\t1\n", "non-positive interval width"),
        ("chr1\t0\t1\tnot-a-number\n", "invalid bedGraph value"),
    ],
)
def test_clip_bedgraph_rejects_unfixable_records(tmp_path, line, message):
    sizes = tmp_path / "chrom.sizes"
    sizes.write_text("chr1\t10\n", encoding="utf-8")
    source = tmp_path / "input.bdg"
    source.write_text(line, encoding="utf-8")
    output = tmp_path / "output.bdg"

    with pytest.raises(ValueError, match=message):
        clip_bedgraph(source, output, sizes)

    assert not output.exists()
