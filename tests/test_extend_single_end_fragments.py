import pytest

from extend_single_end_fragments import extend_intervals


def test_extend_single_end_fragments_is_strand_aware_and_clips_boundaries():
    records = [
        "chr2L\t10\t60\tplus\t30\t+\n",
        "chr2L\t200\t250\tminus\t30\t-\n",
        "chr2L\t950\t1000\tright\t30\t+\n",
        "chr2L\t0\t50\tleft\t30\t-\n",
    ]

    observed = list(
        extend_intervals(
            records,
            chrom_sizes={"chr2L": 1000},
            fragment_length=150,
        )
    )

    assert observed == [
        "chr2L\t10\t160\tplus\t30\t+\n",
        "chr2L\t100\t250\tminus\t30\t-\n",
        "chr2L\t950\t1000\tright\t30\t+\n",
        "chr2L\t0\t50\tleft\t30\t-\n",
    ]


@pytest.mark.parametrize(
    ("record", "message"),
    [
        ("chrX\t0\t50\tread\t0\t+\n", "unknown chromosome"),
        ("chr2L\t0\t50\tread\t0\t.\n", "invalid strand"),
        ("chr2L\t0\t200\tread\t0\t+\n", "longer than fragment"),
    ],
)
def test_extend_single_end_fragments_rejects_invalid_records(record, message):
    with pytest.raises(ValueError, match=message):
        list(
            extend_intervals(
                [record],
                chrom_sizes={"chr2L": 1000},
                fragment_length=150,
            )
        )
