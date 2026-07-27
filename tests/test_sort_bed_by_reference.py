import subprocess
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "src" / "sort_bed_by_reference.sh"


def run_sorter(tmp_path, chrom_sizes, records):
    chrom_sizes_path = tmp_path / "chrom.sizes"
    chrom_sizes_path.write_text(chrom_sizes)
    return subprocess.run(
        ["bash", str(SCRIPT), str(chrom_sizes_path), str(tmp_path / "sort")],
        input=records,
        text=True,
        capture_output=True,
        check=False,
    )


def test_sort_bed_by_reference_uses_reference_then_coordinate_order(tmp_path):
    records = (
        "chr10\t5\t6\td\n"
        "chr2\t20\t25\tc\n"
        "chr1\t10\t11\tb\n"
        "chr2\t2\t3\ta\n"
        "chr1\t10\t11\tstable-first\n"
    )

    result = run_sorter(tmp_path, "chr1\t100\nchr2\t100\nchr10\t100\n", records)

    assert result.returncode == 0, result.stderr
    assert result.stdout == (
        "chr1\t10\t11\tb\n"
        "chr1\t10\t11\tstable-first\n"
        "chr2\t2\t3\ta\n"
        "chr2\t20\t25\tc\n"
        "chr10\t5\t6\td\n"
    )


def test_sort_bed_by_reference_rejects_unknown_chromosome(tmp_path):
    result = run_sorter(
        tmp_path,
        "chr1\t100\n",
        "chr2\t1\t2\n",
    )

    assert result.returncode != 0
    assert "unknown chromosome: chr2" in result.stderr


def test_sort_bed_by_reference_rejects_invalid_interval(tmp_path):
    result = run_sorter(
        tmp_path,
        "chr1\t100\n",
        "chr1\t5\t2\n",
    )

    assert result.returncode != 0
    assert "invalid BED record" in result.stderr
