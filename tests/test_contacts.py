import csv
import gzip
from pathlib import Path
import shutil

import pytest

import numpy as np

from short_read_processing.contacts import (
    balance_cooler_with_retry,
    select_source_resolution,
    standardize_context,
    valid_balanced_diagonal,
)
from short_read_processing.contact_metadata import CONTACT_SOURCE_COLUMNS


def test_select_source_resolution_uses_closest_exact_divisor():
    assert select_source_resolution([100, 200, 400, 800, 1600, 3200], 5000) == 200
    assert select_source_resolution([100, 1000, 2000, 4000], 5000) == 1000


def test_select_source_resolution_rejects_approximate_rebinning():
    with pytest.raises(ValueError, match="evenly divides"):
        select_source_resolution([4000], 5000)


def test_valid_balanced_diagonal_keeps_true_zeros_but_excludes_masked_bins():
    matrix = np.asarray(
        [
            [0.0, 2.0, 3.0, 4.0],
            [2.0, 0.0, 0.0, 5.0],
            [3.0, 0.0, 0.0, 6.0],
            [4.0, 5.0, 6.0, 0.0],
        ]
    )
    weights = np.asarray([1.0, 1.0, 1.0, np.nan])

    assert valid_balanced_diagonal(matrix, weights, 1).tolist() == [2.0, 0.0]


def test_valid_balanced_diagonal_rejects_nonpositive_offsets():
    with pytest.raises(ValueError, match="positive"):
        valid_balanced_diagonal(np.eye(2), np.ones(2), 0)


def test_balance_retries_nonconverged_matrix(monkeypatch):
    calls = []

    class FakeCooler:
        @staticmethod
        def balance_cooler(matrix, **kwargs):
            calls.append((matrix, kwargs))
            converged = len(calls) == 2
            return np.ones(3), {
                "converged": converged,
                "var": 1e-6 if converged else 1e-3,
            }

    monkeypatch.setattr(
        "short_read_processing.contacts._cooler", lambda: FakeCooler
    )

    weights, stats, attempts = balance_cooler_with_retry("matrix")

    assert weights.tolist() == [1.0, 1.0, 1.0]
    assert stats["converged"] is True
    assert [call[1]["max_iters"] for call in calls] == [200, 2_000]
    assert all(call[1]["store"] is True for call in calls)
    assert [attempt["converged"] for attempt in attempts] == [False, True]


def test_contact_environment_pins_atlas_compatible_pandas():
    environment = (
        Path(__file__).resolve().parents[1] / "workflow" / "envs" / "contacts.yaml"
    ).read_text(encoding="utf-8")

    assert "pandas=2.2" in environment


def test_standardize_context_merges_counts_before_balancing(tmp_path):
    cooler = pytest.importorskip("cooler")
    pandas = pytest.importorskip("pandas")
    bins = pandas.DataFrame(
        {
            "chrom": ["chr2L"] * 20,
            "start": list(range(0, 2000, 100)),
            "end": list(range(100, 2100, 100)),
        }
    )
    pixels = pandas.DataFrame(
        [
            {"bin1_id": first, "bin2_id": second, "count": 1}
            for first in range(20)
            for second in range(first, 20)
        ]
    )
    source_root = tmp_path / "repository"
    source_directory = source_root / "data" / "raw" / "contacts"
    source_directory.mkdir(parents=True)
    rows = []
    for replicate in (1, 2):
        source = tmp_path / f"rep{replicate}.cool"
        compressed = source_directory / f"rep{replicate}.cool.gz"
        cooler.create_cooler(str(source), bins, pixels, ordered=True)
        with source.open("rb") as input_handle, gzip.open(
            compressed, "wb"
        ) as output_handle:
            shutil.copyfileobj(input_handle, output_handle)
        rows.append(
            {
                "source_id": f"rep{replicate}",
                "context": "ctx",
                "assay": "Micro-C",
                "replicate": f"rep{replicate}",
                "format": "cool.gz",
                "url": f"https://example.org/rep{replicate}.cool.gz",
                "local_path": f"data/raw/contacts/rep{replicate}.cool.gz",
                "checksum": "",
                "match_quality": "test_exact",
                "biological_context": "synthetic",
                "caveat": "",
            }
        )
    manifest = tmp_path / "contacts.tsv"
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=CONTACT_SOURCE_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)

    output = tmp_path / "merged.cool"
    workdir = tmp_path / "work"
    metrics = standardize_context(
        context="ctx",
        source_manifest=manifest,
        repository_root=source_root,
        target_resolution=100,
        workdir=workdir,
        output=output,
    )

    merged = cooler.Cooler(str(output))
    assert merged.binsize == 100
    assert merged.pixels()[:]["count"].min() == 2
    assert metrics["replicate_count"] == 2
    assert metrics["finite_weight_fraction"] > 0
    assert not workdir.exists()
