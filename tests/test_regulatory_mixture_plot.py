import csv
import gzip
from pathlib import Path

from short_read_processing.activity import _tsv_content, sha256_file, write_deterministic_gzip
from short_read_processing.regulatory_elements import CATALOG_FIELDS, MIXTURE_FIELDS
from short_read_processing.regulatory_mixture_plot import (
    BIN_FIELDS,
    build_mixture_distribution_plot,
)


def test_mixture_plot_uses_positive_context_members_and_is_deterministic(tmp_path):
    catalog = tmp_path / "catalog.tsv.gz"
    catalog_rows = []
    for index, (member, signal) in enumerate(
        ((1, 1.0), (1, 2.0), (1, 20.0), (1, 40.0), (0, 400.0), (1, 0.0)),
        start=1,
    ):
        row = {field: "" for field in CATALOG_FIELDS}
        row.update(
            {
                "master_dhs_id": f"DHS{index:07d}",
                "chrom": "chr1",
                "start": index * 10,
                "end": index * 10 + 5,
                "summit": index * 10 + 2,
                "width_bp": 5,
                "context": "ctx",
                "context_membership": member,
                "h3k27ac_max_500_normalized_cpm_per_kb": signal,
            }
        )
        catalog_rows.append(row)
    write_deterministic_gzip(catalog, _tsv_content(CATALOG_FIELDS, catalog_rows))

    mixtures = tmp_path / "mixtures.tsv"
    mixture = {field: "" for field in MIXTURE_FIELDS}
    mixture.update(
        {
            "context": "ctx",
            "positive_member_n": 4,
            "mixture_supported": 1,
            "support_reason": "supported",
            "delta_bic": 20,
            "low_mean_log10": 0.15,
            "high_mean_log10": 1.45,
            "low_sd_log10": 0.2,
            "high_sd_log10": 0.2,
            "low_weight": 0.5,
            "high_weight": 0.5,
            "ashman_d": 4.0,
            "posterior_crossing_log10": 0.8,
        }
    )
    mixtures.write_text(_tsv_content(MIXTURE_FIELDS, [mixture]))
    svg = tmp_path / "plot.svg"
    bins = tmp_path / "bins.tsv.gz"
    metrics = tmp_path / "metrics.json"

    result = build_mixture_distribution_plot(
        catalog_path=catalog,
        mixture_path=mixtures,
        output_svg=svg,
        output_bins=bins,
        output_metrics=metrics,
        bin_n=20,
    )
    hashes = (sha256_file(svg), sha256_file(bins), sha256_file(metrics))
    build_mixture_distribution_plot(
        catalog_path=catalog,
        mixture_path=mixtures,
        output_svg=svg,
        output_bins=bins,
        output_metrics=metrics,
        bin_n=20,
    )

    assert result["supported_contexts"] == ["ctx"]
    assert hashes == (sha256_file(svg), sha256_file(bins), sha256_file(metrics))
    assert "posterior crossing" in svg.read_text()
    with gzip.open(bins, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assert reader.fieldnames == BIN_FIELDS
        rows = list(reader)
    assert len(rows) == 20
    assert {row["context"] for row in rows} == {"ctx"}


def test_mixture_plot_reports_guard_failure_without_a_fit(tmp_path):
    catalog = tmp_path / "catalog.tsv.gz"
    catalog_rows = []
    for index, signal in enumerate((0.0, 2.0), start=1):
        row = {field: "" for field in CATALOG_FIELDS}
        row.update(
            {
                "master_dhs_id": f"DHS{index:07d}",
                "chrom": "chr1",
                "start": index * 10,
                "end": index * 10 + 5,
                "summit": index * 10 + 2,
                "width_bp": 5,
                "context": "small_context",
                "context_membership": 1,
                "h3k27ac_max_500_normalized_cpm_per_kb": signal,
            }
        )
        catalog_rows.append(row)
    write_deterministic_gzip(catalog, _tsv_content(CATALOG_FIELDS, catalog_rows))

    mixtures = tmp_path / "mixtures.tsv"
    mixture = {field: "" for field in MIXTURE_FIELDS}
    mixture.update(
        {
            "context": "small_context",
            "positive_member_n": 1,
            "mixture_supported": 0,
            "support_reason": "insufficient_positive_members",
        }
    )
    mixtures.write_text(_tsv_content(MIXTURE_FIELDS, [mixture]))

    result = build_mixture_distribution_plot(
        catalog_path=catalog,
        mixture_path=mixtures,
        output_svg=tmp_path / "plot.svg",
        output_bins=tmp_path / "bins.tsv.gz",
        output_metrics=tmp_path / "metrics.json",
        bin_n=20,
    )

    assert result["supported_contexts"] == []
    assert result["unsupported_contexts"] == ["small_context"]
    assert "WARNING: insufficient_positive_members" in (
        tmp_path / "plot.svg"
    ).read_text()
