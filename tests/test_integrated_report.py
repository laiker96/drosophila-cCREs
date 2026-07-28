import json
from pathlib import Path

from short_read_processing.activity import _tsv_content, sha256_file
from short_read_processing.activity_tmm import TMM_FACTOR_FIELDS
from short_read_processing.integrated_report import (
    build_integrated_qc_report,
    discover_report_source_files,
)
from short_read_processing.regulatory_elements import MIXTURE_FIELDS, SUMMARY_FIELDS


def test_integrated_report_discovers_sources_and_writes_deterministic_outputs(tmp_path):
    upstream = tmp_path / "upstream"
    (upstream / "qc").mkdir(parents=True)
    qc = upstream / "qc" / "metrics.json"
    qc.write_text(
        json.dumps(
            [
                {
                    "sample": "atac_1",
                    "assay": "atac",
                    "layout": "paired",
                    "mapped_alignments": 1000,
                    "frip": 0.42,
                    "peak_count": 50,
                }
            ]
        )
    )
    discovered = discover_report_source_files([upstream])
    assert [record["kind"] for record in discovered] == ["qc_metrics"]

    master = tmp_path / "master.json"
    master.write_text(
        json.dumps(
            {
                "master_dhs_count": 100,
                "source_peak_count": 160,
                "multi_context_master_dhs_count": 40,
                "master_width_min": 50,
                "master_width_median": 200,
                "master_width_max": 400,
                "context_peak_counts": {"ctx": 80},
            }
        )
    )
    activity_metrics = tmp_path / "activity.json"
    activity_metrics.write_text(
        json.dumps(
            {
                "master_dhs_count": 100,
                "context_count": 1,
                "library_count": 4,
                "normalization_method": "tmm_background_10kb_v1",
            }
        )
    )
    catalog_metrics = tmp_path / "catalog.json"
    catalog_metrics.write_text(
        json.dumps(
            {"master_dhs_count": 100, "context_count": 1, "catalog_row_count": 100}
        )
    )
    factors = tmp_path / "factors.tsv"
    factor_rows = []
    for library_id, assay in (
        ("atac_1", "atac"),
        ("atac_2", "atac"),
        ("h3_1", "h3k27ac"),
        ("h3_2", "h3k27ac"),
    ):
        factor_rows.append(
            {
                "library_id": library_id,
                "assay": assay,
                "context": "ctx",
                "total_units": 1000,
                "feature_count": 20,
                "tmm_normalization_factor": 1.0,
                "effective_library_size": 1000,
                "normalization_method": "tmm_background_10kb_v1",
            }
        )
    factors.write_text(_tsv_content(TMM_FACTOR_FIELDS, factor_rows))

    mixtures = tmp_path / "mixtures.tsv"
    mixture = {field: "" for field in MIXTURE_FIELDS}
    mixture.update(
        {
            "context": "ctx",
            "positive_member_n": 100,
            "mixture_supported": 0,
            "support_reason": "insufficient_positive_members",
        }
    )
    mixtures.write_text(_tsv_content(MIXTURE_FIELDS, [mixture]))
    summary = tmp_path / "summary.tsv"
    summary.write_text(
        _tsv_content(
            SUMMARY_FIELDS,
            [
                {
                    "context": "ctx",
                    "regulatory_class": "distal_enhancer_like",
                    "mixture_component": "high",
                    "mixture_supported": 0,
                    "guardrail_failures": "insufficient_positive_members",
                    "element_n": 10,
                }
            ],
        )
    )
    mixture_svg = tmp_path / "mixtures.svg"
    mixture_svg.write_text('<svg xmlns="http://www.w3.org/2000/svg"></svg>')

    config = {
        "project": "atlas",
        "run_id": "catalog-v1",
        "assay": "activity",
        "reference": {"name": "dm6"},
        "activity": {
            "contexts": ["ctx"],
            "master": {"stats_json": str(master)},
            "libraries": [
                {
                    "id": library_id,
                    "context": "ctx",
                    "assay": assay,
                    "layout": "paired",
                    "qc_status": "accepted",
                }
                for library_id, assay in (
                    ("atac_1", "atac"),
                    ("atac_2", "atac"),
                    ("h3_1", "h3k27ac"),
                    ("h3_2", "h3k27ac"),
                )
            ],
        },
        "provenance": {},
    }
    current = {
        "activity_metrics": activity_metrics,
        "regulatory_element_metrics": catalog_metrics,
        "normalization_factors": factors,
        "mixture_models": mixtures,
        "regulatory_element_summary": summary,
        "mixture_distributions": mixture_svg,
    }
    output_html = tmp_path / "report.html"
    output_pdf = tmp_path / "report.pdf"
    output_metrics = tmp_path / "report.json"

    def fake_pdf_renderer(html_text: str, path: Path) -> str:
        assert "Integrated atlas regulatory-catalog QC report" in html_text
        path.write_bytes(b"%PDF-1.4\n% test\n")
        return "test-renderer"

    result = build_integrated_qc_report(
        config=config,
        source_files=[
            {
                **discovered[0],
                "path": str(discovered[0]["path"]),
                "source_root": str(discovered[0]["source_root"]),
            }
        ],
        current_files=current,
        output_html=output_html,
        output_pdf=output_pdf,
        output_metrics=output_metrics,
        pdf_renderer=fake_pdf_renderer,
    )
    first_hashes = tuple(
        sha256_file(path) for path in (output_html, output_pdf, output_metrics)
    )
    build_integrated_qc_report(
        config=config,
        source_files=[
            {
                **discovered[0],
                "path": str(discovered[0]["path"]),
                "source_root": str(discovered[0]["source_root"]),
            }
        ],
        current_files=current,
        output_html=output_html,
        output_pdf=output_pdf,
        output_metrics=output_metrics,
        pdf_renderer=fake_pdf_renderer,
    )

    assert result["mixture_warning_count"] == 1
    assert "insufficient_positive_members" in output_html.read_text()
    assert first_hashes == tuple(
        sha256_file(path) for path in (output_html, output_pdf, output_metrics)
    )
