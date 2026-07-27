import csv
import gzip
import json
from pathlib import Path

from short_read_processing.activity import (
    build_activity_outputs,
    read_master_elements,
    sha256_file,
    write_library_signal,
)
from short_read_processing.activity_qc import build_activity_qc_outputs


def _activity_tables(tmp_path: Path):
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
        ("atlas_atac_1", "atac", "atlas", "ctx", [1, 2], 100),
        ("atlas_atac_2", "atac", "atlas", "ctx", [2, 4], 200),
        ("atlas_h3_1", "h3k27ac", "atlas", "ctx", [3, 6], 300),
        ("ref_atac_1", "atac", "reference", "ref", [2, 6], 100),
        ("ref_atac_2", "atac", "reference", "ref", [4, 2], 200),
        ("ref_h3_1", "h3k27ac", "reference", "ref", [1, 8], 100),
        ("ref_h3_2", "h3k27ac", "reference", "ref", [2, 4], 200),
    ]
    signal_paths = {}
    for library_id, assay, cohort, context, counts, total in specifications:
        signal = tmp_path / "signals" / f"{library_id}.tsv.gz"
        summary = tmp_path / "signals" / f"{library_id}.json"
        write_library_signal(
            elements=elements,
            counts=counts,
            total_units=total,
            library_id=library_id,
            assay=assay,
            cohort=cohort,
            context=context,
            output=signal,
            summary=summary,
        )
        signal_paths[library_id] = signal
    outputs = {
        "output_library_signal": tmp_path / "activity" / "library.tsv.gz",
        "output_context_signal": tmp_path / "activity" / "context.tsv.gz",
        "output_reference": tmp_path / "activity" / "reference.tsv.gz",
        "output_activity": tmp_path / "activity" / "activity.tsv.gz",
        "output_context_views": {"ctx": tmp_path / "activity" / "ctx.tsv.gz"},
        "output_metrics": tmp_path / "activity" / "metrics.json",
        "output_provenance": tmp_path / "activity" / "provenance.json",
    }
    build_activity_outputs(
        signal_paths=signal_paths,
        atlas_contexts=["ctx"],
        reference_context="ref",
        provenance={"test": True},
        **outputs,
    )
    return outputs


def test_activity_qc_outputs_are_deterministic_and_descriptive(tmp_path):
    activity = _activity_tables(tmp_path)
    outputs = {
        "output_correlations": tmp_path / "qc" / "correlations.tsv",
        "output_distributions": tmp_path / "qc" / "distributions.tsv",
        "output_metrics": tmp_path / "qc" / "metrics.json",
        "output_report": tmp_path / "qc" / "report.html",
    }

    metrics = build_activity_qc_outputs(
        library_signal=activity["output_library_signal"],
        context_signal=activity["output_context_signal"],
        qnorm_reference=activity["output_reference"],
        activity_table=activity["output_activity"],
        activity_provenance=activity["output_provenance"],
        atlas_contexts=["ctx"],
        reference_context="ref",
        **outputs,
    )
    hashes = {name: sha256_file(path) for name, path in outputs.items()}
    build_activity_qc_outputs(
        library_signal=activity["output_library_signal"],
        context_signal=activity["output_context_signal"],
        qnorm_reference=activity["output_reference"],
        activity_table=activity["output_activity"],
        activity_provenance=activity["output_provenance"],
        atlas_contexts=["ctx"],
        reference_context="ref",
        **outputs,
    )

    assert metrics["status"] == "descriptive_qc_complete"
    assert metrics["automatic_acceptance_thresholds"] is False
    assert metrics["master_dhs_count"] == 2
    assert metrics["library_count"] == 7
    assert len(metrics["replicate_correlations"]) == 3
    assert metrics["groups_without_replicates"] == [
        {
            "cohort": "atlas",
            "context": "ctx",
            "assay": "h3k27ac",
            "library_ids": ["atlas_h3_1"],
        }
    ]
    assert {name: sha256_file(path) for name, path in outputs.items()} == hashes
    assert "Manual review required" in outputs["output_report"].read_text()
    assert "ATAC — post_qnorm" in outputs["output_report"].read_text()
    with outputs["output_distributions"].open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 6
    assert {(row["assay"], row["stage"]) for row in rows} >= {
        ("atac", "reference"),
        ("atac", "pre_qnorm"),
        ("atac", "post_qnorm"),
    }
    assert json.loads(outputs["output_metrics"].read_text()) == metrics


def test_activity_qc_rejects_context_mismatch(tmp_path):
    activity = _activity_tables(tmp_path)

    try:
        build_activity_qc_outputs(
            library_signal=activity["output_library_signal"],
            context_signal=activity["output_context_signal"],
            qnorm_reference=activity["output_reference"],
            activity_table=activity["output_activity"],
            activity_provenance=activity["output_provenance"],
            atlas_contexts=["missing"],
            reference_context="ref",
            output_correlations=tmp_path / "qc" / "correlations.tsv",
            output_distributions=tmp_path / "qc" / "distributions.tsv",
            output_metrics=tmp_path / "qc" / "metrics.json",
            output_report=tmp_path / "qc" / "report.html",
        )
    except ValueError as error:
        assert "context/assay set" in str(error)
    else:
        raise AssertionError("context mismatch was accepted")
