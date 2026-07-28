import csv
import gzip
import json
from pathlib import Path

import pytest

from short_read_processing.activity_tmm import TMM_BACKGROUND_METHOD, TMM_FACTOR_FIELDS
from short_read_processing.catalog_visualization import (
    build_catalog_beds,
    mean_unionbedg_rows,
    read_track_factors,
)


def test_catalog_beds_preserve_membership_classes_and_guardrail_warning(tmp_path):
    master = tmp_path / "master.bed"
    summits = tmp_path / "summits.bed"
    matrix = tmp_path / "matrix.tsv"
    master.write_text(
        "chr1\t10\t40\tDHS0000001\t0\t.\n"
        "chr1\t50\t90\tDHS0000002\t0\t.\n"
    )
    summits.write_text(
        "chr1\t20\t21\tDHS0000001\t0\t.\n"
        "chr1\t70\t71\tDHS0000002\t0\t.\n"
    )
    matrix.write_text(
        "master_dhs_id\tchrom\tstart\tend\tsummit\tcontext_n\tctx\n"
        "DHS0000001\tchr1\t10\t40\t20\t1\t1\n"
        "DHS0000002\tchr1\t50\t90\t70\t0\t0\n"
    )
    active = tmp_path / "ctx.active.tsv.gz"
    fields = [
        "master_dhs_id",
        "chrom",
        "start",
        "end",
        "summit",
        "context",
        "context_membership",
        "regulatory_class",
        "mixture_component",
        "mixture_guardrail_warning",
        "mixture_high_posterior_probability",
    ]
    with gzip.open(active, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerow(
            {
                "master_dhs_id": "DHS0000001",
                "chrom": "chr1",
                "start": 10,
                "end": 40,
                "summit": 20,
                "context": "ctx",
                "context_membership": 1,
                "regulatory_class": "distal_enhancer_like",
                "mixture_component": "high",
                "mixture_guardrail_warning": 1,
                "mixture_high_posterior_probability": 0.75,
            }
        )
    output_master = tmp_path / "bed/master_dhs.bed"
    output_dhs = tmp_path / "bed/ctx.dhs.bed"
    output_active = tmp_path / "bed/ctx.active_elements.bed"
    output_manifest = tmp_path / "bed/bed_tracks.json"

    result = build_catalog_beds(
        master_bed=master,
        summit_bed=summits,
        context_matrix=matrix,
        active_paths={"ctx": active},
        output_master_bed=output_master,
        output_context_dhs={"ctx": output_dhs},
        output_active_beds={"ctx": output_active},
        output_manifest=output_manifest,
    )

    assert output_master.read_text() == master.read_text()
    assert output_dhs.read_text().splitlines() == [
        "chr1\t10\t40\tDHS0000001\t0\t.\t20\t21\t0,145,130"
    ]
    active_fields = output_active.read_text().strip().split("\t")
    assert active_fields[3] == (
        "DHS0000001|distal_enhancer_like|high_mixture_warning"
    )
    assert active_fields[4] == "750"
    assert active_fields[8] == "69,117,180"
    assert result["context_metrics"]["ctx"]["active_element_count"] == 1
    assert json.loads(output_manifest.read_text())["status"] == "ok"


def test_track_factors_and_unionbedg_mean_are_exact(tmp_path):
    factors = tmp_path / "factors.tsv"
    with factors.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=TMM_FACTOR_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for library_id, effective_size in (("rep1", 100.0), ("rep2", 200.0)):
            writer.writerow(
                {
                    "library_id": library_id,
                    "assay": "atac",
                    "context": "ctx",
                    "total_units": 100,
                    "feature_count": 10,
                    "tmm_normalization_factor": effective_size / 100,
                    "effective_library_size": effective_size,
                    "normalization_method": TMM_BACKGROUND_METHOD,
                }
            )

    parsed = read_track_factors(
        factor_path=factors,
        library_ids=["rep2", "rep1"],
        assay="atac",
        context="ctx",
    )
    assert list(parsed) == ["rep2", "rep1"]
    assert parsed["rep1"]["scale"] == 10_000
    assert parsed["rep2"]["scale"] == 5_000

    rows = list(
        mean_unionbedg_rows(
            ["chr1\t0\t10\t2\t0\n", "chr1\t10\t20\t2\t4\n"],
            chromosome_sizes=[("chr1", 20)],
            library_n=2,
        )
    )
    assert rows == [("chr1", 0, 10, 1.0), ("chr1", 10, 20, 3.0)]

    with pytest.raises(ValueError, match="unsorted or overlapping"):
        list(
            mean_unionbedg_rows(
                ["chr1\t10\t20\t1\n", "chr1\t5\t8\t1\n"],
                chromosome_sizes=[("chr1", 20)],
                library_n=1,
            )
        )
