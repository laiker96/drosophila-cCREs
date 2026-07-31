import csv
import gzip
import json
from pathlib import Path

import pytest

from short_read_processing.contact_links import (
    DISTANCE_MODEL_GENE_FIELDS,
    EDGE_FIELDS,
    Element,
    NODE_FIELDS,
    _active_contact_enhancer_gene_rows,
    _active_distance_enhancer_gene_rows,
    _atomic_tsv,
    build_active_contact_enhancer_gene_candidates,
    build_active_distance_enhancer_gene_candidates,
    build_context_links,
    build_promoter_table,
    read_context_elements,
)


CATALOG_FIELDS = (
    "master_dhs_id",
    "chrom",
    "start",
    "end",
    "summit",
    "context",
    "context_membership",
    "regulatory_class",
    "atac_normalized_cpm_per_kb",
    "mixture_high_posterior_probability",
    "activity_state",
    "combined_activity_max_500",
    "blacklist_overlap",
)


def test_compressed_graph_tables_are_byte_deterministic(tmp_path):
    first = tmp_path / "first.tsv.gz"
    second = tmp_path / "second.tsv.gz"

    _atomic_tsv(first, ("node", "score"), [{"node": "x", "score": 1}])
    _atomic_tsv(second, ("node", "score"), [{"node": "x", "score": 1}])

    assert first.read_bytes() == second.read_bytes()


def _read_tsv(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _edge(*, gene_id, promoter_id, enrichment, score, contact_weight=1.0):
    return {
        "context": "ctx",
        "gene_id": gene_id,
        "gene_name": gene_id,
        "promoter_id": promoter_id,
        "distance_bp": 100,
        "promoter_active": 1,
        "promoter_atac_signal": "2",
        "promoter_h3k27ac_posterior": "0.8",
        "contact_strategy": "observed",
        "contact_assay": "Micro-C",
        "contact_match": "matched",
        "contact_resolution_bp": 5000,
        "contact_status": "observed_matrix_pixel",
        "observed_balanced_contact": "1",
        "powerlaw_expected_contact": "1",
        "contact_weight": str(contact_weight),
        "observed_over_expected": enrichment,
        "candidate_link_score": str(score),
    }


def test_active_contact_enhancer_gene_projection_applies_both_thresholds():
    element = Element(
        master_dhs_id="DHS_active",
        chrom="chr2L",
        start=100,
        end=200,
        summit=150,
        regulatory_class="distal_enhancer_like",
        activity_state="active_high_mixture",
        atac_signal=4,
        h3k27ac_posterior=0.5,
        combined_activity=4,
        blacklist_overlap=0,
    )
    rows = _active_contact_enhancer_gene_rows(
        element,
        [
            _edge(
                gene_id="gene_x",
                promoter_id="promoter_below",
                enrichment="0.99",
                score=20,
            ),
            _edge(
                gene_id="gene_x",
                promoter_id="promoter_qualifying",
                enrichment="1",
                score=10,
            ),
            _edge(
                gene_id="gene_y",
                promoter_id="promoter_modelled",
                enrichment="",
                score=30,
            ),
        ],
        element_posterior_threshold=0.5,
        observed_over_expected_threshold=1.0,
    )

    assert [row["gene_id"] for row in rows] == ["gene_x"]
    assert rows[0]["best_promoter_id"] == "promoter_qualifying"
    assert rows[0]["best_observed_over_expected"] == "1"

    below_threshold = Element(
        **{**element.__dict__, "h3k27ac_posterior": 0.499}
    )
    assert not _active_contact_enhancer_gene_rows(
        below_threshold,
        [_edge(gene_id="gene_x", promoter_id="promoter", enrichment="2", score=1)],
        element_posterior_threshold=0.5,
        observed_over_expected_threshold=1.0,
    )


def test_focused_candidate_table_is_derived_from_existing_graph_tables(tmp_path):
    nodes = tmp_path / "nodes.tsv.gz"
    edges = tmp_path / "edges.tsv.gz"
    output = tmp_path / "candidates.tsv.gz"
    metrics = tmp_path / "candidates.metrics.json"
    _atomic_tsv(
        nodes,
        NODE_FIELDS,
        [
            {
                "context": "ctx",
                "node_id": "element:DHS_active",
                "node_type": "element",
                "master_dhs_id": "DHS_active",
                "chrom": "chr2L",
                "start": 100,
                "end": 200,
                "anchor": 150,
                "regulatory_class": "distal_enhancer_like",
                "atac_signal": 4,
                "h3k27ac_posterior": 0.5,
                "combined_activity": 4,
                "blacklist_overlap": 0,
            },
            {
                "context": "ctx",
                "node_id": "promoter:P1",
                "node_type": "promoter",
                "promoter_id": "P1",
                "chrom": "chr2L",
                "start": 900,
                "end": 1100,
                "anchor": 1000,
                "regulatory_class": "promoter",
                "atac_signal": 2,
                "h3k27ac_posterior": 0.8,
                "combined_activity": 2,
            },
        ],
    )
    _atomic_tsv(
        edges,
        EDGE_FIELDS,
        [
            _edge(
                gene_id="gene_x",
                promoter_id="promoter_below",
                enrichment="0.99",
                score=20,
            )
            | {"master_dhs_id": "DHS_active"},
            _edge(
                gene_id="gene_x",
                promoter_id="promoter_qualifying",
                enrichment="1",
                score=10,
            )
            | {"master_dhs_id": "DHS_active"},
            _edge(
                gene_id="gene_y",
                promoter_id="promoter_modelled",
                enrichment="",
                score=30,
            )
            | {"master_dhs_id": "DHS_active"},
        ],
    )

    result = build_active_contact_enhancer_gene_candidates(
        context="ctx",
        nodes_path=nodes,
        edges_path=edges,
        element_posterior_threshold=0.5,
        observed_over_expected_threshold=1.0,
        output=output,
        metrics_output=metrics,
    )

    rows = _read_tsv(output)
    assert [row["gene_id"] for row in rows] == ["gene_x"]
    assert rows[0]["best_promoter_id"] == "promoter_qualifying"
    assert result["active_contact_enhancer_gene_candidate_count"] == 1
    assert result["element_with_candidate_count"] == 1


def _distance_edge(
    *,
    gene_id,
    promoter_id,
    distance,
    score,
    promoter_active=1,
    contact_weight=1.0,
):
    return _edge(
        gene_id=gene_id,
        promoter_id=promoter_id,
        enrichment="",
        score=score,
        contact_weight=contact_weight,
    ) | {
        "distance_bp": distance,
        "promoter_active": promoter_active,
        "contact_strategy": "powerlaw",
        "contact_assay": "distance_model",
        "contact_match": "no_exact_map",
        "contact_status": "distance_model",
        "observed_balanced_contact": "",
        "observed_over_expected": "",
    }


def test_distance_projection_ranks_active_promoters_and_reports_tss_baselines():
    element = Element(
        master_dhs_id="DHS_active",
        chrom="chr2L",
        start=100,
        end=200,
        summit=150,
        regulatory_class="distal_enhancer_like",
        activity_state="active_high_mixture",
        atac_signal=4,
        h3k27ac_posterior=0.5,
        combined_activity=4,
        blacklist_overlap=0,
    )
    rows = _active_distance_enhancer_gene_rows(
        element,
        [
            _distance_edge(
                gene_id="gene_x",
                promoter_id="promoter_x",
                distance=500,
                score=10,
            ),
            _distance_edge(
                gene_id="gene_y",
                promoter_id="promoter_y",
                distance=50,
                score=9,
            ),
            _distance_edge(
                gene_id="gene_z",
                promoter_id="promoter_z",
                distance=10,
                score=30,
                promoter_active=0,
            ),
        ],
        element_posterior_threshold=0.5,
    )

    assert [row["gene_id"] for row in rows] == ["gene_x", "gene_y"]
    assert rows[0]["is_primary_candidate"] == 1
    assert rows[1]["is_primary_candidate"] == 0
    assert {row["evidence_type"] for row in rows} == {
        "distance_model_active_promoter"
    }
    assert {row["nearest_active_tss_gene_ids"] for row in rows} == {"gene_y"}
    assert {row["nearest_active_tss_distance_bp"] for row in rows} == {"50"}
    assert {row["nearest_tss_gene_ids"] for row in rows} == {"gene_z"}
    assert {row["nearest_tss_distance_bp"] for row in rows} == {"10"}
    assert {row["best_observed_over_expected"] for row in rows} == {""}


def test_distance_candidate_table_is_derived_from_existing_graph_tables(tmp_path):
    nodes = tmp_path / "nodes.tsv.gz"
    edges = tmp_path / "edges.tsv.gz"
    output = tmp_path / "candidates.tsv.gz"
    metrics = tmp_path / "candidates.metrics.json"
    _atomic_tsv(
        nodes,
        NODE_FIELDS,
        [
            {
                "context": "ctx",
                "node_id": "element:DHS_active",
                "node_type": "element",
                "master_dhs_id": "DHS_active",
                "chrom": "chr2L",
                "start": 100,
                "end": 200,
                "anchor": 150,
                "regulatory_class": "distal_enhancer_like",
                "atac_signal": 4,
                "h3k27ac_posterior": 0.5,
                "combined_activity": 4,
                "blacklist_overlap": 0,
            }
        ],
    )
    _atomic_tsv(
        edges,
        EDGE_FIELDS,
        [
            _distance_edge(
                gene_id="gene_x",
                promoter_id="promoter_x",
                distance=100,
                score=2,
            )
            | {"master_dhs_id": "DHS_active"},
            _distance_edge(
                gene_id="gene_y",
                promoter_id="promoter_y",
                distance=50,
                score=1,
                promoter_active=0,
            )
            | {"master_dhs_id": "DHS_active"},
        ],
    )

    result = build_active_distance_enhancer_gene_candidates(
        context="ctx",
        nodes_path=nodes,
        edges_path=edges,
        element_posterior_threshold=0.5,
        output=output,
        metrics_output=metrics,
    )

    rows = _read_tsv(output)
    assert tuple(rows[0]) == DISTANCE_MODEL_GENE_FIELDS
    assert [row["gene_id"] for row in rows] == ["gene_x"]
    assert rows[0]["nearest_tss_gene_ids"] == "gene_y"
    assert rows[0]["best_observed_balanced_contact"] == ""
    assert result["active_distance_enhancer_gene_candidate_count"] == 1
    assert result["element_with_candidate_count"] == 1


def test_legacy_catalog_defaults_missing_blacklist_overlap_to_zero(tmp_path):
    elements = tmp_path / "legacy.tsv"
    legacy_fields = tuple(
        field for field in CATALOG_FIELDS if field != "blacklist_overlap"
    )
    row = {
        "master_dhs_id": "DHS_legacy",
        "chrom": "chr2L",
        "start": 100,
        "end": 200,
        "summit": 150,
        "context": "ctx",
        "context_membership": 1,
        "regulatory_class": "distal_enhancer_like",
        "atac_normalized_cpm_per_kb": 4,
        "mixture_high_posterior_probability": 0.8,
        "activity_state": "active_high_mixture",
        "combined_activity_max_500": 4,
    }
    with elements.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=legacy_fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerow(row)

    parsed, defaulted = read_context_elements(elements, "ctx")

    assert defaulted is True
    assert parsed[0].blacklist_overlap == 0


def test_legacy_catalog_still_requires_non_blacklist_fields(tmp_path):
    elements = tmp_path / "invalid.tsv"
    fields = tuple(
        field
        for field in CATALOG_FIELDS
        if field not in {"blacklist_overlap", "summit"}
    )
    elements.write_text("\t".join(fields) + "\n", encoding="utf-8")

    with pytest.raises(ValueError, match="summit"):
        read_context_elements(elements, "ctx")


def test_powerlaw_links_prioritize_an_accessible_active_promoter(tmp_path):
    annotation = tmp_path / "genes.gtf"
    annotation.write_text(
        'chr2L\ttest\ttranscript\t1001\t1300\t.\t+\t.\tgene_id "gene_x"; '
        'gene_name "x"; transcript_id "tx1";\n'
        'chr2L\ttest\ttranscript\t5001\t5300\t.\t+\t.\tgene_id "gene_y"; '
        'gene_name "y"; transcript_id "ty1";\n',
        encoding="utf-8",
    )
    chrom_sizes = tmp_path / "chrom.sizes"
    chrom_sizes.write_text("chr2L\t10000\n", encoding="utf-8")
    promoters = tmp_path / "promoters.tsv.gz"
    promoter_metrics = tmp_path / "promoters.json"
    build_promoter_table(
        annotation=annotation,
        chrom_sizes=chrom_sizes,
        canonical_chromosomes=["chr2L"],
        promoter_width=500,
        annotation_checksum=None,
        output=promoters,
        metrics_output=promoter_metrics,
    )

    elements = tmp_path / "ctx.elements.tsv.gz"
    rows = [
        {
            "master_dhs_id": "DHS_element",
            "chrom": "chr2L",
            "start": 1900,
            "end": 2000,
            "summit": 1950,
            "context": "ctx",
            "context_membership": 1,
            "regulatory_class": "distal_enhancer_like",
            "atac_normalized_cpm_per_kb": 4,
            "mixture_high_posterior_probability": 0.8,
            "activity_state": "active_high_mixture",
            "combined_activity_max_500": 4,
            "blacklist_overlap": 0,
        },
        {
            "master_dhs_id": "DHS_promoter_x",
            "chrom": "chr2L",
            "start": 900,
            "end": 1100,
            "summit": 1000,
            "context": "ctx",
            "context_membership": 1,
            "regulatory_class": "promoter_associated",
            "atac_normalized_cpm_per_kb": 10,
            "mixture_high_posterior_probability": 0.9,
            "activity_state": "active_high_mixture",
            "combined_activity_max_500": 10,
            "blacklist_overlap": 0,
        },
    ]
    with gzip.open(elements, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=CATALOG_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    powerlaw = tmp_path / "powerlaw.json"
    powerlaw.write_text(
        json.dumps(
            {
                "contexts": {},
                "atlas_powerlaw": {"gamma": -1.0, "scale": 1000.0},
            }
        ),
        encoding="utf-8",
    )
    nodes = tmp_path / "nodes.tsv.gz"
    edges = tmp_path / "edges.tsv.gz"
    genes = tmp_path / "genes.tsv.gz"
    metrics = tmp_path / "metrics.json"

    result = build_context_links(
        context="ctx",
        context_elements=elements,
        promoters_path=promoters,
        powerlaw_path=powerlaw,
        strategy="powerlaw",
        contact_assay="distance_model",
        contact_match="no_exact_map",
        configured_resolution=100,
        maximum_distance=10_000,
        pseudocount_fraction=0.01,
        posterior_threshold=0.5,
        observed_contact_path=None,
        nodes_output=nodes,
        edges_output=edges,
        genes_output=genes,
        metrics_output=metrics,
    )

    element_genes = [
        row for row in _read_tsv(genes) if row["master_dhs_id"] == "DHS_element"
    ]
    element_edges = [
        row for row in _read_tsv(edges) if row["master_dhs_id"] == "DHS_element"
    ]
    by_gene = {row["gene_id"]: row for row in element_genes}
    assert by_gene["gene_x"]["best_promoter_active"] == "1"
    assert by_gene["gene_x"]["candidate_gene_rank"] == "1"
    assert by_gene["gene_y"]["best_promoter_active"] == "0"
    assert by_gene["gene_y"]["active_candidate_gene_rank"] == ""
    assert all(row["contact_strategy"] == "powerlaw" for row in element_genes)
    assert element_edges[0]["source_node_id"] == "element:DHS_element"
    assert "element_chrom" not in element_edges[0]
    assert result["element_node_count"] == 2
    assert result["promoter_node_count"] == 2
    assert result["active_promoter_count"] == 1


def test_promoter_activity_accepts_overlapping_master_dhs_boundaries(tmp_path):
    annotation = tmp_path / "genes.gtf"
    annotation.write_text(
        'chr2L\ttest\ttranscript\t1001\t1300\t.\t+\t.\tgene_id "gene_x"; '
        'gene_name "x"; transcript_id "tx1";\n',
        encoding="utf-8",
    )
    chrom_sizes = tmp_path / "chrom.sizes"
    chrom_sizes.write_text("chr2L\t10000\n", encoding="utf-8")
    promoters = tmp_path / "promoters.tsv.gz"
    build_promoter_table(
        annotation=annotation,
        chrom_sizes=chrom_sizes,
        canonical_chromosomes=["chr2L"],
        promoter_width=500,
        annotation_checksum=None,
        output=promoters,
        metrics_output=tmp_path / "promoters.json",
    )

    # The summit-aware master registry may retain distinct sites with overlapping
    # boundaries, so promoter annotation must not assume disjoint intervals.
    elements = tmp_path / "ctx.elements.tsv.gz"
    rows = []
    rows.append(
        {
            "master_dhs_id": "DHS_other",
            "chrom": "chr2L",
            "start": 100,
            "end": 200,
            "summit": 150,
            "context": "other",
            "context_membership": 0,
            "regulatory_class": "distal",
            "atac_normalized_cpm_per_kb": 0,
            "mixture_high_posterior_probability": "",
            "activity_state": "not_accessible",
            "combined_activity_max_500": 0,
            "blacklist_overlap": 0,
        }
    )
    for identifier, start, end, summit in (
        ("DHS_1", 850, 1050, 950),
        ("DHS_2", 950, 1150, 1050),
    ):
        rows.append(
            {
                "master_dhs_id": identifier,
                "chrom": "chr2L",
                "start": start,
                "end": end,
                "summit": summit,
                "context": "ctx",
                "context_membership": 1,
                "regulatory_class": "promoter_associated",
                "atac_normalized_cpm_per_kb": 2,
                "mixture_high_posterior_probability": 0.8,
                "activity_state": "active_high_mixture",
                "combined_activity_max_500": 2,
                "blacklist_overlap": 0,
            }
        )
    with gzip.open(elements, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=CATALOG_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    powerlaw = tmp_path / "powerlaw.json"
    powerlaw.write_text(
        json.dumps(
            {
                "contexts": {},
                "atlas_powerlaw": {"gamma": -1.0, "scale": 1000.0},
            }
        ),
        encoding="utf-8",
    )

    build_context_links(
        context="ctx",
        context_elements=elements,
        promoters_path=promoters,
        powerlaw_path=powerlaw,
        strategy="powerlaw",
        contact_assay="distance_model",
        contact_match="no_exact_map",
        configured_resolution=100,
        maximum_distance=10_000,
        pseudocount_fraction=0.01,
        posterior_threshold=0.5,
        observed_contact_path=None,
        nodes_output=tmp_path / "nodes.tsv.gz",
        edges_output=tmp_path / "edges.tsv.gz",
        genes_output=tmp_path / "genes.tsv.gz",
        metrics_output=tmp_path / "metrics.json",
    )

    promoter_node = next(
        row for row in _read_tsv(tmp_path / "nodes.tsv.gz") if row["node_type"] == "promoter"
    )
    assert set(promoter_node["master_dhs_id"].split(";")) == {"DHS_1", "DHS_2"}
