import csv
import gzip
import json
from pathlib import Path

from short_read_processing.contact_links import (
    _atomic_tsv,
    build_context_links,
    build_promoter_table,
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
