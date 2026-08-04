"""Aggregate contact-link metrics and freeze their provenance."""

from pathlib import Path

from short_read_processing.contact_links import aggregate_link_metrics


metrics, _provenance = aggregate_link_metrics(
    context_metric_paths=[Path(str(path)) for path in snakemake.input.context_metrics],
    candidate_metric_paths=[
        Path(str(path)) for path in snakemake.input.candidate_metrics
    ],
    nearest_candidate_metric_paths=[
        Path(str(path)) for path in snakemake.input.nearest_candidate_metrics
    ],
    distance_candidate_metric_paths=[
        Path(str(path)) for path in snakemake.input.distance_candidate_metrics
    ],
    source_manifest=Path(str(snakemake.input.manifest)),
    promoter_metrics=Path(str(snakemake.input.promoter_metrics)),
    contact_metrics=[Path(str(path)) for path in snakemake.input.contact_metrics],
    powerlaw=Path(str(snakemake.input.powerlaw)),
    output_metrics=Path(str(snakemake.output.metrics)),
    output_provenance=Path(str(snakemake.output.provenance)),
)
log = Path(str(snakemake.log[0]))
log.parent.mkdir(parents=True, exist_ok=True)
log.write_text(
    f"Aggregated {metrics['context_count']} contact-graph contexts\n",
    encoding="utf-8",
)
