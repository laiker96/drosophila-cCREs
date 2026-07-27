"""Build descriptive QC outputs from the completed activity tables."""

from pathlib import Path

from short_read_processing.activity_qc import build_activity_qc_outputs


metrics = build_activity_qc_outputs(
    library_signal=Path(str(snakemake.input.library_signal)),
    context_signal=Path(str(snakemake.input.context_signal)),
    qnorm_reference=Path(str(snakemake.input.reference)),
    activity_table=Path(str(snakemake.input.activity)),
    activity_provenance=Path(str(snakemake.input.provenance)),
    atlas_contexts=list(snakemake.params.atlas_contexts),
    reference_context=str(snakemake.params.reference_context),
    output_correlations=Path(str(snakemake.output.correlations)),
    output_distributions=Path(str(snakemake.output.distributions)),
    output_metrics=Path(str(snakemake.output.metrics)),
    output_report=Path(str(snakemake.output.report)),
)
Path(str(snakemake.log[0])).parent.mkdir(parents=True, exist_ok=True)
Path(str(snakemake.log[0])).write_text(
    "Built descriptive activity QC for "
    f"{metrics['master_dhs_count']} master DHSs, "
    f"{metrics['atlas_context_count']} contexts, and "
    f"{metrics['library_count']} libraries\n",
    encoding="utf-8",
)
