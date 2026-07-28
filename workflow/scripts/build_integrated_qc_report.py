"""Build the integrated catalog QC report and render its PDF."""

from pathlib import Path

from short_read_processing.integrated_report import build_integrated_qc_report


metrics = build_integrated_qc_report(
    config=dict(snakemake.params.workflow_config),
    source_files=[dict(source) for source in snakemake.params.report_sources],
    current_files={
        str(name): Path(str(path))
        for name, path in dict(snakemake.params.current_files).items()
    },
    output_html=Path(str(snakemake.output.html)),
    output_pdf=Path(str(snakemake.output.pdf)),
    output_metrics=Path(str(snakemake.output.metrics)),
)
Path(str(snakemake.log[0])).parent.mkdir(parents=True, exist_ok=True)
Path(str(snakemake.log[0])).write_text(
    f"Built integrated QC report with {metrics['source_file_count']} upstream files\n",
    encoding="utf-8",
)
