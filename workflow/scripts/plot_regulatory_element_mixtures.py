"""Plot the regulatory-catalog H3K27ac mixture distributions."""

from pathlib import Path

from short_read_processing.regulatory_mixture_plot import build_mixture_distribution_plot


metrics = build_mixture_distribution_plot(
    catalog_path=Path(str(snakemake.input.catalog)),
    mixture_path=Path(str(snakemake.input.mixtures)),
    output_svg=Path(str(snakemake.output.svg)),
    output_bins=Path(str(snakemake.output.bins)),
    output_metrics=Path(str(snakemake.output.metrics)),
)
Path(str(snakemake.log[0])).parent.mkdir(parents=True, exist_ok=True)
Path(str(snakemake.log[0])).write_text(
    f"Plotted {len(metrics['contexts'])} contexts; "
    f"unsupported={','.join(metrics['unsupported_contexts']) or 'none'}\n",
    encoding="utf-8",
)
