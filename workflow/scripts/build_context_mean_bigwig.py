"""Build one context/assay mean background-TMM normalized BigWig."""

from pathlib import Path

from short_read_processing.catalog_visualization import build_context_mean_bigwig


metrics = build_context_mean_bigwig(
    unit_paths={
        str(library): Path(str(path))
        for library, path in dict(snakemake.params.unit_paths).items()
    },
    factor_path=Path(str(snakemake.input.factors)),
    chrom_sizes_path=Path(str(snakemake.input.chrom_sizes)),
    assay=str(snakemake.params.assay),
    context=str(snakemake.params.context),
    output_bigwig=Path(str(snakemake.output.bigwig)),
    output_metrics=Path(str(snakemake.output.metrics)),
    threads=int(snakemake.threads),
)
log = Path(str(snakemake.log[0]))
log.parent.mkdir(parents=True, exist_ok=True)
log.write_text(
    f"Built mean {metrics['assay']} track from {metrics['library_n']} libraries\n",
    encoding="utf-8",
)
