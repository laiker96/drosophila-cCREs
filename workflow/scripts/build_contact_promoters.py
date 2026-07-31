"""Build the promoter-node registry used by contact links."""

from pathlib import Path

from short_read_processing.contact_links import build_promoter_table


metrics = build_promoter_table(
    annotation=Path(str(snakemake.input.annotation)),
    chrom_sizes=Path(str(snakemake.input.chrom_sizes)),
    canonical_chromosomes=list(snakemake.params.canonical_chromosomes),
    promoter_width=int(snakemake.params.promoter_width),
    annotation_checksum=str(snakemake.params.annotation_checksum),
    output=Path(str(snakemake.output.promoters)),
    metrics_output=Path(str(snakemake.output.metrics)),
)
log = Path(str(snakemake.log[0]))
log.parent.mkdir(parents=True, exist_ok=True)
log.write_text(
    f"Built {metrics['promoter_count']} distinct promoter nodes\n",
    encoding="utf-8",
)
