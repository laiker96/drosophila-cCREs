"""Build one portable context-specific regulatory-catalog IGV session."""

from pathlib import Path

from build_igv_session import build_catalog_session


track_count = build_catalog_session(
    context=str(snakemake.params.context),
    genome=str(snakemake.params.genome),
    atac_bigwig=Path(str(snakemake.input.atac)),
    h3k27ac_bigwig=Path(str(snakemake.input.h3k27ac)),
    context_dhs_bed=Path(str(snakemake.input.context_dhs)),
    master_dhs_bed=Path(str(snakemake.input.master_dhs)),
    active_elements_bed=Path(str(snakemake.input.active_elements)),
    output=Path(str(snakemake.output.session)),
)
log = Path(str(snakemake.log[0]))
log.parent.mkdir(parents=True, exist_ok=True)
log.write_text(
    f"Built {snakemake.params.context} IGV session with {track_count} tracks\n",
    encoding="utf-8",
)
