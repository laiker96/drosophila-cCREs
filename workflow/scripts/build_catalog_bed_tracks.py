"""Build portable BED tracks for the context-resolved regulatory catalog."""

from pathlib import Path

from short_read_processing.catalog_visualization import build_catalog_beds


manifest = build_catalog_beds(
    master_bed=Path(str(snakemake.input.master)),
    summit_bed=Path(str(snakemake.input.summits)),
    context_matrix=Path(str(snakemake.input.context_matrix)),
    active_paths={
        str(context): Path(str(path))
        for context, path in dict(snakemake.params.active_paths).items()
    },
    output_master_bed=Path(str(snakemake.output.master)),
    output_context_dhs={
        str(context): Path(str(path))
        for context, path in dict(snakemake.params.context_dhs_paths).items()
    },
    output_active_beds={
        str(context): Path(str(path))
        for context, path in dict(snakemake.params.active_bed_paths).items()
    },
    output_manifest=Path(str(snakemake.output.manifest)),
)
log = Path(str(snakemake.log[0]))
log.parent.mkdir(parents=True, exist_ok=True)
log.write_text(
    f"Built catalog BED tracks for {len(manifest['context_metrics'])} contexts\n",
    encoding="utf-8",
)
