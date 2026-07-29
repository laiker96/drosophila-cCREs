"""Build one portable IGV session containing every regulatory-catalog context."""

from pathlib import Path

from build_igv_session import build_all_contexts_catalog_session


contexts = list(snakemake.params.contexts)
track_count = build_all_contexts_catalog_session(
    contexts=contexts,
    genome=str(snakemake.params.genome),
    atac_bigwigs={
        context: Path(str(snakemake.params.atac_paths[context]))
        for context in contexts
    },
    h3k27ac_bigwigs={
        context: Path(str(snakemake.params.h3k27ac_paths[context]))
        for context in contexts
    },
    context_dhs_beds={
        context: Path(str(snakemake.params.context_dhs_paths[context]))
        for context in contexts
    },
    master_dhs_bed=Path(str(snakemake.input.master_dhs)),
    active_elements_beds={
        context: Path(str(snakemake.params.active_element_paths[context]))
        for context in contexts
    },
    output=Path(str(snakemake.output.session)),
)
log = Path(str(snakemake.log[0]))
log.parent.mkdir(parents=True, exist_ok=True)
log.write_text(
    f"Built all-context IGV session with {len(contexts)} contexts and "
    f"{track_count} tracks\n",
    encoding="utf-8",
)
