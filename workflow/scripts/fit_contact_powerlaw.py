"""Fit context and atlas-wide dm6 contact-decay models."""

from pathlib import Path

from short_read_processing.contacts import fit_atlas_powerlaw


result = fit_atlas_powerlaw(
    contacts={
        str(context): Path(str(path))
        for context, path in dict(snakemake.params.contacts).items()
    },
    canonical_chromosomes=list(snakemake.params.canonical_chromosomes),
    maximum_distance=int(snakemake.params.maximum_distance),
    output=Path(str(snakemake.output.powerlaw)),
)
log = Path(str(snakemake.log[0]))
log.parent.mkdir(parents=True, exist_ok=True)
log.write_text(
    f"Fitted contact decay for {len(result['contexts'])} observed contexts\n",
    encoding="utf-8",
)
