"""Normalize one context's processed contact sources into a balanced Cooler."""

from pathlib import Path

from short_read_processing.contacts import standardize_context, write_json_atomic


metrics = standardize_context(
    context=str(snakemake.wildcards.context),
    source_manifest=Path(str(snakemake.input.manifest)),
    repository_root=Path(str(snakemake.params.repository_root)),
    target_resolution=int(snakemake.params.resolution),
    workdir=Path(str(snakemake.params.workdir)),
    output=Path(str(snakemake.output.cool)),
)
write_json_atomic(Path(str(snakemake.output.metrics)), metrics)
log = Path(str(snakemake.log[0]))
log.parent.mkdir(parents=True, exist_ok=True)
log.write_text(
    f"Normalized {metrics['replicate_count']} contact sources for "
    f"{metrics['context']} at {metrics['target_resolution_bp']} bp\n",
    encoding="utf-8",
)
