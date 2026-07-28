"""Build the human-editable per-library QC review table."""

import json
from pathlib import Path

from short_read_processing.qc_review import build_review_table


library_qc = json.loads(str(snakemake.params.library_qc))
count = build_review_table(
    final_bam_manifest=Path(str(snakemake.input.manifest)),
    metrics_tsv=Path(str(snakemake.input.metrics)),
    library_qc=library_qc,
    multiqc_report=Path(str(snakemake.input.multiqc)),
    output=Path(str(snakemake.output.table)),
)
log = Path(str(snakemake.log[0]))
log.parent.mkdir(parents=True, exist_ok=True)
log.write_text(
    f"Wrote QC review table for {count} library/libraries\n", encoding="utf-8"
)
