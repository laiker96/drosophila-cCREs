"""Build the canonical background-TMM master-DHS activity table."""

from pathlib import Path

from short_read_processing.activity import sha256_file
from short_read_processing.activity_tmm import (
    TMM_BACKGROUND_METHOD,
    build_tmm_activity_outputs,
)


activity = dict(snakemake.params.activity_config)
signal_paths = {
    str(library_id): Path(str(path))
    for library_id, path in dict(snakemake.params.signal_paths).items()
}
workflow_provenance = dict(snakemake.params.workflow_provenance)
provenance = {
    "workflow_semantic_sha256": workflow_provenance.get("semantic_sha256", ""),
    "master": {
        key: activity["master"][key]
        for key in (
            "genome",
            "method",
            "source_project",
            "source_run_id",
            "master_bed_sha256",
            "summits_bed_sha256",
            "membership_tsv_sha256",
            "context_matrix_tsv_sha256",
            "stats_json_sha256",
        )
    },
    "libraries": [
        {
            **{
                key: library[key]
                for key in (
                    "id",
                    "assay",
                    "context",
                    "layout",
                    "bam_sha256",
                    "bai_sha256",
                    "filtering_contract",
                    "qc_status",
                    "source_project",
                    "source_run_id",
                )
            },
            "estimated_fragment_length_bp": library.get(
                "estimated_fragment_length_bp"
            ),
        }
        for library in activity["libraries"]
    ],
    "excluded_libraries": workflow_provenance.get(
        "excluded_activity_libraries", []
    ),
    "input_signal_sha256": {
        library_id: sha256_file(path)
        for library_id, path in sorted(signal_paths.items())
    },
    "tmm_estimation": {
        "count_matrix": {
            "path": str(Path(str(snakemake.input.counts)).resolve()),
            "sha256": sha256_file(Path(str(snakemake.input.counts))),
        },
        "library_metadata": {
            "path": str(Path(str(snakemake.input.metadata)).resolve()),
            "sha256": sha256_file(Path(str(snakemake.input.metadata))),
        },
        "software_receipt": {
            "path": str(Path(str(snakemake.input.receipt)).resolve()),
            "sha256": sha256_file(Path(str(snakemake.input.receipt))),
        },
    },
}
metrics = build_tmm_activity_outputs(
    method=TMM_BACKGROUND_METHOD,
    signal_paths=signal_paths,
    factor_path=Path(str(snakemake.input.factors)),
    contexts=list(snakemake.params.contexts),
    output_context_signal=Path(str(snakemake.output.context_signal)),
    output_activity=Path(str(snakemake.output.activity)),
    output_metrics=Path(str(snakemake.output.metrics)),
    output_provenance=Path(str(snakemake.output.provenance)),
    provenance=provenance,
)
Path(str(snakemake.log[0])).parent.mkdir(parents=True, exist_ok=True)
Path(str(snakemake.log[0])).write_text(
    f"Built background-TMM activity for {metrics['master_dhs_count']} master DHSs\n",
    encoding="utf-8",
)
