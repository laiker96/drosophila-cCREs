"""Aggregate activity-library signals and write the final normalized tables."""

from pathlib import Path

from short_read_processing.activity import build_activity_outputs, sha256_file


activity = dict(snakemake.params.activity_config)
signal_paths = {
    library_id: Path(str(path))
    for library_id, path in dict(snakemake.params.signal_paths).items()
}
context_paths = {
    context: Path(str(path))
    for context, path in dict(snakemake.params.context_paths).items()
}
provenance = {
    "workflow_semantic_sha256": dict(
        snakemake.params.workflow_provenance
    ).get("semantic_sha256", ""),
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
                    "cohort",
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
    "input_signal_sha256": {
        library_id: sha256_file(path)
        for library_id, path in sorted(signal_paths.items())
    },
    "atac_fragment_maximum": activity["atac_fragment_maximum"],
    "excluded_libraries": dict(snakemake.params.workflow_provenance).get(
        "excluded_activity_libraries", []
    ),
}
metrics = build_activity_outputs(
    signal_paths=signal_paths,
    atlas_contexts=list(snakemake.params.atlas_contexts),
    reference_context=str(snakemake.params.reference_context),
    output_library_signal=Path(str(snakemake.output.library_signal)),
    output_context_signal=Path(str(snakemake.output.context_signal)),
    output_reference=Path(str(snakemake.output.reference)),
    output_activity=Path(str(snakemake.output.activity)),
    output_context_views=context_paths,
    output_metrics=Path(str(snakemake.output.metrics)),
    output_provenance=Path(str(snakemake.output.provenance)),
    provenance=provenance,
)
Path(str(snakemake.log[0])).parent.mkdir(parents=True, exist_ok=True)
Path(str(snakemake.log[0])).write_text(
    f"Built activity table for {metrics['master_dhs_count']} master DHSs\n",
    encoding="utf-8",
)
