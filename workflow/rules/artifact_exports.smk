"""Export automatically reusable artifact manifests from completed stages."""


if SAMPLE_IDS:
    rule export_final_bam_manifest:
        input:
            bams=list(FINAL_BAMS.values()),
            bais=list(FINAL_BAIS.values()),
            validations=list(EXTERNAL_BAM_VALIDATIONS.values()),
            metrics=METRICS_JSON
        output:
            manifest=FINAL_BAM_EXPORT_MANIFEST
        params:
            rows=FINAL_BAM_EXPORT_ROWS
        resources:
            mem_mb=1000
        conda:
            "../envs/reporting.yaml"
        log:
            f"{RESULT_ROOT}/logs/provenance/export-final-bams.log"
        script:
            "../scripts/write_final_bam_manifest.py"


if ATAC_MASTER_ENABLED:
    rule export_master_manifest:
        input:
            master_bed=ATAC_MASTER_BED,
            summits_bed=ATAC_MASTER_SUMMITS,
            membership_tsv=ATAC_MASTER_MEMBERSHIP,
            context_matrix_tsv=ATAC_MASTER_CONTEXT_MATRIX,
            stats_json=ATAC_MASTER_STATS,
            validation=(
                [EXTERNAL_MASTER_VALIDATION]
                if EXTERNAL_MASTER_VALIDATION
                else []
            )
        output:
            manifest=MASTER_EXPORT_MANIFEST
        params:
            metadata=MASTER_EXPORT_METADATA
        resources:
            mem_mb=1000
        conda:
            "../envs/reporting.yaml"
        log:
            f"{RESULT_ROOT}/logs/provenance/export-master-dhs.log"
        script:
            "../scripts/write_master_manifest.py"
