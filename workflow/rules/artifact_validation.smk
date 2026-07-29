"""Validate immutable external artifacts before scientific downstream work."""


if EXTERNAL_BAM_SAMPLES:
    rule validate_external_final_bam:
        input:
            bam=lambda wc: FINAL_BAMS[wc.sample],
            bai=lambda wc: FINAL_BAIS[wc.sample],
            chrom_sizes=str(REFERENCE["chrom_sizes"])
        output:
            validation=f"{RESULT_ROOT}/provenance/external_bams/{{sample}}.validated.json"
        params:
            expected=lambda wc: SAMPLES[wc.sample]["final_bam"]
        wildcard_constraints:
            sample=EXTERNAL_BAM_SAMPLE_RE
        threads: 2
        resources:
            mem_mb=2000
        conda:
            "../envs/alignment.yaml"
        log:
            f"{RESULT_ROOT}/logs/provenance/external_bams/{{sample}}.validate.log"
        script:
            "../scripts/validate_external_bam.py"


if EXTERNAL_QC_PEAK_SAMPLES:
    rule validate_external_qc_peak:
        input:
            peak=lambda wc: SAMPLES[wc.sample]["qc_peak"]["path"]
        output:
            validation=(
                f"{RESULT_ROOT}/provenance/external_qc_peaks/"
                "{sample}.validated.json"
            )
        params:
            expected=lambda wc: SAMPLES[wc.sample]["qc_peak"]
        wildcard_constraints:
            sample=EXTERNAL_QC_PEAK_SAMPLE_RE
        resources:
            mem_mb=1000
        conda:
            "../envs/reporting.yaml"
        log:
            f"{RESULT_ROOT}/logs/provenance/external_qc_peaks/{{sample}}.validate.log"
        script:
            "../scripts/validate_external_qc_peak.py"


if EXTERNAL_MASTER:
    rule validate_external_master:
        input:
            master_bed=ATAC_MASTER_BED,
            summits_bed=ATAC_MASTER_SUMMITS,
            membership_tsv=ATAC_MASTER_MEMBERSHIP,
            context_matrix_tsv=ATAC_MASTER_CONTEXT_MATRIX,
            stats_json=ATAC_MASTER_STATS,
            chrom_sizes=str(REFERENCE["chrom_sizes"])
        output:
            validation=EXTERNAL_MASTER_VALIDATION
        params:
            expected=EXTERNAL_MASTER
        resources:
            mem_mb=2000
        conda:
            "../envs/reporting.yaml"
        log:
            f"{RESULT_ROOT}/logs/provenance/external_master.validate.log"
        script:
            "../scripts/validate_external_master.py"
