rule callpeak_narrow:
    input:
        bam=lambda wc: FINAL_BAMS[wc.sample],
        control=peak_control_input,
        validated=final_bam_validation_input,
        control_validated=peak_control_validation_input
    output:
        peaks=f"{RESULT_ROOT}/peaks/{{sample}}/{{sample}}_peaks.narrowPeak",
        summits=f"{RESULT_ROOT}/peaks/{{sample}}/{{sample}}_summits.bed",
        xls=f"{RESULT_ROOT}/peaks/{{sample}}/{{sample}}_peaks.xls",
        treat_bdg=f"{RESULT_ROOT}/peaks/{{sample}}/{{sample}}_treat_pileup.bdg",
        control_bdg=f"{RESULT_ROOT}/peaks/{{sample}}/{{sample}}_control_lambda.bdg"
    params:
        command=callpeak_argv
    wildcard_constraints:
        sample=CALLPEAK_NARROW_RE
    conda:
        "../envs/peaks.yaml"
    resources:
        mem_mb=8000
    log:
        f"{RESULT_ROOT}/logs/peaks/{{sample}}.callpeak.log"
    shell:
        "mkdir -p $(dirname {output.peaks:q}) $(dirname {log:q}) && "
        "{params.command} > {log:q} 2>&1"


rule callpeak_broad:
    input:
        bam=lambda wc: FINAL_BAMS[wc.sample],
        control=peak_control_input,
        validated=final_bam_validation_input,
        control_validated=peak_control_validation_input
    output:
        peaks=f"{RESULT_ROOT}/peaks/{{sample}}/{{sample}}_peaks.broadPeak",
        gapped=f"{RESULT_ROOT}/peaks/{{sample}}/{{sample}}_peaks.gappedPeak",
        xls=f"{RESULT_ROOT}/peaks/{{sample}}/{{sample}}_peaks.xls",
        treat_bdg=f"{RESULT_ROOT}/peaks/{{sample}}/{{sample}}_treat_pileup.bdg",
        control_bdg=f"{RESULT_ROOT}/peaks/{{sample}}/{{sample}}_control_lambda.bdg"
    params:
        command=callpeak_argv
    wildcard_constraints:
        sample=CALLPEAK_BROAD_RE
    conda:
        "../envs/peaks.yaml"
    resources:
        mem_mb=8000
    log:
        f"{RESULT_ROOT}/logs/peaks/{{sample}}.callpeak.log"
    shell:
        "mkdir -p $(dirname {output.peaks:q}) $(dirname {log:q}) && "
        "{params.command} > {log:q} 2>&1"
