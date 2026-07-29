NEXTERA_ADAPTER = "CTGTCTCTTATACACATCT"
TRUSEQ_R1_ADAPTER = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
TRUSEQ_R2_ADAPTER = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

rule fastqc_raw:
    input:
        raw_fastq
    output:
        html=f"{RESULT_ROOT}/qc/fastqc/raw/{{unit}}.html",
        zip=f"{RESULT_ROOT}/qc/fastqc/raw/{{unit}}.zip"
    log:
        f"{RESULT_ROOT}/logs/fastqc/raw/{{unit}}.log"
    threads: 2
    resources:
        mem_mb=1024
    conda:
        "../envs/read_qc.yaml"
    script:
        "../scripts/fastqc.py"


rule trim_pe:
    input:
        r1=lambda wc: lane_input(wc, "r1"),
        r2=lambda wc: lane_input(wc, "r2")
    output:
        r1=(
            f"{WORK_ROOT}/trimmed/{{sample}}/{{lane}}_R1.fastq.gz"
            if OUTPUT_STAGE == "trimming"
            else temp(f"{WORK_ROOT}/trimmed/{{sample}}/{{lane}}_R1.fastq.gz")
        ),
        r2=(
            f"{WORK_ROOT}/trimmed/{{sample}}/{{lane}}_R2.fastq.gz"
            if OUTPUT_STAGE == "trimming"
            else temp(f"{WORK_ROOT}/trimmed/{{sample}}/{{lane}}_R2.fastq.gz")
        ),
        json=f"{RESULT_ROOT}/qc/cutadapt/{{sample}}.{{lane}}.json"
    log:
        f"{RESULT_ROOT}/logs/cutadapt/{{sample}}.{{lane}}.log"
    wildcard_constraints:
        sample=PE_SAMPLE_RE
    params:
        adapter_r1=lambda wc: adapter_value(wc, "r1"),
        adapter_r2=lambda wc: adapter_value(wc, "r2"),
        quality=lambda wc: trim_parameter(wc, "trimming", "quality_cutoff"),
        minimum_length=lambda wc: trim_parameter(wc, "trimming", "minimum_length"),
        error_rate=lambda wc: trim_parameter(wc, "trimming", "error_rate"),
        minimum_overlap=lambda wc: trim_parameter(wc, "trimming", "minimum_overlap")
    threads: 4
    resources:
        mem_mb=2048
    conda:
        "../envs/read_qc.yaml"
    shell:
        "mkdir -p $(dirname {output.r1:q}) $(dirname {output.json:q}) $(dirname {log:q}) && "
        "cutadapt -j {threads} -a {params.adapter_r1:q} -A {params.adapter_r2:q} "
        "-q {params.quality} -m {params.minimum_length} -e {params.error_rate} "
        "-O {params.minimum_overlap} --json {output.json:q} "
        "-o {output.r1:q} -p {output.r2:q} {input.r1:q} {input.r2:q} > {log:q} 2>&1"


rule trim_se:
    input:
        r1=lambda wc: lane_input(wc, "r1")
    output:
        r1=(
            f"{WORK_ROOT}/trimmed/{{sample}}/{{lane}}_SE.fastq.gz"
            if OUTPUT_STAGE == "trimming"
            else temp(f"{WORK_ROOT}/trimmed/{{sample}}/{{lane}}_SE.fastq.gz")
        ),
        json=f"{RESULT_ROOT}/qc/cutadapt/{{sample}}.{{lane}}.json"
    log:
        f"{RESULT_ROOT}/logs/cutadapt/{{sample}}.{{lane}}.log"
    wildcard_constraints:
        sample=SE_SAMPLE_RE
    params:
        adapter=lambda wc: adapter_value(wc, "r1"),
        quality=lambda wc: trim_parameter(wc, "trimming", "quality_cutoff"),
        minimum_length=lambda wc: trim_parameter(wc, "trimming", "minimum_length"),
        error_rate=lambda wc: trim_parameter(wc, "trimming", "error_rate"),
        minimum_overlap=lambda wc: trim_parameter(wc, "trimming", "minimum_overlap")
    threads: 4
    resources:
        mem_mb=2048
    conda:
        "../envs/read_qc.yaml"
    shell:
        "mkdir -p $(dirname {output.r1:q}) $(dirname {output.json:q}) $(dirname {log:q}) && "
        "cutadapt -j {threads} -a {params.adapter:q} -q {params.quality} "
        "-m {params.minimum_length} -e {params.error_rate} -O {params.minimum_overlap} "
        "--json {output.json:q} -o {output.r1:q} {input.r1:q} > {log:q} 2>&1"


rule fastqc_trimmed:
    input:
        trimmed_fastq
    output:
        html=f"{RESULT_ROOT}/qc/fastqc/trimmed/{{unit}}.html",
        zip=f"{RESULT_ROOT}/qc/fastqc/trimmed/{{unit}}.zip"
    log:
        f"{RESULT_ROOT}/logs/fastqc/trimmed/{{unit}}.log"
    threads: 2
    resources:
        mem_mb=1024
    conda:
        "../envs/read_qc.yaml"
    script:
        "../scripts/fastqc.py"
