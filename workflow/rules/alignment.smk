rule align_lane:
    input:
        reads=trimmed_lane_reads,
        index=BT2_INDEX
    output:
        bam=temp(f"{WORK_ROOT}/alignment/lanes/{{sample}}.{{lane}}.coordsort.bam")
    log:
        f"{RESULT_ROOT}/logs/alignment/{{sample}}.{{lane}}.bowtie2.log"
    params:
        reads=bowtie_lane_arguments,
        layout=bowtie_layout_arguments,
        preset=bowtie_preset,
        index=lambda wc, input: str(input.index[0]).removesuffix(".1.bt2"),
        workers=worker_threads
    threads: 4
    resources:
        mem_mb=4000
    conda:
        "../envs/alignment.yaml"
    wildcard_constraints:
        sample=BUILD_SAMPLE_RE,
        lane=LANE_RE
    shell:
        r"""
        mkdir -p $(dirname {output.bam:q}) $(dirname {log:q})
        temporary=$(mktemp -d "${{TMPDIR:-/tmp}}/align.{wildcards.sample}.{wildcards.lane}.XXXXXX")
        staged={output.bam:q}.partial.$$
        trap 'rm -rf "$temporary"; rm -f "$staged"' EXIT
        bowtie2 {params.preset} {params.layout} -x {params.index:q} {params.reads} \
          --rg-id {wildcards.sample:q}.{wildcards.lane:q} --rg SM:{wildcards.sample:q} \
          -p {params.workers} 2> {log:q} \
          | samtools view -u -o "$temporary/unsorted.bam" - 2>> {log:q}
        samtools sort -n -@ {params.workers} \
          -o "$temporary/namesort.bam" "$temporary/unsorted.bam" 2>> {log:q}
        samtools fixmate -@ {params.workers} -m \
          "$temporary/namesort.bam" "$temporary/fixmate.bam" 2>> {log:q}
        samtools sort -@ {params.workers} \
          -o "$temporary/coordsort.bam" "$temporary/fixmate.bam" 2>> {log:q}
        samtools quickcheck -v "$temporary/coordsort.bam" 2>> {log:q}
        cp "$temporary/coordsort.bam" "$staged"
        samtools quickcheck -v "$staged" 2>> {log:q}
        mv "$staged" {output.bam:q}
        """


rule merge_and_mark_duplicates:
    input:
        bams=sample_lane_bams
    output:
        bam=temp(f"{WORK_ROOT}/alignment/{{sample}}.marked.bam")
    log:
        f"{RESULT_ROOT}/logs/alignment/{{sample}}.merge-markdup.log"
    params:
        workers=worker_threads
    threads: 4
    resources:
        mem_mb=6000
    conda:
        "../envs/alignment.yaml"
    wildcard_constraints:
        sample=BUILD_SAMPLE_RE
    shell:
        r"""
        mkdir -p $(dirname {output.bam:q}) $(dirname {log:q})
        temporary=$(mktemp -d "${{TMPDIR:-/tmp}}/markdup.{wildcards.sample}.XXXXXX")
        staged={output.bam:q}.partial.$$
        trap 'rm -rf "$temporary"; rm -f "$staged"' EXIT
        samtools merge -f -@ {params.workers} \
          -o "$temporary/merged.bam" {input.bams:q} > {log:q} 2>&1
        samtools markdup -@ {params.workers} \
          "$temporary/merged.bam" "$temporary/marked.bam" >> {log:q} 2>&1
        samtools quickcheck -v "$temporary/marked.bam" 2>> {log:q}
        cp "$temporary/marked.bam" "$staged"
        samtools quickcheck -v "$staged" 2>> {log:q}
        mv "$staged" {output.bam:q}
        """


rule filter_bam:
    input:
        bam=f"{WORK_ROOT}/alignment/{{sample}}.marked.bam"
    output:
        bam=f"{RESULT_ROOT}/bam/{{sample}}.final.bam",
        bai=f"{RESULT_ROOT}/bam/{{sample}}.final.bam.bai"
    params:
        required=required_flags,
        excluded=excluded_flags,
        mapq=lambda wc: SAMPLES[wc.sample]["parameters"]["alignment"]["mapq_minimum"],
        mitochondrial=str(REFERENCE["mitochondrial_contig"]),
        remove_mito=lambda wc: int(
            SAMPLES[wc.sample]["parameters"]["filtering"]["remove_mitochondrial"]
        )
    threads: 6
    resources:
        mem_mb=6000
    conda:
        "../envs/alignment.yaml"
    log:
        f"{RESULT_ROOT}/logs/alignment/{{sample}}.filter.log"
    wildcard_constraints:
        sample=BUILD_SAMPLE_RE
    shell:
        r"""
        mkdir -p $(dirname {output.bam:q}) $(dirname {log:q})
        samtools view -@ {threads} -h -q {params.mapq} {params.required} -F {params.excluded} {input.bam:q} \
          | awk -v mt={params.mitochondrial:q} -v remove_mt={params.remove_mito} 'BEGIN{{OFS="\t"}} /^@/ {{print; next}} remove_mt == 0 || $3 != mt {{print}}' \
          | samtools view -u - \
          | samtools sort -@ {threads} -o {output.bam:q} - 2> {log:q}
        samtools index -@ {threads} {output.bam:q} {output.bai:q}
        samtools quickcheck -v {output.bam:q}
        """


rule alignment_stats:
    input:
        bam=lambda wc: FINAL_BAMS[wc.sample],
        bai=lambda wc: FINAL_BAIS[wc.sample],
        validated=final_bam_validation_input
    output:
        flagstat=f"{RESULT_ROOT}/qc/alignment/{{sample}}.flagstat.txt",
        stats=f"{RESULT_ROOT}/qc/alignment/{{sample}}.stats.txt",
        idxstats=f"{RESULT_ROOT}/qc/alignment/{{sample}}.idxstats.txt"
    threads: 2
    wildcard_constraints:
        sample=SAMPLE_RE
    conda:
        "../envs/alignment.yaml"
    log:
        f"{RESULT_ROOT}/logs/alignment/{{sample}}.stats.log"
    shell:
        "mkdir -p $(dirname {output.flagstat:q}) $(dirname {log:q}) && "
        "samtools flagstat -@ {threads} {input.bam:q} > {output.flagstat:q} 2> {log:q} && "
        "samtools stats -@ {threads} {input.bam:q} > {output.stats:q} 2>> {log:q} && "
        "samtools idxstats {input.bam:q} > {output.idxstats:q} 2>> {log:q}"
