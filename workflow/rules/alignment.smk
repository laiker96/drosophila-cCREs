if CHUNKED_ALIGNMENT_ENABLED:
    checkpoint split_lane_fastq:
        input:
            reads=chunked_lane_reads
        output:
            manifest=temp(
                f"{ALIGNMENT_CHUNK_ROOT}/reads/{{unit}}/manifest.json"
            )
        log:
            f"{RESULT_ROOT}/logs/alignment/{{unit}}.split.log"
        params:
            layout=chunked_lane_layout,
            records_per_chunk=ALIGNMENT_CHUNK_PAIRS,
            output_dir=lambda wc, output: str(Path(output.manifest).parent)
        threads: 1
        resources:
            mem_mb=1000,
            fastq_split_jobs=1
        conda:
            "../envs/alignment.yaml"
        wildcard_constraints:
            unit=CHUNK_ALIGNMENT_UNIT_RE
        shell:
            r"""
            mkdir -p $(dirname {output.manifest:q}) $(dirname {log:q})
            python workflow/scripts/split_fastq_chunks.py \
              --unit {wildcards.unit:q} --layout {params.layout:q} \
              --records-per-chunk {params.records_per_chunk} \
              --output-dir {params.output_dir:q} {input.reads:q} \
              > {log:q} 2>&1
            """


    rule align_lane_chunk:
        input:
            reads=alignment_chunk_reads,
            index=BT2_INDEX
        output:
            bam=temp(
                f"{ALIGNMENT_CHUNK_ROOT}/bams/"
                "{unit}.{chunk}.coordsort.bam"
            ),
            bowtie=temp(
                f"{ALIGNMENT_CHUNK_ROOT}/logs/"
                "{unit}.{chunk}.bowtie2.log"
            )
        params:
            reads=alignment_chunk_bowtie_arguments,
            layout=alignment_chunk_layout_arguments,
            preset=alignment_chunk_preset,
            sample=alignment_chunk_sample,
            index=lambda wc, input: str(input.index[0]).removesuffix(".1.bt2"),
            workers=worker_threads
        threads: 4
        resources:
            mem_mb=4000
        conda:
            "../envs/alignment.yaml"
        wildcard_constraints:
            unit=CHUNK_ALIGNMENT_UNIT_RE,
            chunk=ALIGNMENT_CHUNK_RE
        shell:
            r"""
            mkdir -p $(dirname {output.bam:q}) $(dirname {output.bowtie:q})
            temporary=$(mktemp -d "${{TMPDIR:-/tmp}}/align.{wildcards.unit}.{wildcards.chunk}.XXXXXX")
            staged_bam={output.bam:q}.partial.$$
            staged_log={output.bowtie:q}.partial.$$
            trap 'rm -rf "$temporary"; rm -f "$staged_bam" "$staged_log"' EXIT
            bowtie2 {params.preset} {params.layout} -x {params.index:q} {params.reads} \
              --rg-id {wildcards.unit:q} --rg SM:{params.sample:q} \
              -p {params.workers} 2> "$temporary/bowtie2.log" \
              | samtools view -u -o "$temporary/unsorted.bam" - \
                  2>> "$temporary/bowtie2.log"
            samtools sort -n -@ {params.workers} \
              -o "$temporary/namesort.bam" "$temporary/unsorted.bam" \
              2>> "$temporary/bowtie2.log"
            samtools fixmate -@ {params.workers} -m \
              "$temporary/namesort.bam" "$temporary/fixmate.bam" \
              2>> "$temporary/bowtie2.log"
            samtools sort -@ {params.workers} \
              -o "$temporary/coordsort.bam" "$temporary/fixmate.bam" \
              2>> "$temporary/bowtie2.log"
            samtools quickcheck -v "$temporary/coordsort.bam" \
              2>> "$temporary/bowtie2.log"
            python workflow/scripts/copy_verified.py \
              "$temporary/coordsort.bam" "$staged_bam" \
              >> "$temporary/bowtie2.log" 2>&1
            python workflow/scripts/copy_verified.py \
              "$temporary/bowtie2.log" "$staged_log" >/dev/null
            samtools quickcheck -v "$staged_bam" \
              2>> "$temporary/bowtie2.log"
            mv "$staged_bam" {output.bam:q}
            mv "$staged_log" {output.bowtie:q}
            """


    rule align_lane:
        input:
            bams=alignment_chunk_bams,
            bowtie=alignment_chunk_logs
        output:
            bam=temp(f"{WORK_ROOT}/alignment/lanes/{{unit}}.coordsort.bam")
        log:
            f"{RESULT_ROOT}/logs/alignment/{{unit}}.bowtie2.log"
        params:
            layout=chunked_lane_layout,
            workers=worker_threads,
            chunk_root=lambda wc: (
                f"{ALIGNMENT_CHUNK_ROOT}/reads/{wc.unit}"
            )
        threads: 4
        resources:
            mem_mb=4000
        conda:
            "../envs/alignment.yaml"
        wildcard_constraints:
            unit=CHUNK_ALIGNMENT_UNIT_RE
        shell:
            r"""
            mkdir -p $(dirname {output.bam:q}) $(dirname {log:q})
            temporary=$(mktemp -d "${{TMPDIR:-/tmp}}/merge-align.{wildcards.unit}.XXXXXX")
            staged_bam={output.bam:q}.partial.$$
            staged_log={log:q}.partial.$$
            trap 'rm -rf "$temporary"; rm -f "$staged_bam" "$staged_log"' EXIT
            python workflow/scripts/aggregate_bowtie2_logs.py \
              --layout {params.layout:q} {input.bowtie:q} \
              > "$temporary/bowtie2.log" 2>&1
            samtools merge -f -@ {params.workers} \
              -o "$temporary/coordsort.bam" {input.bams:q} \
              >> "$temporary/bowtie2.log" 2>&1
            samtools quickcheck -v "$temporary/coordsort.bam" \
              2>> "$temporary/bowtie2.log"
            python workflow/scripts/copy_verified.py \
              "$temporary/coordsort.bam" "$staged_bam" \
              >> "$temporary/bowtie2.log" 2>&1
            python workflow/scripts/copy_verified.py \
              "$temporary/bowtie2.log" "$staged_log" >/dev/null
            samtools quickcheck -v "$staged_bam" 2>> "$temporary/bowtie2.log"
            mv "$staged_bam" {output.bam:q}
            mv "$staged_log" {log:q}
            rm -rf -- {params.chunk_root:q}
            """

else:
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
            python workflow/scripts/copy_verified.py \
              "$temporary/coordsort.bam" "$staged" >> {log:q} 2>&1
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
        python workflow/scripts/copy_verified.py \
          "$temporary/marked.bam" "$staged" >> {log:q} 2>&1
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
