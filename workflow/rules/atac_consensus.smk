"""Pool ATAC replicates and retain condition peaks with replicate support."""


rule pool_atac_condition_insertions:
    input:
        beds=atac_condition_insertion_inputs,
        counts=atac_condition_insertion_count_inputs
    output:
        bed=f"{ATAC_WORK}/conditions/{{condition}}.tn5-insertions.bed.gz",
        insertion_count=f"{ATAC_WORK}/conditions/{{condition}}.tn5-insertions.count.txt"
    wildcard_constraints:
        condition=ATAC_QPOIS_CONDITION_RE
    threads: 4
    resources:
        mem_mb=2000
    conda:
        "../envs/alignment.yaml"
    log:
        f"{RESULT_ROOT}/logs/atac/conditions/{{condition}}.pool-insertions.log"
    shell:
        r"""
        mkdir -p $(dirname {output.bed:q}) $(dirname {log:q})
        pigz -p {threads} -dc {input.beds:q} | pigz -p {threads} -c > {output.bed:q}.tmp 2> {log:q}
        awk '{{total += $1}} END {{print total + 0}}' {input.counts:q} > {output.insertion_count:q}.tmp
        pigz -t {output.bed:q}.tmp
        mv {output.bed:q}.tmp {output.bed:q}
        mv {output.insertion_count:q}.tmp {output.insertion_count:q}
        """


rule call_atac_condition_qpois:
    input:
        insertions=lambda wc: ATAC_CONDITION_INSERTIONS[wc.condition],
        chrom_sizes=str(REFERENCE["chrom_sizes"]),
        clip_script=str(REPO_ROOT / "src" / "clip_bedgraph.py")
    output:
        peaks=f"{ATAC_ROOT}/conditions/{{condition}}/peaks/{{condition}}.candidates.narrowPeak",
        pileup=temp(f"{ATAC_WORK}/conditions/{{condition}}.pileup.bdg"),
        lambda_bdg=temp(f"{ATAC_WORK}/conditions/{{condition}}.lambda.bdg"),
        qpois=temp(f"{ATAC_WORK}/conditions/{{condition}}.qpois.bdg")
    params:
        command=atac_qpois_condition_callpeak_argv
    wildcard_constraints:
        condition=ATAC_QPOIS_CONDITION_RE
    resources:
        mem_mb=8000
    conda:
        "../envs/peaks.yaml"
    log:
        f"{RESULT_ROOT}/logs/atac/conditions/{{condition}}.macs3-qpois.log"
    shell:
        r"""
        mkdir -p $(dirname {output.peaks:q}) $(dirname {output.qpois:q}) $(dirname {log:q})
        temporary=$(mktemp -d $(dirname {output.qpois:q})/.{wildcards.condition}.macs3.XXXXXX)
        trap 'rm -rf "$temporary"' EXIT
        (cd "$temporary" && {params.command}) > {log:q} 2>&1
        for bedgraph in \
          "$temporary/{wildcards.condition}_treat_pileup.bdg" \
          "$temporary/{wildcards.condition}_control_lambda.bdg"; do
            python {input.clip_script:q} --bedgraph "$bedgraph" \
              --chrom-sizes {input.chrom_sizes:q} --output "$bedgraph.clipped" >> {log:q} 2>&1
            mv "$bedgraph.clipped" "$bedgraph"
        done
        macs3 bdgcmp -t "$temporary/{wildcards.condition}_treat_pileup.bdg" \
          -c "$temporary/{wildcards.condition}_control_lambda.bdg" -m qpois \
          -o "$temporary/{wildcards.condition}_qpois.bdg" >> {log:q} 2>&1
        python {input.clip_script:q} \
          --bedgraph "$temporary/{wildcards.condition}_qpois.bdg" \
          --chrom-sizes {input.chrom_sizes:q} \
          --output "$temporary/{wildcards.condition}_qpois.bdg.clipped" >> {log:q} 2>&1
        mv "$temporary/{wildcards.condition}_qpois.bdg.clipped" \
          "$temporary/{wildcards.condition}_qpois.bdg"
        mv "$temporary/{wildcards.condition}_peaks.narrowPeak" {output.peaks:q}
        mv "$temporary/{wildcards.condition}_treat_pileup.bdg" {output.pileup:q}
        mv "$temporary/{wildcards.condition}_control_lambda.bdg" {output.lambda_bdg:q}
        mv "$temporary/{wildcards.condition}_qpois.bdg" {output.qpois:q}
        """


rule refine_atac_condition_qpois:
    input:
        peaks=lambda wc: ATAC_CONDITION_CANDIDATES[wc.condition],
        qpois=lambda wc: ATAC_CONDITION_QPOIS_BDG[wc.condition],
        script=str(REPO_ROOT / "src" / "refine_atac_qpois_peaks.py"),
        implementation=str(REPO_ROOT / "src" / "short_read_processing" / "qpois_refinement.py")
    output:
        refined=f"{ATAC_ROOT}/conditions/{{condition}}/peaks/{{condition}}.qpois-refined.bed",
        excluded=f"{ATAC_ROOT}/conditions/{{condition}}/peaks/{{condition}}.qpois-excluded.bed",
        stats=f"{ATAC_ROOT}/conditions/{{condition}}/peaks/{{condition}}.qpois-refinement.json"
    params:
        minimum_exponent=int(ATAC_QPOIS["minimum_exponent"]),
        maximum_exponent=int(ATAC_QPOIS["maximum_exponent"]),
        minimum_length=int(ATAC_QPOIS["minimum_length"]),
        maximum_length=int(ATAC_QPOIS["maximum_length"]),
        merge_gap=int(ATAC_QPOIS["merge_gap"])
    wildcard_constraints:
        condition=ATAC_QPOIS_CONDITION_RE
    resources:
        mem_mb=4000
    conda:
        "../envs/atac_qc.yaml"
    log:
        f"{RESULT_ROOT}/logs/atac/conditions/{{condition}}.qpois-refinement.log"
    shell:
        "mkdir -p $(dirname {output.refined:q}) $(dirname {log:q}) && "
        "python {input.script:q} --qpois {input.qpois:q} --peaks {input.peaks:q} "
        "--output {output.refined:q} --excluded {output.excluded:q} "
        "--stats {output.stats:q} --name-prefix {wildcards.condition:q} "
        "--minimum-exponent {params.minimum_exponent} "
        "--maximum-exponent {params.maximum_exponent} "
        "--minimum-length {params.minimum_length} "
        "--maximum-length {params.maximum_length} "
        "--merge-gap {params.merge_gap} > {log:q} 2>&1"


rule atac_condition_pileup_bigwig:
    input:
        bedgraph=lambda wc: ATAC_CONDITION_PILEUP_BDG[wc.condition],
        chrom_sizes=str(REFERENCE["chrom_sizes"]),
        script=str(REPO_ROOT / "src" / "bedgraph_to_bigwig.py")
    output:
        bigwig=f"{ATAC_ROOT}/conditions/{{condition}}/tracks/{{condition}}.MACS3-pileup.unscaled.bw"
    wildcard_constraints:
        condition=ATAC_QPOIS_CONDITION_RE
    resources:
        mem_mb=2000
    conda:
        "../envs/atac_qc.yaml"
    log:
        f"{RESULT_ROOT}/logs/atac/conditions/{{condition}}.pileup-bigwig.log"
    shell:
        "mkdir -p $(dirname {output.bigwig:q}) $(dirname {log:q}) && "
        "python {input.script:q} --bedgraph {input.bedgraph:q} "
        "--chrom-sizes {input.chrom_sizes:q} --output {output.bigwig:q} > {log:q} 2>&1"


rule atac_condition_qpois_bigwig:
    input:
        bedgraph=lambda wc: ATAC_CONDITION_QPOIS_BDG[wc.condition],
        chrom_sizes=str(REFERENCE["chrom_sizes"]),
        script=str(REPO_ROOT / "src" / "bedgraph_to_bigwig.py")
    output:
        bigwig=f"{ATAC_ROOT}/conditions/{{condition}}/tracks/{{condition}}.qpois.bw"
    wildcard_constraints:
        condition=ATAC_QPOIS_CONDITION_RE
    resources:
        mem_mb=2000
    conda:
        "../envs/atac_qc.yaml"
    log:
        f"{RESULT_ROOT}/logs/atac/conditions/{{condition}}.qpois-bigwig.log"
    shell:
        "mkdir -p $(dirname {output.bigwig:q}) $(dirname {log:q}) && "
        "python {input.script:q} --bedgraph {input.bedgraph:q} "
        "--chrom-sizes {input.chrom_sizes:q} --output {output.bigwig:q} > {log:q} 2>&1"


rule merge_atac_condition_bams:
    input:
        bams=atac_condition_bam_inputs,
        validated=atac_condition_validation_inputs
    output:
        bam=temp(f"{ATAC_WORK}/conditions/{{condition}}.pooled.bam"),
        bai=temp(f"{ATAC_WORK}/conditions/{{condition}}.pooled.bam.bai")
    wildcard_constraints:
        condition=ATAC_HMMRATAC_CONDITION_RE
    threads: 6
    resources:
        mem_mb=8000
    conda:
        "../envs/alignment.yaml"
    log:
        f"{RESULT_ROOT}/logs/atac/conditions/{{condition}}.merge-bams.log"
    shell:
        "mkdir -p $(dirname {output.bam:q}) $(dirname {log:q}) && "
        "samtools merge -f -@ {threads} -o {output.bam:q}.tmp {input.bams:q} > {log:q} 2>&1 && "
        "samtools quickcheck -v {output.bam:q}.tmp 2>> {log:q} && "
        "samtools index -@ {threads} {output.bam:q}.tmp {output.bai:q}.tmp 2>> {log:q} && "
        "mv {output.bam:q}.tmp {output.bam:q} && mv {output.bai:q}.tmp {output.bai:q}"


rule hmmratac_condition:
    input:
        bam=lambda wc: ATAC_CONDITION_HMM_BAMS[wc.condition],
        bai=lambda wc: ATAC_CONDITION_HMM_BAIS[wc.condition],
        blacklist=str(REFERENCE["blacklist_bed"])
    output:
        peaks=f"{ATAC_ROOT}/conditions/{{condition}}/peaks/{{condition}}.hmmratac.narrowPeak"
    params:
        lower=lambda wc: atac_condition_peak_config(wc.condition)["lower"],
        upper=lambda wc: atac_condition_peak_config(wc.condition)["upper"],
        prescan=lambda wc: atac_condition_peak_config(wc.condition)["prescan_cutoff"],
        outdir=lambda wc: f"{ATAC_ROOT}/conditions/{wc.condition}/peaks"
    wildcard_constraints:
        condition=ATAC_HMMRATAC_CONDITION_RE
    resources:
        mem_mb=12000
    conda:
        "../envs/peaks.yaml"
    log:
        f"{RESULT_ROOT}/logs/atac/conditions/{{condition}}.hmmratac.log"
    shell:
        "mkdir -p {params.outdir:q} $(dirname {log:q}) && "
        "macs3 hmmratac -i {input.bam:q} -f BAMPE -n {wildcards.condition:q} "
        "--outdir {params.outdir:q} -l {params.lower} -u {params.upper} "
        "-c {params.prescan} -e {input.blacklist:q} > {log:q} 2>&1 && "
        "mv {params.outdir:q}/{wildcards.condition}_accessible_regions.narrowPeak {output.peaks:q}"


rule atac_hmmratac_condition_bigwig:
    input:
        bam=lambda wc: ATAC_CONDITION_HMM_BAMS[wc.condition],
        bai=lambda wc: ATAC_CONDITION_HMM_BAIS[wc.condition]
    output:
        bigwig=f"{ATAC_ROOT}/conditions/{{condition}}/tracks/{{condition}}.CPM.bw"
    wildcard_constraints:
        condition=ATAC_HMMRATAC_CONDITION_RE
    threads: 6
    resources:
        mem_mb=6000
    conda:
        "../envs/atac_qc.yaml"
    log:
        f"{RESULT_ROOT}/logs/atac/conditions/{{condition}}.CPM-bigwig.log"
    shell:
        "mkdir -p $(dirname {output.bigwig:q}) $(dirname {log:q}) && "
        "bamCoverage -b {input.bam:q} -o {output.bigwig:q} --outFileFormat bigwig "
        "--normalizeUsing CPM --binSize 10 --numberOfProcessors {threads} > {log:q} 2>&1"


rule filter_atac_qpois_replicate_support:
    input:
        pooled=lambda wc: ATAC_CONDITION_REFINED[wc.condition],
        replicates=lambda wc: [
            ATAC_REPLICATE_REFINED[sample]
            for sample in ATAC_CONDITIONS[wc.condition].samples
        ],
        script=str(REPO_ROOT / "src" / "build_atac_consensus.py"),
        implementation=str(REPO_ROOT / "src" / "short_read_processing" / "consensus.py")
    output:
        bed=f"{ATAC_ROOT}/conditions/{{condition}}/peaks/{{condition}}.replicate-supported.bed",
        support=f"{ATAC_ROOT}/conditions/{{condition}}/peaks/{{condition}}.replicate-support.tsv",
        stats=f"{ATAC_ROOT}/conditions/{{condition}}/peaks/{{condition}}.replicate-support.json"
    params:
        replicates=atac_consensus_replicate_arguments,
        minimum_replicates=int(ATAC_CONSENSUS.get("minimum_replicates", 2)),
        overlap_fraction=float(ATAC_CONSENSUS.get("replicate_overlap_fraction", 0.5))
    wildcard_constraints:
        condition=ATAC_QPOIS_CONDITION_RE
    resources:
        mem_mb=4000
    conda:
        "../envs/atac_qc.yaml"
    log:
        f"{RESULT_ROOT}/logs/atac/conditions/{{condition}}.replicate-support.log"
    shell:
        "mkdir -p $(dirname {output.bed:q}) $(dirname {log:q}) && "
        "python {input.script:q} --condition-id {wildcards.condition:q} --peak-method qpois "
        "--pooled-peaks {input.pooled:q} {params.replicates} "
        "--minimum-replicates {params.minimum_replicates} "
        "--overlap-fraction {params.overlap_fraction} "
        "--output-bed {output.bed:q} --support-tsv {output.support:q} "
        "--stats-json {output.stats:q} > {log:q} 2>&1"


rule filter_atac_hmmratac_replicate_support:
    input:
        pooled=lambda wc: ATAC_CONDITION_HMM_PEAKS[wc.condition],
        replicates=lambda wc: [
            ATAC_REPLICATE_HMM_PEAKS[sample]
            for sample in ATAC_CONDITIONS[wc.condition].samples
        ],
        script=str(REPO_ROOT / "src" / "build_atac_consensus.py"),
        implementation=str(REPO_ROOT / "src" / "short_read_processing" / "consensus.py")
    output:
        bed=f"{ATAC_ROOT}/conditions/{{condition}}/peaks/{{condition}}.replicate-supported.bed",
        support=f"{ATAC_ROOT}/conditions/{{condition}}/peaks/{{condition}}.replicate-support.tsv",
        stats=f"{ATAC_ROOT}/conditions/{{condition}}/peaks/{{condition}}.replicate-support.json"
    params:
        replicates=atac_consensus_replicate_arguments,
        minimum_replicates=int(ATAC_CONSENSUS.get("minimum_replicates", 2)),
        overlap_fraction=float(ATAC_CONSENSUS.get("replicate_overlap_fraction", 0.5))
    wildcard_constraints:
        condition=ATAC_HMMRATAC_CONDITION_RE
    resources:
        mem_mb=4000
    conda:
        "../envs/atac_qc.yaml"
    log:
        f"{RESULT_ROOT}/logs/atac/conditions/{{condition}}.replicate-support.log"
    shell:
        "mkdir -p $(dirname {output.bed:q}) $(dirname {log:q}) && "
        "python {input.script:q} --condition-id {wildcards.condition:q} --peak-method hmmratac "
        "--pooled-peaks {input.pooled:q} {params.replicates} "
        "--minimum-replicates {params.minimum_replicates} "
        "--overlap-fraction {params.overlap_fraction} "
        "--output-bed {output.bed:q} --support-tsv {output.support:q} "
        "--stats-json {output.stats:q} > {log:q} 2>&1"


if ATAC_MASTER_BUILD_ENABLED:
    rule build_atac_master_dhs:
        input:
            peaks=list(ATAC_CONDITION_CONSENSUS.values()),
            signals=list(ATAC_CONDITION_SIGNAL_BW.values()),
            chrom_sizes=str(REFERENCE["chrom_sizes"]),
            script=str(REPO_ROOT / "src" / "build_atac_master.py"),
            implementation=str(REPO_ROOT / "src" / "short_read_processing" / "master_dhs.py")
        output:
            bed=ATAC_MASTER_BED,
            summits=ATAC_MASTER_SUMMITS,
            membership=ATAC_MASTER_MEMBERSHIP,
            context_matrix=ATAC_MASTER_CONTEXT_MATRIX,
            stats=ATAC_MASTER_STATS
        params:
            inputs=atac_master_arguments,
            summit_max_distance=int(ATAC_MASTER.get("summit_max_distance", 150)),
            minimum_summit_separation=int(
                ATAC_MASTER.get("minimum_summit_separation", 50)
            )
        threads: 4
        resources:
            mem_mb=4000
        conda:
            "../envs/atac_qc.yaml"
        log:
            f"{RESULT_ROOT}/logs/atac/master-dhs.log"
        shell:
            "mkdir -p $(dirname {output.bed:q}) $(dirname {log:q}) && "
            "python {input.script:q} {params.inputs} --chrom-sizes {input.chrom_sizes:q} "
            "--summit-max-distance {params.summit_max_distance} --workers {threads} "
            "--minimum-summit-separation {params.minimum_summit_separation} "
            "--output-bed {output.bed:q} --summit-bed {output.summits:q} "
            "--membership-tsv {output.membership:q} "
            "--context-matrix-tsv {output.context_matrix:q} "
            "--stats-json {output.stats:q} > {log:q} 2>&1"
