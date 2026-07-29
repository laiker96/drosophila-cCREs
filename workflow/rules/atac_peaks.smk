"""Per-replicate ATAC qpois and optional HMMRATAC peak calling."""


rule filter_atac_short_fragments:
    input:
        bam=lambda wc: FINAL_BAMS[wc.sample],
        bai=lambda wc: FINAL_BAIS[wc.sample],
        validated=final_bam_validation_input
    output:
        bam=temp(f"{ATAC_WORK}/replicates/{{sample}}.fragments-lt{ATAC_FRAGMENT_MAXIMUM}.bam"),
        bai=temp(f"{ATAC_WORK}/replicates/{{sample}}.fragments-lt{ATAC_FRAGMENT_MAXIMUM}.bam.bai")
    params:
        maximum=ATAC_FRAGMENT_MAXIMUM
    wildcard_constraints:
        sample=ATAC_QPOIS_PE_RE
    threads: 4
    resources:
        mem_mb=4000
    conda:
        "../envs/alignment.yaml"
    log:
        f"{RESULT_ROOT}/logs/atac/replicates/{{sample}}.short-fragments.log"
    shell:
        r"""
        mkdir -p $(dirname {output.bam:q}) $(dirname {log:q})
        samtools view -@ {threads} -b -f 2 -F 3852 \
          -e 'tlen != 0 && tlen > -{params.maximum} && tlen < {params.maximum}' \
          -o {output.bam:q}.tmp {input.bam:q} 2> {log:q}
        samtools index -@ {threads} {output.bam:q}.tmp {output.bai:q}.tmp 2>> {log:q}
        samtools quickcheck -v {output.bam:q}.tmp 2>> {log:q}
        mv {output.bam:q}.tmp {output.bam:q}
        mv {output.bai:q}.tmp {output.bai:q}
        """


rule prepare_atac_tn5_insertions:
    input:
        bam=atac_insertion_bam,
        bai=atac_insertion_bai,
        validated=final_bam_validation_input
    output:
        bed=temp(f"{ATAC_WORK}/replicates/{{sample}}.tn5-insertions.bed.gz"),
        insertion_count=temp(f"{ATAC_WORK}/replicates/{{sample}}.tn5-insertions.count.txt")
    wildcard_constraints:
        sample=ATAC_QPOIS_RE
    threads: 4
    resources:
        mem_mb=6000
    conda:
        "../envs/atac_qc.yaml"
    log:
        f"{RESULT_ROOT}/logs/atac/replicates/{{sample}}.tn5-insertions.log"
    shell:
        r"""
        mkdir -p $(dirname {output.bed:q}) $(dirname {log:q})
        shifted={output.bed:q}.shifted.tmp.bam
        temporary_bed={output.bed:q}.tmp
        temporary_count={output.insertion_count:q}.tmp
        rm -f "$shifted" "$temporary_bed" "$temporary_count"
        alignmentSieve --ATACshift -b {input.bam:q} -o "$shifted" \
          --numberOfProcessors {threads} > {log:q} 2>&1
        bedtools bamtobed -i "$shifted" \
          | awk -v countfile="$temporary_count" 'BEGIN {{OFS="\t"}}
              $6 == "+" {{print $1,$2,$2+1,$4,$5,$6; n++; next}}
              $6 == "-" {{print $1,$3-1,$3,$4,$5,$6; n++; next}}
              END {{print n+0 > countfile}}' \
          | pigz -p 2 -c > "$temporary_bed"
        test "$(samtools view -c "$shifted")" -eq "$(cat "$temporary_count")"
        pigz -t "$temporary_bed"
        mv "$temporary_bed" {output.bed:q}
        mv "$temporary_count" {output.insertion_count:q}
        rm -f "$shifted"
        """


rule call_atac_replicate_qpois:
    input:
        insertions=lambda wc: ATAC_INSERTIONS[wc.sample],
        chrom_sizes=str(REFERENCE["chrom_sizes"]),
        clip_script=str(REPO_ROOT / "src" / "clip_bedgraph.py")
    output:
        peaks=temp(f"{ATAC_WORK}/replicates/{{sample}}/peaks/{{sample}}.candidates.narrowPeak"),
        pileup=temp(f"{ATAC_WORK}/replicates/{{sample}}.pileup.bdg"),
        lambda_bdg=temp(f"{ATAC_WORK}/replicates/{{sample}}.lambda.bdg"),
        qpois=temp(f"{ATAC_WORK}/replicates/{{sample}}.qpois.bdg")
    params:
        command=atac_qpois_replicate_callpeak_argv
    wildcard_constraints:
        sample=ATAC_QPOIS_RE
    resources:
        mem_mb=8000
    conda:
        "../envs/peaks.yaml"
    log:
        f"{RESULT_ROOT}/logs/atac/replicates/{{sample}}.macs3-qpois.log"
    shell:
        r"""
        mkdir -p $(dirname {output.peaks:q}) $(dirname {output.qpois:q}) $(dirname {log:q})
        temporary=$(mktemp -d $(dirname {output.qpois:q})/.{wildcards.sample}.macs3.XXXXXX)
        trap 'rm -rf "$temporary"' EXIT
        (cd "$temporary" && {params.command}) > {log:q} 2>&1
        for bedgraph in \
          "$temporary/{wildcards.sample}_treat_pileup.bdg" \
          "$temporary/{wildcards.sample}_control_lambda.bdg"; do
            python {input.clip_script:q} --bedgraph "$bedgraph" \
              --chrom-sizes {input.chrom_sizes:q} --output "$bedgraph.clipped" >> {log:q} 2>&1
            mv "$bedgraph.clipped" "$bedgraph"
        done
        macs3 bdgcmp -t "$temporary/{wildcards.sample}_treat_pileup.bdg" \
          -c "$temporary/{wildcards.sample}_control_lambda.bdg" -m qpois \
          -o "$temporary/{wildcards.sample}_qpois.bdg" >> {log:q} 2>&1
        python {input.clip_script:q} \
          --bedgraph "$temporary/{wildcards.sample}_qpois.bdg" \
          --chrom-sizes {input.chrom_sizes:q} \
          --output "$temporary/{wildcards.sample}_qpois.bdg.clipped" >> {log:q} 2>&1
        mv "$temporary/{wildcards.sample}_qpois.bdg.clipped" \
          "$temporary/{wildcards.sample}_qpois.bdg"
        mv "$temporary/{wildcards.sample}_peaks.narrowPeak" {output.peaks:q}
        mv "$temporary/{wildcards.sample}_treat_pileup.bdg" {output.pileup:q}
        mv "$temporary/{wildcards.sample}_control_lambda.bdg" {output.lambda_bdg:q}
        mv "$temporary/{wildcards.sample}_qpois.bdg" {output.qpois:q}
        """


rule refine_atac_replicate_qpois:
    input:
        peaks=lambda wc: ATAC_REPLICATE_CANDIDATES[wc.sample],
        qpois=lambda wc: ATAC_REPLICATE_QPOIS_BDG[wc.sample],
        script=str(REPO_ROOT / "src" / "refine_atac_qpois_peaks.py"),
        implementation=str(REPO_ROOT / "src" / "short_read_processing" / "qpois_refinement.py")
    output:
        refined=f"{ATAC_ROOT}/replicates/{{sample}}/peaks/{{sample}}.qpois-refined.bed",
        excluded=f"{ATAC_ROOT}/replicates/{{sample}}/peaks/{{sample}}.qpois-excluded.bed",
        stats=f"{ATAC_ROOT}/replicates/{{sample}}/peaks/{{sample}}.qpois-refinement.json"
    params:
        minimum_exponent=int(ATAC_QPOIS["minimum_exponent"]),
        maximum_exponent=int(ATAC_QPOIS["maximum_exponent"]),
        minimum_length=int(ATAC_QPOIS["minimum_length"]),
        maximum_length=int(ATAC_QPOIS["maximum_length"]),
        merge_gap=int(ATAC_QPOIS["merge_gap"])
    wildcard_constraints:
        sample=ATAC_QPOIS_RE
    resources:
        mem_mb=4000
    conda:
        "../envs/atac_qc.yaml"
    log:
        f"{RESULT_ROOT}/logs/atac/replicates/{{sample}}.qpois-refinement.log"
    shell:
        "mkdir -p $(dirname {output.refined:q}) $(dirname {log:q}) && "
        "python {input.script:q} --qpois {input.qpois:q} --peaks {input.peaks:q} "
        "--output {output.refined:q} --excluded {output.excluded:q} "
        "--stats {output.stats:q} --name-prefix {wildcards.sample:q} "
        "--minimum-exponent {params.minimum_exponent} "
        "--maximum-exponent {params.maximum_exponent} "
        "--minimum-length {params.minimum_length} "
        "--maximum-length {params.maximum_length} "
        "--merge-gap {params.merge_gap} > {log:q} 2>&1"


rule hmmratac_replicate:
    input:
        bam=lambda wc: FINAL_BAMS[wc.sample],
        bai=lambda wc: FINAL_BAIS[wc.sample],
        blacklist=str(REFERENCE["blacklist_bed"]),
        validated=final_bam_validation_input
    output:
        peaks=f"{ATAC_ROOT}/replicates/{{sample}}/peaks/{{sample}}.hmmratac.narrowPeak"
    params:
        lower=lambda wc: SAMPLES[wc.sample]["peak_caller"]["lower"],
        upper=lambda wc: SAMPLES[wc.sample]["peak_caller"]["upper"],
        prescan=lambda wc: SAMPLES[wc.sample]["peak_caller"]["prescan_cutoff"],
        outdir=lambda wc: f"{ATAC_ROOT}/replicates/{wc.sample}/peaks"
    wildcard_constraints:
        sample=HMMRATAC_RE
    resources:
        mem_mb=12000
    conda:
        "../envs/peaks.yaml"
    log:
        f"{RESULT_ROOT}/logs/atac/replicates/{{sample}}.hmmratac.log"
    shell:
        "mkdir -p {params.outdir:q} $(dirname {log:q}) && "
        "macs3 hmmratac -i {input.bam:q} -f BAMPE -n {wildcards.sample:q} "
        "--outdir {params.outdir:q} -l {params.lower} -u {params.upper} "
        "-c {params.prescan} -e {input.blacklist:q} > {log:q} 2>&1 && "
        "mv {params.outdir:q}/{wildcards.sample}_accessible_regions.narrowPeak {output.peaks:q}"
