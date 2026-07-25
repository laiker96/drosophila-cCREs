"""Master-DHS activity counting and quantile normalization."""


if ACTIVITY:
    rule validate_activity_bam:
        input:
            bam=lambda wc: ACTIVITY_LIBRARIES[wc.library]["bam"],
            bai=lambda wc: ACTIVITY_LIBRARIES[wc.library]["bai"],
            chrom_sizes=str(REFERENCE["chrom_sizes"])
        output:
            validation=f"{RESULT_ROOT}/provenance/activity/bams/{{library}}.validated.json"
        params:
            expected=lambda wc: ACTIVITY_LIBRARIES[wc.library]
        wildcard_constraints:
            library=ACTIVITY_LIBRARY_RE
        threads: 2
        resources:
            mem_mb=2000
        conda:
            "../envs/alignment.yaml"
        log:
            f"{RESULT_ROOT}/logs/activity/validation/{{library}}.log"
        script:
            "../scripts/validate_external_bam.py"


    rule validate_activity_master:
        input:
            master_bed=str(ACTIVITY["master"]["master_bed"]),
            summits_bed=str(ACTIVITY["master"]["summits_bed"]),
            membership_tsv=str(ACTIVITY["master"]["membership_tsv"]),
            context_matrix_tsv=str(ACTIVITY["master"]["context_matrix_tsv"]),
            stats_json=str(ACTIVITY["master"]["stats_json"]),
            chrom_sizes=str(REFERENCE["chrom_sizes"])
        output:
            validation=ACTIVITY_MASTER_VALIDATION
        params:
            expected=ACTIVITY["master"]
        resources:
            mem_mb=2000
        conda:
            "../envs/reporting.yaml"
        log:
            f"{RESULT_ROOT}/logs/activity/validation/master.log"
        script:
            "../scripts/validate_external_master.py"


    rule prepare_activity_atac_insertions:
        input:
            bam=lambda wc: ACTIVITY_LIBRARIES[wc.library]["bam"],
            bai=lambda wc: ACTIVITY_LIBRARIES[wc.library]["bai"],
            validated=lambda wc: ACTIVITY_BAM_VALIDATIONS[wc.library],
            chrom_sizes=str(REFERENCE["chrom_sizes"])
        output:
            bed=f"{ACTIVITY_WORK}/units/{{library}}.bed.gz",
            unit_count=f"{ACTIVITY_WORK}/units/{{library}}.count.txt"
        params:
            maximum=int(ACTIVITY["atac_fragment_maximum"])
        wildcard_constraints:
            library=ACTIVITY_ATAC_LIBRARY_RE
        threads: 4
        resources:
            mem_mb=6000
        conda:
            "../envs/atac_qc.yaml"
        log:
            f"{RESULT_ROOT}/logs/activity/units/{{library}}.atac.log"
        shell:
            r"""
            mkdir -p $(dirname {output.bed:q}) $(dirname {log:q})
            temporary=$(mktemp -d $(dirname {output.bed:q})/.{wildcards.library}.XXXXXX)
            trap 'rm -rf "$temporary"' EXIT
            samtools view -@ {threads} -b -f 2 -F 3852 \
              -e 'tlen != 0 && tlen > -{params.maximum} && tlen < {params.maximum}' \
              -o "$temporary/short.bam" {input.bam:q} 2> {log:q}
            samtools index -@ {threads} "$temporary/short.bam" \
              "$temporary/short.bam.bai" 2>> {log:q}
            alignmentSieve --ATACshift -b "$temporary/short.bam" \
              -o "$temporary/shifted.bam" --numberOfProcessors {threads} \
              >> {log:q} 2>&1
            bedtools bamtobed -i "$temporary/shifted.bam" \
              | awk -v countfile="$temporary/count" 'BEGIN {{OFS="\t"}}
                  $6 == "+" {{print $1,$2,$2+1,$4,$5,$6; n++; next}}
                  $6 == "-" {{print $1,$3-1,$3,$4,$5,$6; n++; next}}
                  END {{print n+0 > countfile}}' \
              | bedtools sort -g {input.chrom_sizes:q} -i - \
              | pigz -p 2 -c > "$temporary/units.bed.gz"
            test "$(samtools view -c "$temporary/shifted.bam")" \
              -eq "$(cat "$temporary/count")"
            pigz -t "$temporary/units.bed.gz"
            mv "$temporary/units.bed.gz" {output.bed:q}
            mv "$temporary/count" {output.unit_count:q}
            """


    rule prepare_activity_h3k27ac_fragments:
        input:
            bam=lambda wc: ACTIVITY_LIBRARIES[wc.library]["bam"],
            bai=lambda wc: ACTIVITY_LIBRARIES[wc.library]["bai"],
            validated=lambda wc: ACTIVITY_BAM_VALIDATIONS[wc.library],
            chrom_sizes=str(REFERENCE["chrom_sizes"])
        output:
            bed=f"{ACTIVITY_WORK}/units/{{library}}.bed.gz",
            unit_count=f"{ACTIVITY_WORK}/units/{{library}}.count.txt"
        wildcard_constraints:
            library=ACTIVITY_H3K27AC_LIBRARY_RE
        threads: 4
        resources:
            mem_mb=4000
        conda:
            "../envs/alignment.yaml"
        log:
            f"{RESULT_ROOT}/logs/activity/units/{{library}}.h3k27ac.log"
        shell:
            r"""
            mkdir -p $(dirname {output.bed:q}) $(dirname {log:q})
            temporary=$(mktemp -d $(dirname {output.bed:q})/.{wildcards.library}.XXXXXX)
            trap 'rm -rf "$temporary"' EXIT
            samtools view -@ {threads} -f 2 -F 3852 -e 'tlen > 0' \
              {input.bam:q} 2> {log:q} \
              | awk -v countfile="$temporary/count" 'BEGIN {{OFS="\t"}}
                  {{start=$4-1; print $3,start,start+$9,$1,0,"."; n++}}
                  END {{print n+0 > countfile}}' \
              | bedtools sort -g {input.chrom_sizes:q} -i - \
              | pigz -p 2 -c > "$temporary/units.bed.gz"
            pigz -t "$temporary/units.bed.gz"
            mv "$temporary/units.bed.gz" {output.bed:q}
            mv "$temporary/count" {output.unit_count:q}
            """


    rule count_activity_library:
        input:
            units=lambda wc: ACTIVITY_UNIT_BEDS[wc.library],
            total=lambda wc: ACTIVITY_UNIT_COUNTS[wc.library],
            master=str(ACTIVITY["master"]["master_bed"]),
            summits=str(ACTIVITY["master"]["summits_bed"]),
            chrom_sizes=str(REFERENCE["chrom_sizes"]),
            master_validated=ACTIVITY_MASTER_VALIDATION,
            script=str(REPO_ROOT / "src" / "count_activity_units.py"),
            implementation=str(
                REPO_ROOT / "src" / "short_read_processing" / "activity.py"
            )
        output:
            signal=f"{ACTIVITY_ROOT}/libraries/{{library}}.signal.tsv.gz",
            summary=f"{ACTIVITY_ROOT}/libraries/{{library}}.summary.json"
        params:
            assay=lambda wc: ACTIVITY_LIBRARIES[wc.library]["assay"],
            cohort=lambda wc: ACTIVITY_LIBRARIES[wc.library]["cohort"],
            context=lambda wc: ACTIVITY_LIBRARIES[wc.library]["context"]
        wildcard_constraints:
            library=ACTIVITY_LIBRARY_RE
        resources:
            mem_mb=4000
        conda:
            "../envs/reporting.yaml"
        log:
            f"{RESULT_ROOT}/logs/activity/counts/{{library}}.log"
        shell:
            "mkdir -p $(dirname {output.signal:q}) $(dirname {log:q}) && "
            "python {input.script:q} --master-bed {input.master:q} "
            "--summit-bed {input.summits:q} --units-bed {input.units:q} "
            "--chrom-sizes {input.chrom_sizes:q} "
            "--total-units {input.total:q} --library-id {wildcards.library:q} "
            "--assay {params.assay:q} --cohort {params.cohort:q} "
            "--context {params.context:q} --output {output.signal:q} "
            "--summary {output.summary:q} > {log:q} 2>&1"


    rule build_master_dhs_activity_table:
        input:
            signals=list(ACTIVITY_LIBRARY_SIGNALS.values()),
            summaries=list(ACTIVITY_LIBRARY_SUMMARIES.values()),
            master=str(ACTIVITY["master"]["master_bed"]),
            summits=str(ACTIVITY["master"]["summits_bed"]),
            master_validated=ACTIVITY_MASTER_VALIDATION,
            bam_validations=list(ACTIVITY_BAM_VALIDATIONS.values()),
            implementation=str(
                REPO_ROOT / "src" / "short_read_processing" / "activity.py"
            )
        output:
            library_signal=ACTIVITY_LIBRARY_SIGNAL,
            context_signal=ACTIVITY_CONTEXT_SIGNAL,
            reference=ACTIVITY_QNORM_REFERENCE,
            activity=ACTIVITY_TABLE,
            context_views=list(ACTIVITY_CONTEXT_VIEWS.values()),
            metrics=ACTIVITY_METRICS,
            provenance=ACTIVITY_PROVENANCE
        params:
            signal_paths=ACTIVITY_LIBRARY_SIGNALS,
            context_paths=ACTIVITY_CONTEXT_VIEWS,
            atlas_contexts=ACTIVITY_CONTEXTS,
            reference_context=str(ACTIVITY["reference_context"]),
            activity_config=ACTIVITY,
            workflow_provenance=config.get("provenance", {})
        resources:
            mem_mb=8000
        conda:
            "../envs/reporting.yaml"
        log:
            f"{RESULT_ROOT}/logs/activity/build-table.log"
        script:
            "../scripts/build_activity_table.py"
