# Implementation plan: master-DHS activity table for ABC

## 1. Objective

Extend the canonical workflow so that its final atlas endpoint is a
quantile-normalized ATAC/H3K27ac activity table over the existing
variable-width master DHS registry.

The workflow will:

1. Build the summit-aware, non-overlapping master DHS registry from the atlas
   ATAC contexts.
2. Quantify ATAC and H3K27ac signal for every biological library over exactly
   those master intervals.
3. Convert library-level counts to counts per million per kilobase
   (`CPM_per_kb`), then aggregate biological libraries with equal weight by
   context.
4. Quantile-normalize each context, separately by assay, to the corresponding
   frozen, QC-qualified reference distribution selected explicitly for the
   activity run.
5. Return a deterministic context-by-master-DHS table containing raw,
   depth/length-normalized, quantile-normalized, and combined activity values.

The output is an activity input for a downstream Drosophila ABC workflow. Gene
annotation, forced promoter inclusion, contact estimation, enhancer-gene
pairing, ABC scoring, and score calibration remain downstream responsibilities.
The activity table must not be described as a drop-in replacement for the
current Broad ABC `EnhancerList.txt` until that downstream interface is
implemented and tested.

## 2. Scientific decisions

### 2.1 Element universe

- Use `atac/master/master_dhs.bed` as the sole interval universe.
- Preserve its variable-width, non-overlapping coordinates and stable
  `master_dhs_id` values.
- Do not resize master DHSs to 500 bp.
- Do not merge nearby master elements for activity quantification.
- Do not add reference-cell peaks to the atlas master set. Reference libraries
  are used for normalization only and are quantified on the atlas-derived
  master coordinates.
- Require strict BED6, positive interval widths, unique IDs, reference
  contigs, chromosome-order sorting, and no overlapping master intervals.
- Treat BED coordinates as zero-based, half-open intervals throughout.
- Reject an empty master registry; an activity table with no elements is not a
  valid successful output.
- Preserve short elements, including boundary-clipped elements below 50 bp.
  Record their widths so their increased sampling variance remains visible.

This is an adapted ABC element definition. Published thresholds calibrated on
approximately 500-bp human ABC candidates must not be assumed to transfer to
these Drosophila elements.

### 2.2 ATAC counting unit

For paired-end ATAC, quantify the same two-ended Tn5 insertion representation
used by the qpois branch:

1. Start from the already filtered, indexed final BAM.
2. Retain proper pairs with `0 < abs(TLEN) < 150`.
3. Apply the existing Tn5 offsets.
4. Emit one one-base insertion record from each mate.
5. Count insertion records overlapping each master DHS.

An insertion at genomic position `p` belongs to interval `[start, end)` only
when `start <= p < end`. Because master intervals do not overlap, an insertion
can contribute to at most one master element.

Definitions for ATAC library `l` and master DHS `e`:

```text
atac_raw_count[e,l] = number of retained Tn5 insertion records in e
atac_total_units[l] = total retained Tn5 insertion records genome-wide
```

The denominator is the insertion count emitted by the existing
`prepare_atac_tn5_insertions` logic, not total BAM alignments and not total
sequenced read pairs.

Single-end ATAC, if supported by an input atlas in the future, will use its
existing one-ended qpois insertion representation and the corresponding total
insertion count. Paired- and single-end libraries must not be mixed silently
within one context.

### 2.3 H3K27ac counting unit

For paired-end H3K27ac, quantify one fragment per proper template:

1. Start from the already filtered, indexed final IP BAM.
2. Select one representative alignment per proper pair using positive TLEN.
3. Construct the genomic fragment interval from alignment start and TLEN.
4. Count fragments that overlap each master DHS.
5. Count every usable fragment once in the genome-wide denominator.

Use zero-based, half-open fragment intervals and require a non-empty
intersection for an element-level count. A fragment may count once in each of
multiple master elements that it genuinely overlaps; it is not arbitrarily
assigned to one element. It still contributes exactly once to the genome-wide
denominator. Record this distinction in the output documentation.

Definitions:

```text
h3k27ac_raw_count[e,l] = number of retained IP fragments overlapping e
h3k27ac_total_units[l] = total retained IP fragments genome-wide
```

The activity workflow consumes H3K27ac IP BAMs only. ChIP input/control
libraries are deliberately omitted from the selected activity-reference table
and are not processed, subtracted, ratioed, or included in the default ABC
activity signal. Any future input-adjusted activity must be a new, explicitly
named method and `run_id`, not a silent change to this one.

For single-end H3K27ac, infer one fragment per retained alignment using the
reviewed library-specific `estimated_fragment_length_bp`:

1. Retain mapped primary alignments while excluding QC-failed, duplicate, and
   supplementary records from the already filtered final BAM.
2. Convert each alignment to its zero-based reference-aligned BED interval.
3. Anchor the inferred fragment at the read's five-prime coordinate: extend
   right from the aligned start for `+` reads and left from the aligned end for
   `-` reads.
4. Clip inferred fragments at chromosome boundaries and sort them in reference
   chromosome order.
5. Count every inferred fragment once in the genome-wide denominator and count
   its non-empty overlaps with master DHS intervals.

The estimate is the primary phantompeakqualtools cross-correlation estimate
selected during QC review. Record the library layout and estimate in validation
receipts and activity provenance. Never substitute another library's estimate
or silently apply a cohort-wide default.

### 2.4 Depth and length normalization

For assay `a`, element `e`, and biological library `l`:

```text
CPM_per_kb[a,e,l]
    = raw_count[a,e,l] / total_units[a,l] * 1,000,000
      / (width_bp[e] / 1,000)

    = raw_count[a,e,l] * 1,000,000,000
      / (total_units[a,l] * width_bp[e])
```

Requirements:

- Use the original master interval width `end - start`.
- Do not round before context aggregation or quantile normalization.
- Write numeric values with a stable documented precision.
- Reject zero-width elements and zero-unit libraries.
- Retain raw counts and denominators so every normalized value can be
  recomputed.
- Call the value `CPM_per_kb` throughout the code and output. Do not label it
  simply `CPM`, because it is also length normalized and is therefore
  RPKM-like.

### 2.5 Biological-replicate aggregation

Technical runs continue to merge into one biological `library_id` before
duplicate marking. Activity aggregation occurs only after biological-library
quantification.

For each assay, element, and context:

```text
context_CPM_per_kb[a,e,c]
    = arithmetic mean of CPM_per_kb[a,e,l]
      over biological libraries l assigned to context c
```

This gives every biological library equal weight rather than allowing the
deepest library to dominate a pooled count.

Also save:

- number of contributing biological libraries;
- sum and arithmetic mean of raw counts;
- arithmetic mean and sample standard deviation of library `CPM_per_kb`;
- the contributing library IDs in deterministic order.

Define the sample standard deviation as zero when only one biological library
contributes, and use the usual `n - 1` denominator when two or more contribute.
Do not emit `NaN` for a valid single-library context.

Do not average ATAC and H3K27ac together. Do not average across contexts.
Contexts missing either assay are invalid for the final ABC activity table and
must cause a clear configuration error rather than receiving zero-filled
values.

### 2.6 Reference qualification and selection

The leading external reference candidate is untreated S2 at T0 from GSE95689.
It is the only reviewed option with two ATAC-seq replicates and two conventional
H3K27ac ChIP-seq replicates from the same named cell line, paper, laboratory,
state, and time point:

```text
ATAC:     GSM2521726/SRR5312675, GSM2521727/SRR5312676
H3K27ac:  GSM2521732/SRR5312681, GSM2521733/SRR5312682
```

The submitted processed files use dm3. Use the paired-end raw reads and process
them through this repository against dm6 before qualification. Do not include
GSM2521738/GSM2521739 input libraries in the activity sample sheet or activity
DAG.

Retain the following as explicitly named comparators, not interchangeable
replicates:

- ML-DmD17-c3 from GSE245076/GSE245077 is a modern, within-study matched
  ATAC/H3K27ac candidate. Its ATAC libraries require re-evaluation with the
  corrected dm6 TSS reference before acceptance.
- Untreated S2 ATAC from GSE119708 (GSM3381126-GSM3381129) plus S2-DRSC
  H3K27ac from GSE245077 (GSM7835795-GSM7835796) is a cross-study and
  not-proven-identical-line sensitivity reference.
- Untreated Kc ATAC from GSE119708 (GSM3381113-GSM3381116) plus Kc167
  H3K27ac from GSE245077 (GSM7835793-GSM7835794) is likewise a cross-study
  sensitivity reference. Kc and Kc167 must not be asserted to be identical.
- A frozen internal atlas reference, such as a predeclared consensus
  distribution, remains a valid alternative if no external candidate passes
  QC. It is a different normalization method and requires a distinct
  `run_id`.

Selection is configuration, not fallback logic. A production activity run
names exactly one reference cohort and records its accepted QC decision. If
that cohort fails validation or an input is absent, stop; never silently switch
to D17, a cross-study pair, or an atlas context.

Process any candidate taken forward for QC through the same dm6
alignment/filtering and activity-counting definitions used for the atlas.
Quantify the selected reference on the exact atlas master BED and calculate one
vector per assay:

```text
reference_CPM_per_kb[a,e]
    = arithmetic mean of selected-reference replicate CPM_per_kb[a,e,l]
```

Before accepting the reference:

- require both replicate BAMs and indexes;
- require all BAMs to use the same dm6 contig names and lengths as the master
  registry;
- require nonzero total units;
- require the reference QC review to be recorded as accepted in the activity
  configuration or manifest;
- report replicate correlation and library-depth summaries;
- evaluate ATAC TSS enrichment, ATAC FRiP, H3K27ac FRiP, mapping, duplication,
  insert-size, and replicate agreement under one documented QC contract;
- use the same reference cohort for both assays in the primary analysis when a
  matched cohort such as GSE95689 passes QC;
- do not silently replace a failed reference assay with an atlas sample.

ATAC and H3K27ac quantile normalization are mathematically separate, so a
cross-study assay pair can be used for a declared sensitivity analysis.
However, it must not be presented as a matched biological reference. Changing
reference libraries, reference method, or accepted-QC state is a scientific
parameter change and requires a new `run_id`.

### 2.7 Quantile normalization

Quantile-normalize ATAC and H3K27ac independently. Never combine their values
before normalization.

For each assay `a` and atlas context `c`:

1. Take the complete vector `context_CPM_per_kb[a,*,c]` in master-DHS order.
2. Rank it from lowest to highest.
3. Replace its ordered values with the ordered values of
   `reference_CPM_per_kb[a,*]`.
4. Map the replacement values back to the original master-DHS order.

Tie handling must be explicit and deterministic:

- all target elements in one tied group receive the arithmetic mean of the
  reference quantiles covered by that group;
- NaN and infinite inputs are forbidden;
- stable element IDs provide deterministic ordering for diagnostics, but must
  not break a tie into different normalized values.

The atlas and reference vectors must have identical element IDs and length.
Any mismatch is a hard error.

Do not add a pseudocount by default. A zero reference quantile remains zero.
The implementation must verify that each normalized assay/context vector has
the expected reference distribution, subject only to the documented averaging
of reference quantiles for target ties.

### 2.8 Combined activity

For every master DHS and context:

```text
activity = sqrt(atac_qnorm * h3k27ac_qnorm)
```

Save the two assay-specific quantile-normalized values as separate columns.
Do not quantile-normalize `activity` itself. Do not add a pseudocount unless a
future calibrated method introduces one under a new run namespace.

## 3. Input contracts

### 3.1 Existing accession sample sheets

Continue to use the canonical atlas accession table for atlas metadata:

- atlas: `resources/atlas_samples_ip_only.tsv` or the explicitly selected
  controlled alternative.

The selected atlas table is a scientific input and must be recorded in
provenance. The two atlas tables are alternative selections and must never be
combined implicitly.

The production reference accession table is frozen only after candidate QC.
It contains treatment/IP libraries only and follows the canonical sample-sheet
schema. Detailed candidate selection and QC metadata do not belong in
`resources/`; store runtime candidate sheets and QC artifacts under ignored
run namespaces. The currently reviewed D17 candidate table is not implicit
authorization to use D17 as the production reference.

### 3.2 Explicit final-BAM manifest

Add a runtime final-BAM manifest for activity quantification. It is the only
supported way for final-BAM reuse to consume BAMs that were copied from an
earlier run or are stored outside the current canonical result namespace.

Proposed tab-separated columns:

```text
library_id
assay
context
role
layout
bam
bai
genome
filtering_contract
source_project
source_run_id
bam_sha256
bai_sha256
qc_status
estimated_fragment_length_bp
notes
```

Rules:

- `library_id`, `assay`, `context`, and `role` must agree with the selected
  sample sheet.
- Paths are resolved relative to the manifest location or repository root by
  one documented convention; do not mix conventions.
- `genome` must be `dm6`.
- `filtering_contract` must identify the expected final-BAM semantics, for
  example `short-read-processing-final-v1`.
- Activity review manifests must cover every selected treatment library with
  `qc_status=accepted` or `qc_status=rejected`; `pending_review` is an error.
- Rejected rows remain in the reviewed manifest, require a reason in `notes`,
  and are skipped explicitly while their decision is retained in provenance.
- Reviewed single-end histone libraries require the primary
  phantompeakqualtools estimate in `estimated_fragment_length_bp`. Paired-end
  rows leave this field blank.
- Control libraries may be present for provenance but are excluded from
  activity.
- Every required treatment library must occur exactly once, including
  explicitly rejected rows.
- Duplicate paths, missing libraries, unexpected treatment libraries, and
  inconsistent contexts are hard errors.
- Full BAM and BAI SHA-256 values are recorded once when the manifest is
  created and verified before a scientific run.
- Validate BAMs with `samtools quickcheck`, verify readable indexes, and compare
  `@SQ` names and lengths with `references/dm6/dm6.chrom.sizes`.

Do not commit runtime filesystem paths or generated BAM manifests under
`resources/`. Store them under the ignored raw/config or results provenance
namespace.

### 3.3 No-realignment contract

The implemented `--from-stage final-bam --final-bam-manifest PATH` mode is the
no-realignment entry point for downstream peak calling and master construction.
The activity endpoint is the stricter
`--from-stage activity` mode: it consumes accepted atlas/reference final-BAM
manifests plus a frozen master manifest and uses the same artifact validators.

When a final-BAM manifest is supplied:

- final BAMs and BAIs are external, immutable Snakemake inputs;
- no alignment, trimming, download, merge, duplicate-marking, or `filter_bam`
  rule may be in the activity DAG;
- a missing or invalid BAM causes the command to fail;
- there is no automatic fallback to FASTQ download or realignment;
- there is no filename-based recursive BAM discovery;
- input BAMs are never copied, rewritten, indexed in place, or timestamped.

Normal accession-to-results mode may still create missing canonical BAMs.
When those BAMs already exist as complete Snakemake outputs, they are reused.
External reuse is deliberately strict and assay-complete: one run never mixes
manifest BAMs with newly aligned libraries from the same represented assay.
The CLI prints the selected reuse stage and manifest, and per-library
validation receipts identify every BAM actually consumed.
Completed sample-processing runs automatically export
`provenance/manifests/final-bams.tsv`. Pipeline-produced BAMs begin at
`qc_status=pending_review`; QC acceptance is an explicit later decision.

### 3.4 Frozen master-DHS bundle

`--from-stage master --master-manifest PATH` registers a complete immutable
five-file master bundle: BED6 elements, one-base summits, membership table,
context matrix, and statistics JSON. All files carry full SHA-256 values and
must agree structurally. This mode never schedules BAM processing or
`build_atac_master_dhs`; at the current processing endpoint it performs
validation and provenance registration only. The activity branch consumes
these validated paths directly.
Every successful build or validation also exports
`provenance/manifests/master-dhs.tsv`, so the next run can start from the
complete bundle without a manual manifest-construction step.

## 4. Configuration and orchestration

### 4.1 Canonical entry point

Keep `src/run_pipeline.py` as the canonical accession-to-results entry point
and `workflow/Snakefile` as the only production Snakefile.

Activity is implemented in the canonical CLI rather than a second production
workflow:

```text
--from-stage {accessions,final-bam,master,activity}
--final-bam-manifest PATH
--master-manifest PATH
--activity-atlas-bam-manifest PATH
--activity-reference-sheet PATH
--activity-reference-bam-manifest PATH
--activity-reference-context ID
```

- `--from-stage activity` consumes only complete artifact manifests and never
  enters accession processing;
- the positional sample sheet supplies atlas library/context metadata;
- activity requires the selected master, accepted atlas BAMs, accepted
  reference BAMs, and reference metadata;
- either BAM-manifest option is repeatable for assay-specific manifests;
- incompatible or incomplete combinations fail during argument validation.

### 4.2 Combined activity configuration

Generate a resolved activity YAML after validating all inputs. It should
contain:

- project and activity `run_id`;
- genome and reference paths;
- master BED and summit BED paths;
- master BED SHA-256;
- selected atlas and reference sample-sheet paths and SHA-256 values;
- final-BAM manifest path and SHA-256;
- ordered atlas contexts;
- ordered biological libraries per context and assay;
- ordered reference libraries per assay;
- ATAC TLEN and Tn5 definitions;
- H3K27ac fragment definition;
- normalization formula and scale constants;
- replicate aggregation method;
- quantile-normalization and tie method;
- activity formula;
- output root;
- implementation/schema version.

Do not regenerate or replace an identical resolved configuration merely to
change a timestamp. Wall-clock timestamps belong in execution logs, not in the
semantic configuration hash.

### 4.3 Workflow integration

Add shared production rules under `workflow/rules/activity.smk` and include
them from `workflow/Snakefile`. Do not create an experimental parallel
Snakefile.

The activity configuration should select an activity branch in the canonical
Snakefile. The branch treats the master BED and every final BAM/BAI as input
files and targets only activity artifacts.

Proposed dependency flow:

```text
validated master BED
        |
        +------------------------------+
        |                              |
atlas/reference ATAC BAMs       atlas/reference H3K27ac BAMs
        |                              |
short-fragment Tn5 insertions   paired-fragment intervals
        |                              |
per-library raw counts          per-library raw counts
        |                              |
per-library CPM_per_kb          per-library CPM_per_kb
        +---------------+--------------+
                        |
              context aggregation
                        |
           selected reference vectors
                        |
       assay-specific quantile normalization
                        |
             combined activity table
                        |
            validation and provenance
```

Implemented rules:

1. `validate_activity_bam`
2. `validate_activity_master`
3. `prepare_activity_atac_insertions`
4. `prepare_activity_h3k27ac_fragments`
5. `count_activity_library`
6. `build_master_dhs_activity_table`

The final aggregation rule performs deterministic context aggregation,
reference construction, tie-aware quantile normalization, table assembly,
metrics, and provenance as one atomic multi-output operation. The originally
proposed finer-grained names were:

1. `validate_activity_inputs`
2. `prepare_activity_atac_insertions`
3. `prepare_activity_h3k27ac_fragments`
4. `count_activity_library`
5. `summarize_activity_library`
6. `aggregate_activity_context`
7. `build_activity_qnorm_reference`
8. `quantile_normalize_activity`
9. `assemble_master_dhs_activity_table`
10. `validate_master_dhs_activity_table`
11. `write_activity_provenance`

Reuse the existing Tn5 insertion implementation rather than reimplementing
ATAC shifting in a second script. Extract shared logic if necessary, with
regression tests showing byte-equivalent insertion coordinates and counts.

## 5. Output contract

All outputs are namespaced by project and activity run ID:

```text
results/<project>/<activity-run-id>/activity/
```

### 5.1 Canonical final table

Canonical output:

```text
master_dhs_activity.tsv.gz
```

It is long format: one row per `(master_dhs_id, context)`, ordered by reference
chromosome order, start, end, master ID, then configured context order.

Required columns:

```text
master_dhs_id
chrom
start
end
summit
width_bp
context
atac_library_n
atac_library_ids
atac_raw_count_sum
atac_raw_count_mean
atac_total_units_sum
atac_cpm_per_kb
atac_cpm_per_kb_sd
atac_qnorm
h3k27ac_library_n
h3k27ac_library_ids
h3k27ac_raw_count_sum
h3k27ac_raw_count_mean
h3k27ac_total_units_sum
h3k27ac_cpm_per_kb
h3k27ac_cpm_per_kb_sd
h3k27ac_qnorm
activity
```

The canonical quantile-normalized activity values are `atac_qnorm`,
`h3k27ac_qnorm`, and `activity`. Raw and `CPM_per_kb` values remain in the same
table for auditability. In this context-level table, each `*_cpm_per_kb`
column is the equal-weight arithmetic mean across biological libraries; the
replicate-level values remain in `library_signal.tsv.gz`.

### 5.2 Supporting tables

Also save:

```text
library_signal.tsv.gz
```

One row per `(master_dhs_id, library_id)` with assay, cohort
(`atlas`/`reference`), context, raw count, total units, width, and
`CPM_per_kb`. This preserves replicate-level values and is the source for all
context aggregation.

```text
context_signal.pre_qnorm.tsv.gz
```

One row per `(master_dhs_id, context, assay)` with the raw-count summaries and
equal-weight replicate aggregate used as quantile-normalization input.

```text
qnorm_reference.tsv.gz
```

One row per `(master_dhs_id, assay)` with the selected cohort ID, all
reference-replicate values, and their arithmetic mean.

```text
contexts/<context>.activity.tsv.gz
```

Deterministic per-context views of the canonical final table for convenient
downstream ABC ingestion. These are derived outputs, not independent
calculations.

```text
activity_metrics.json
activity_provenance.json
qc/replicate_correlations.tsv
qc/distribution_quantiles.tsv
qc/activity_qc_metrics.json
qc/activity_qc_report.html
```

Metrics include element counts, width distribution, zero fractions, library
depths, replicate correlations, normalization summaries, tie counts, and
min/max/quantiles for every assay/context. Provenance contains input hashes,
resolved parameters, software versions, commands, and output hashes.

Compression must be deterministic. Use a method with fixed metadata (for
example gzip with `mtime=0`) so identical inputs produce byte-identical tables.

## 6. Idempotency and file safety

Every activity rule must satisfy all of the following:

- Declare every scientific input, parameter, code dependency, and output to
  Snakemake.
- Write to a temporary file in the destination filesystem.
- Validate the temporary file before atomic `os.replace`.
- Never modify the master BED, final BAMs, BAIs, sample sheets, or BAM
  manifest.
- Never write an index beside an external BAM. Require the declared BAI.
- Keep library-level intermediates addressable by library ID so independent
  jobs can run concurrently and be reused.
- Use `rerun-incomplete: true`.
- A second identical run must schedule no compute jobs and must not alter
  output mtimes or content hashes.
- A changed master BED, BAM manifest, reference selection, counting method,
  normalization method, or activity formula requires a new `run_id`.
- Refuse to write a scientifically different attempt into an existing result
  namespace.
- Fail on partial, empty, malformed, unsorted, duplicate-ID, or
  chromosome-incompatible inputs.
- Do not delete or overwrite imported BAMs if validation fails.

Before launching activity from final-BAM mode, print and record the resolved BAM for
every required library. This preflight list is the primary defense against
accidental realignment or quantifying the wrong file.

## 7. Reproducibility requirements

- Keep orchestration dependencies in root `.venv`.
- Put any new bioinformatics dependency in the smallest applicable
  `workflow/envs/` file.
- Prefer existing samtools, bedtools, and Python dependencies.
- Set `LC_ALL=C` for sorting and text processing.
- Use `references/dm6/dm6.chrom.sizes` for chromosome order.
- Use stable numeric formatting and deterministic tie handling.
- Record schema/implementation versions in every table header sidecar and
  provenance JSON.
- Record complete input and output SHA-256 hashes.
- Record library denominators and exact normalization constants.
- Record the selected atlas table, reference table and cohort ID, BAM manifest,
  master BED, and resolved activity YAML.
- Preserve per-library values; never make the context mean the only surviving
  representation.
- Regenerate `docs/workflow-dag.svg` after adding activity rule dependencies.
- Update `README.md`, `resources/README.md`, validators, tests, and provenance
  together because this changes the documented scientific endpoint.

## 8. Validation and tests

### 8.1 Unit tests

Add focused tests for:

- `CPM_per_kb` formula with multiple widths, including a 31-bp element;
- zero-width rejection;
- zero-denominator rejection;
- ATAC insertion counting at interval boundaries;
- correct two-ended Tn5 denominator;
- paired H3K27ac fragment construction and exactly-one-fragment accounting;
- fragments overlapping no element and boundary-adjacent elements;
- equal-weight replicate means with unequal library depths;
- sample standard deviation for one and multiple replicates;
- independent ATAC and H3K27ac reference construction;
- quantile normalization with no ties;
- target ties, reference ties, all-zero vectors, and mixed zero/nonzero values;
- deterministic output under permuted input row order;
- geometric-mean activity with zeros;
- empty master BED behavior;
- missing context/assay rejection;
- mismatched atlas/reference master IDs;
- duplicate BAM-manifest libraries and paths;
- wrong contig names or lengths;
- stable numeric formatting and deterministic gzip output.

### 8.2 Small integration fixture

Create only minimal synthetic, coordinate-sorted paired-end BAM fixtures.
Include:

- two atlas contexts;
- two biological libraries in at least one context;
- two selected-reference ATAC and two selected-reference H3K27ac libraries;
- variable-width, non-overlapping master elements;
- a short boundary-clipped element;
- known ATAC insertion and H3 fragment counts;
- zero-signal elements and tied values.

Assert exact raw counts, denominators, `CPM_per_kb`, context means, qnorm
values, activity values, row order, and output hashes.

### 8.3 No-realignment regression

For activity dry runs starting with `--from-stage activity`:

- assert that no download, trim, align, merge, mark-duplicate, or `filter_bam`
  rule appears in the DAG;
- assert that all BAMs and BAIs appear only as input files;
- remove one BAM from the manifest fixture and assert a preflight failure,
  not an alignment plan;
- run the completed fixture a second time and assert “Nothing to be done”;
- record output mtimes and hashes before and after the second run and assert
  they are unchanged.

### 8.4 Workflow verification

Run:

```bash
export MAMBA_ROOT_PREFIX="$PWD/.micromamba"
export XDG_CACHE_HOME="$PWD/.cache"

mamba run --prefix "$PWD/.venv" pytest -q

mamba run --prefix "$PWD/.venv" \
  snakemake --snakefile workflow/Snakefile \
  --configfile tests/fixtures/activity_config.yaml --lint

mamba run --prefix "$PWD/.venv" \
  snakemake --snakefile workflow/Snakefile \
  --configfile tests/fixtures/activity_config.yaml \
  --cores 8 --dry-run

git diff --check
```

The dry run must terminate at the activity table and provenance artifacts
without scheduling upstream read processing.

## 9. Implementation phases

### Phase 1: schemas and input validation

1. Define the final-BAM and master-bundle manifest schemas. **Implemented.**
2. Add typed parsing and cross-field validation. **Implemented.**
3. Add BAM/index/reference and master-bundle validation. **Implemented.**
4. Add the resolved activity configuration model. **Implemented.**
5. Add positive and negative tests. **Implemented for configuration and QC
   acceptance; production-artifact validation remains part of the first run.**

Acceptance criterion: the current atlas and selected-reference BAMs can be
represented and validated without changing, copying, or aligning any BAM.

### Phase 2: per-library counting

1. Reuse/extract the existing ATAC insertion logic. **Implemented with the
   same samtools/alignmentSieve/bedtools definition.**
2. Implement paired H3K27ac fragment generation. **Implemented.**
3. Count assay units over the master BED. **Implemented.**
4. Emit deterministic per-library tables and denominator summaries.
   **Implemented.**
5. Add exact synthetic-count tests. **Implemented for interval counting;
   paired-BAM execution remains to be covered by the integration fixture.**

Acceptance criterion: every library-level normalized value can be reproduced
from its raw count, denominator, and master width.

### Phase 3: context aggregation and selected reference

1. Aggregate biological-library values with equal weight by context.
   **Implemented.**
2. Build separate selected-reference vectors for ATAC and H3K27ac.
   **Implemented.**
3. Emit replicate correlation and zero-fraction metrics. **Implemented.**
4. Enforce reference QC acceptance and complete assay coverage.
   **Implemented.**

Acceptance criterion: aggregation is invariant to input row order and is not
weighted by library depth.

### Phase 4: quantile normalization and final table

1. Implement deterministic tie-aware quantile normalization. **Implemented.**
2. Normalize assays independently. **Implemented.**
3. Calculate geometric-mean activity. **Implemented.**
4. Assemble the canonical long table and per-context views. **Implemented.**
5. Validate distributions, completeness, row order, and finite values.
   **Implemented.**
6. Generate deterministic downstream replicate and normalization QC with raw,
   log1p, and rank correlations; pre/post/reference quantile profiles; zero
   fractions; and tie summaries. **Implemented; scientific acceptance remains
   a manual review decision.**

Acceptance criterion: normalized vectors match their assay-specific selected
reference distributions under the documented tie semantics.

### Phase 5: canonical workflow and CLI integration

1. Add `workflow/rules/activity.smk`. **Implemented.**
2. Add the activity branch to `workflow/Snakefile`. **Implemented.**
3. Add the strict `--from-stage activity` entry in `src/run_pipeline.py`.
   **Implemented.**
4. Generate stable activity configs and provenance. **Implemented.**
5. Add activity targets to the appropriate top-level endpoint.
   **Implemented.**
6. Add the no-realignment dry-run regression. **Implemented.**
7. Add explicit `trimming`, `alignment`, `qc`, `master`, `activity`, and
   `activity-qc` stopping points while keeping the stopping point outside the
   scientific semantic hash. **Implemented.**
8. Bound MACS3 ATAC treatment, lambda, and qpois bedGraphs to the reference
   chromosome sizes before downstream use. **Implemented.**

Acceptance criterion: an activity run from complete final-BAM inputs never
schedules alignment and an identical rerun performs no work.

### Phase 6: documentation and production validation

1. Update the repository endpoint description in `README.md`.
2. Document the selected normalization-reference identifier and QC decision in
   the run provenance; keep detailed candidate QC in the atlas-analysis
   repository.
3. Document formulas, units, missing-data policy, tie behavior, and output
   columns.
4. Regenerate the workflow DAG.
5. Run the full practical verification suite.
6. Execute a production dry run against the actual atlas/reference BAM
   manifest.
7. Review disk and memory estimates before starting production quantification.

Acceptance criterion: a user can reproduce the table from documented inputs
and commands without relying on an undocumented BAM location or manual
post-processing step.

### Phase 7 (deferred): numbered workflow navigation

Organize the workflow rule modules by stable conceptual phase so that a new
maintainer can follow the code from inputs to the activity endpoint. This is a
navigation-only change: Snakemake's DAG remains the authority for execution
order, and reuse modes may skip earlier phases.

Use numeric gaps so later phases can be inserted without renumbering the whole
workflow:

| Prefix | Conceptual phase |
|---|---|
| `00` | shared workflow definitions |
| `10` | reference preparation and external-artifact validation |
| `20` | read processing and read-level QC |
| `30` | alignment, duplicate marking, filtering, and final BAMs |
| `40` | assay signal and replicate peak calling |
| `50` | ATAC context support and master-DHS construction |
| `60` | QC aggregation, reporting, and artifact export |
| `70` | activity counting, CPM-per-kb normalization, and quantile normalization |
| `80` | descriptive activity QC and normalization review report |

Implementation constraints:

1. Prefix only the rule-module filenames under `workflow/rules/` and update the
   canonical `workflow/Snakefile` includes.
2. Keep individual Snakemake rule names, Python module and command names,
   schemas, public result paths, manifest formats, and scientific parameters
   unchanged.
3. Document in `README.md` that the numbers aid navigation and do not impose a
   linear execution order.
4. Make this a separate mechanical change after the current functional work is
   verified; do not mix it with scientific or behavioral modifications.
5. Regenerate the documented DAG and verify that the dry-run rule graph and
   targets are unchanged apart from rule-module paths.
6. Run the full practical verification suite and `git diff --check`.

Acceptance criterion: the numbered module layout makes the logical phases
obvious while producing the same DAG, targets, outputs, and restart behavior.

## 10. Final acceptance criteria

The implementation is complete only when:

- the master DHS registry is built exclusively from atlas ATAC contexts;
- all atlas and selected-reference signals are quantified on exactly the same
  non-overlapping variable-width master elements;
- raw per-library counts and total assay units are retained;
- per-library and per-context `CPM_per_kb` values are retained;
- ATAC and H3K27ac are quantile-normalized separately to the explicitly
  selected reference;
- the canonical table contains raw, `CPM_per_kb`, qnorm, and activity values;
- the downstream activity-QC report records replicate agreement, sparsity,
  tie behavior, and assay-specific reference matching without silently
  imposing an acceptance threshold;
- every atlas context has both assays or fails explicitly;
- ChIP input/control libraries are absent from the selected activity-reference
  table and activity DAG;
- the activity DAG entered from final BAMs contains no alignment or download
  jobs;
- all existing BAMs are immutable declared inputs;
- rerunning identical inputs performs no work and does not change outputs;
- changed scientific parameters require a new run namespace;
- provenance is sufficient to identify and verify every input BAM, interval
  set, reference library, parameter, implementation version, and output;
- unit, integration, negative, dry-run, and idempotency tests pass;
- documentation and the workflow DAG describe the new endpoint accurately.
