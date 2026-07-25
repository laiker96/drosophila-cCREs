# short-read-processing

Snakemake pipeline for accession-to-peaks ATAC-seq and ChIP-seq processing.
It downloads public runs, prepares `dm6` or `hg38`, processes technical runs in
parallel, and produces reproducible peak, signal, and QC outputs.

![Simplified workflow DAG](docs/workflow-dag.svg)

## Current endpoints

- **ATAC-seq (default):** two-ended Tn5 insertion sites from proper paired
  fragments shorter than 150 bp, lenient MACS3 candidates, unscaled qpois
  signal refinement, context-level pooling, biological-replicate support, and
  a summit-aware master DHS registry across contexts.
- **ATAC-seq (optional):** MACS3 HMMRATAC, followed by the same pooled-context
  biological-replicate support step. HMMRATAC requires paired-end data.
- **TF ChIP-seq:** MACS3 narrow peaks and CPM BigWigs.
- **Histone ChIP-seq:** MACS3 broad peaks and CPM BigWigs.
- **ChIP controls:** a matched input/control can be named explicitly; IP-only
  ChIP is also valid and runs without `-c`.
- **QC:** FastQC before and after trimming, Cutadapt reports, alignment metrics,
  fragment-aware FRiP, ATAC TSS/fragment profiles, ChIP QC, and MultiQC.

The master DHS registry is the final ATAC processing endpoint. H3K27ac
integration, fixed-width ABC candidate generation, ABC scoring, and Micro-C
integration remain in the downstream atlas repository.

## Install

All environments remain inside the repository. Create the orchestration
environment at `.venv`:

```bash
git clone git@github.com:laiker96/short-read-processing.git
cd short-read-processing

export MAMBA_ROOT_PREFIX="$PWD/.micromamba"
export XDG_CACHE_HOME="$PWD/.cache"
micromamba create --prefix "$PWD/.venv" --file environment.yml -y
```

Run commands without activation:

```bash
micromamba run --prefix "$PWD/.venv" python src/run_pipeline.py --help
```

Snakemake creates rule-specific environments below `.snakemake/conda` through
the local profile. Bioinformatics packages are not installed globally.

## Input files

The primary public-data input is a CSV or TSV following
[`schemas/sample-sheet.schema.yaml`](schemas/sample-sheet.schema.yaml). One row
is one public `SRR`, `SRX`, `ERR`, or `ERX` accession.

Required columns:

| Column | Meaning |
|---|---|
| `accession` | Public run or experiment accession |
| `library_id` | Biological-library identifier; repeat it for technical runs |
| `assay` | `atac`, `h3k27ac`, `chip_tf`, or the `chip_histone` alias |
| `context` | Tissue, stage, or cell-type ID used for ATAC pooling |

The smallest valid ATAC table is:

```tsv
accession	library_id	assay	context
SRR100001	eye_atac_rep1	atac	eye
SRR100002	eye_atac_rep2	atac	eye
```

Distinct ATAC `library_id` values in the same `context` are biological
replicates. The pipeline calls replicate peaks, pools the context, and retains
pooled peaks supported by the configured number of libraries. Multiple
accessions sharing one `library_id` are technical runs and merge before
duplicate marking.

The genome is supplied once with `--genome` and defaults to `dm6`. FASTQ URLs
and paired/single-end layout are resolved from ENA/SRA. ATAC defaults to the
two-ended Tn5/MACS3-qpois method; add the optional `peak_caller` column and set
it to `hmmratac` to choose HMMRATAC.

H3K27ac is IP-only with the same four columns. For matched inputs, add `role`
and `control_library`:

```tsv
accession	library_id	assay	context	role	control_library
SRR200001	eye_h3_rep1	h3k27ac	eye	treatment	eye_input_rep1
SRR200002	eye_input_rep1	h3k27ac	eye	control
```

For IP-only ChIP, omit both optional columns. H3K27ac/histone ChIP defaults to
broad peaks; TF ChIP defaults to narrow peaks. The schema also defines optional
typed trimming, alignment, MACS3, and HMMRATAC override columns.

## Run

The canonical command validates the tables, downloads FASTQs concurrently,
writes a resolved YAML, prepares the reference, and starts Snakemake:

```bash
micromamba run --prefix "$PWD/.venv" \
  python src/run_pipeline.py samples.tsv \
  --project chromatin-study \
  --run-id baseline \
  --genome dm6 \
  --output-dir data/raw/chromatin-study \
  --reference-root references \
  --cores 24 \
  --file-jobs 8 \
  --connections 8
```

One mixed table produces a separate resolved workflow config for each assay.
ATAC contexts require at least two biological libraries by default, each
covering at least 50% of a pooled peak. Override these thresholds with
`--atac-minimum-replicates` and `--atac-overlap-fraction`.
For ATAC, the default `all` target automatically ends by building the master
DHS registry; no separate master-building command is required.

Useful boundaries:

```bash
# Download only
python src/run_pipeline.py samples.tsv --download-only --output-dir data/raw/project

# Reuse completed downloads and only write the resolved YAML
python src/run_pipeline.py samples.tsv --skip-download \
  --manifest data/raw/project/download_manifest.tsv --config-only

# Build the DAG without executing jobs
python src/run_pipeline.py samples.tsv --skip-download \
  --manifest data/raw/project/download_manifest.tsv --snakemake-dry-run
```

Run these through `micromamba run --prefix "$PWD/.venv"` as in the main
example.

## Reuse existing artifacts

The workflow starting stage is explicit:

| `--from-stage` | Required artifact | First scientific processing step |
|---|---|---|
| `accessions` | accession sample sheet | FASTQ acquisition |
| `final-bam` | complete final-BAM manifest | ATAC/ChIP downstream processing |
| `master` | complete master-DHS bundle manifest | bundle validation and registration |
| `activity` | master plus accepted atlas/reference BAM manifests | assay-unit counting and activity normalization |

`accessions` is the default and preserves the original behavior. Reuse modes
are strict: a missing, mismatched, rejected, or corrupt artifact is an error.
They never fall back to FASTQ download, trimming, alignment, duplicate marking,
or general BAM filtering.

### Start from filtered BAMs

The final-BAM manifest follows
[`schemas/final-bam-manifest.schema.yaml`](schemas/final-bam-manifest.schema.yaml).
Paths are relative to the manifest directory unless absolute. Every biological
library in each represented assay must occur exactly once. A manifest may
select only ATAC from a mixed accession sheet, but it may not contain only some
ATAC libraries.

Every successful sample-processing run automatically exports:

```text
results/<project>/<run-id>/provenance/manifests/final-bams.tsv
```

Pipeline-produced BAMs are recorded as `qc_status=pending_review`; completing
the workflow is not itself scientific QC acceptance. Both `pending_review` and
`accepted` BAMs can be reused for downstream peak calling and QC, while
`rejected` BAMs fail immediately. Activity quantification requires explicit
`accepted` status for both atlas and normalization-reference libraries.

For BAMs imported from outside this workflow, create the initial manifest from
exact `<library_id>.final.bam` filenames:

```bash
micromamba run --prefix "$PWD/.venv" \
  python src/create_final_bam_manifest.py \
  resources/atlas_samples_ip_only.tsv \
  --assay atac \
  --bam-dir results/atac \
  --output data/raw/drosophila-atlas/atlas-atac.final-bams.tsv \
  --genome dm6 --layout paired \
  --qc-status accepted \
  --source-project atlas-atac-dm6 --source-run-id imported-final-bams
```

The command computes full BAM/BAI SHA-256 values and writes the manifest
atomically without replacing identical content. Then build the qpois master
directly from those BAMs:

```bash
micromamba run --prefix "$PWD/.venv" \
  python src/run_pipeline.py resources/atlas_samples_ip_only.tsv \
  --from-stage final-bam \
  --final-bam-manifest data/raw/drosophila-atlas/atlas-atac.final-bams.tsv \
  --project drosophila-atlas --run-id qpois-master-v1 --genome dm6 \
  --config-dir data/raw/drosophila-atlas/configs/qpois-master-v1 \
  --cores 24
```

Before any peak or QC job consumes an imported BAM, a validation rule checks:

- the declared full BAM and BAI SHA-256 values;
- `samtools quickcheck`;
- readability of the explicitly declared index;
- coordinate sort order;
- exact chromosome names, order, and lengths against the selected reference;
- the `short-read-processing-final-v1` filtering contract and
  `qc_status=accepted`.

Validation receipts are written under
`provenance/external_bams/<library_id>.validated.json`. The imported BAMs and
indexes remain immutable and are never copied, indexed, touched, or rewritten.

### Reuse a frozen master set

The complete bundle, rather than only `master_dhs.bed`, is frozen with
[`schemas/master-manifest.schema.yaml`](schemas/master-manifest.schema.yaml):

An ATAC run that builds or validates a master registry automatically exports:

```text
results/<project>/<run-id>/provenance/manifests/master-dhs.tsv
```

Use `src/create_master_manifest.py` only to import a complete master bundle
that was created outside the current workflow namespace:

```bash
micromamba run --prefix "$PWD/.venv" \
  python src/create_master_manifest.py \
  results/drosophila-atlas/qpois-master-v1/atac/master \
  --output data/raw/drosophila-atlas/qpois-master-v1.master.tsv \
  --genome dm6 \
  --method reciprocal_summit_complete_linkage_v2 \
  --source-project drosophila-atlas --source-run-id qpois-master-v1
```

Validate and register it in another run namespace:

```bash
micromamba run --prefix "$PWD/.venv" \
  python src/run_pipeline.py resources/atlas_samples_ip_only.tsv \
  --from-stage master \
  --master-manifest data/raw/drosophila-atlas/qpois-master-v1.master.tsv \
  --project drosophila-atlas --run-id reuse-master-v1 --genome dm6 \
  --config-dir data/raw/drosophila-atlas/configs/reuse-master-v1 \
  --cores 2
```

Master mode performs bundle validation and provenance registration only.
`build_atac_master_dhs` is deliberately absent from this DAG.

### Build the ABC activity input

The activity stage consumes the frozen master set, accepted paired-end atlas
BAMs, and accepted paired-end normalization-reference BAMs. All artifacts are
explicit and checksummed; this stage has no download, trimming, alignment,
peak-calling, or master-reconstruction fallback.

```bash
micromamba run --prefix "$PWD/.venv" \
  python src/run_pipeline.py resources/atlas_samples_ip_only.tsv \
  --from-stage activity \
  --master-manifest \
    results/drosophila-atlas/qpois-master-v1/provenance/manifests/master-dhs.tsv \
  --activity-atlas-bam-manifest data/raw/drosophila-atlas/atlas-atac.accepted.tsv \
  --activity-atlas-bam-manifest data/raw/drosophila-atlas/atlas-h3k27ac.accepted.tsv \
  --activity-reference-sheet \
    data/raw/drosophila-s2-t0-reference/s2_t0_gse95689.tsv \
  --activity-reference-bam-manifest \
    data/raw/drosophila-s2-t0-reference/s2-t0-atac.accepted.tsv \
  --activity-reference-bam-manifest \
    data/raw/drosophila-s2-t0-reference/s2-t0-h3k27ac.accepted.tsv \
  --activity-reference-context s2_t0 \
  --project drosophila-atlas --run-id activity-s2-t0-v1 --genome dm6 \
  --config-dir data/raw/drosophila-atlas/configs/activity-s2-t0-v1 \
  --cores 24
```

Repeat either BAM-manifest option when ATAC and H3K27ac are stored in separate
manifests. The atlas accession sheet remains the positional input and supplies
library-to-context metadata; ChIP controls are ignored for activity
quantification.

For each accepted ATAC BAM, the workflow retains proper pairs with
`0 < abs(TLEN) < 150`, applies the standard two-ended Tn5 shift, and counts both
one-base insertion records. For each accepted H3K27ac BAM, it counts one
positive-TLEN fragment per proper pair. Raw counts are converted to
`CPM_per_kb = count * 10^9 / (total assay units * element width_bp)`, preserving
the variable-width master elements. Biological libraries are averaged with
equal weight within context and assay. ATAC and H3K27ac distributions are then
tie-aware quantile-normalized separately to the corresponding mean reference
profile. The combined value is `sqrt(atac_qnorm * h3k27ac_qnorm)`.

The canonical endpoint is:

```text
results/<project>/<run-id>/activity/master_dhs_activity.tsv.gz
```

It contains one row per master DHS and atlas context, including coordinates,
the original element width, contributing libraries, raw count sums and means,
library-depth/length-normalized values, replicate standard deviations,
assay-specific quantile-normalized values, and the combined activity value.
Supporting outputs are:

```text
activity/library_signal.tsv.gz             per-library raw and CPM_per_kb values
activity/context_signal.pre_qnorm.tsv.gz   context/assay aggregates before qnorm
activity/qnorm_reference.tsv.gz            frozen assay-specific reference profiles
activity/contexts/<context>.activity.tsv.gz
activity/activity_metrics.json
activity/activity_provenance.json
```

Every resolved configuration records sample-sheet, schema, and artifact
manifest hashes plus a timestamp-independent semantic SHA-256. Changing an
artifact, parameter, or scientific selection requires a new `run-id`; the
pipeline does not overwrite a scientifically different result namespace.

### Curated dm6 inputs

```bash
# Current atlas selection: ATAC plus IP-only H3K27ac
micromamba run --prefix "$PWD/.venv" \
  python src/run_pipeline.py resources/atlas_samples_ip_only.tsv \
  --project drosophila-atlas --run-id ip-only --genome dm6 \
  --output-dir data/raw/drosophila-atlas \
  --cores 24

# Alternative table containing available matched H3K27ac inputs
micromamba run --prefix "$PWD/.venv" \
  python src/run_pipeline.py resources/atlas_samples_with_inputs.tsv \
  --project drosophila-atlas --run-id matched-inputs --genome dm6 \
  --output-dir data/raw/drosophila-atlas \
  --cores 24

# Reprocess the D17 ATAC/H3K27ac comparator (not an automatic qnorm choice)
micromamba run --prefix "$PWD/.venv" \
  python src/run_pipeline.py resources/hq_cell_line_samples.tsv \
  --project drosophila-cell-line-reference --run-id d17-reference --genome dm6 \
  --output-dir data/raw/drosophila-cell-line-reference \
  --cores 24
```

Selection provenance is documented in
[`resources/README.md`](resources/README.md).

## ATAC default method

For paired-end ATAC, each biological library is processed as follows:

1. Retain proper, nonduplicate alignments with `0 < |TLEN| < 150`.
2. Apply the Tn5 offsets with `alignmentSieve --ATACshift`.
3. Convert both shifted mates to one-base insertion records.
4. Run MACS3 `callpeak -f BED -q 0.10 --nomodel --shift -75 --extsize
   150 --keep-dup all -B`.
5. Run `macs3 bdgcmp -m qpois` on the unscaled treatment pileup and local
   lambda. `--SPMR` is intentionally not used in this branch.
6. Progress from qpois exponent 2 through 325 and retain components 50–400 bp;
   broader components split as the threshold rises.
7. Concatenate replicate insertion records within each context and repeat
   candidate calling and refinement on the pool.
8. Retain a pooled peak when the configured number of replicate peak sets each
   cover the configured fraction of its bases.
9. Find each retained peak's summit in its pooled signal track and reconcile
   peaks across contexts into a variable-width master DHS registry.

Single-end ATAC follows the same insertion/qpois path without the unavailable
paired-fragment-length filter. HMMRATAC is an explicit paired-end alternative.

## Outputs

Each run is namespaced below `results/<project>/<run-id>/`.

ATAC context endpoints (the directory remains named `conditions` internally):

```text
atac/conditions/<context>/
  peaks/
    <context>.candidates.narrowPeak
    <context>.qpois-refined.bed
    <context>.qpois-excluded.bed
    <context>.qpois-refinement.json
    <context>.replicate-supported.bed       final context-level peak set
    <context>.replicate-support.tsv
    <context>.replicate-support.json
  tracks/
    <context>.MACS3-pileup.unscaled.bw
    <context>.qpois.bw
```

For HMMRATAC contexts, the pooled files are
`<context>.hmmratac.narrowPeak`, `<context>.CPM.bw`, and the same
`replicate-supported` BED/TSV/JSON outputs.

Qpois-refined BEDs contain BED6 followed by maximum qpois score and selection
exponent. Replicate-supported BEDs contain BED6 followed by `condition_id`,
`support_n`, `replicate_n`, `support_fraction`, comma-separated supporting
library IDs, and `peak_method`.

The final cross-context ATAC outputs are:

```text
atac/master/
  master_dhs.bed                  strict BED6 variable-width master intervals
  master_dhs_summits.bed          one-base representative summits
  master_dhs_membership.tsv       every contributing context peak
  master_dhs_context_matrix.tsv   context presence for each master DHS
  master_dhs.json                 parameters and summary statistics
```

For qpois contexts, each source summit is the center of the maximum plateau in
the pooled unscaled MACS3 pileup within that refined peak. HMMRATAC contexts use
their pooled CPM BigWig. If an interval contains no finite signal, its midpoint
is used and recorded as a fallback. A source interval extending beyond a
reference contig is clipped to the contig boundary; the original coordinates
and clipping flag remain in `master_dhs_membership.tsv`.

Source peaks are considered the same DHS only when each peak contains the
other's summit, their complete summit span is at most 150 bp (recorded as
`atac_master.summit_max_distance`), and the cluster does not already contain a
peak from that context. Narrow peaks are considered first, so a broad peak from
one context is assigned only to the narrow DHS containing its maximum and
cannot collapse two sites resolved in another. The representative summit is
the observed source summit nearest the median of the contributing source
summits (with deterministic ties).

After this initial clustering, adjacent clusters with representative summits
less than 50 bp apart are treated as context-shifted calls of the same DHS and
merged when their context sets are disjoint. They remain separate when at least
one context contributes a source peak to both clusters, because that context
independently resolved two sites. The closest eligible pair merges first, and
the combined source-summit span must still be at most 150 bp. The 50 bp rule is
recorded as `atac_master.minimum_summit_separation`. Consequently, the default
qpois workflow does not pad a boundary-clipped master DHS merely to reach 50
bp. This setting is a minimum separation between representative summits, not a
minimum final interval width: a sub-50-bp interval may remain when
shared-context evidence resolves two nearby sites, or when midpoint clipping
trims an asymmetric source-peak envelope even though neighboring summits are
at least 50 bp apart.

Final boundaries are the envelope of contributing refined peaks and are
clipped at the midpoint between adjacent master summits only when their
envelopes overlap. This step never resizes DHSs to 500 bp; standardized ABC
windows are constructed downstream.

ChIP endpoints:

```text
peaks/<sample>/
  <sample>_peaks.narrowPeak       TF ChIP
  <sample>_peaks.broadPeak        histone ChIP
  <sample>_treat_pileup.bdg
  <sample>_control_lambda.bdg
tracks/<sample>.CPM.bw
```

ChIP `callpeak` uses `-B --SPMR`; `-c` is added only when `control_library` is
present. Alignment BAMs are retained for reproducible downstream reruns, while
replicate-only ATAC peak evidence and insertion files live below `work/`.

Shared outputs include:

```text
bam/                         filtered, indexed alignments
qc/fastqc/                   raw and trimmed FastQC
qc/cutadapt/                 trimming reports
qc/alignment/                SAMtools statistics
qc/frip/                     numerator, denominator, and FRiP
qc/tss/ and qc/fragments/    ATAC QC
qc/chip/                     ChIP fingerprint/cross-correlation
qc/metrics.tsv and .json     stable machine-readable summary
qc/multiqc/                  aggregate HTML report
provenance/resolved_config.json
logs/
```

## Restartability and parallelism

Re-run the identical command to resume:

- aria2 resumes partial ENA downloads and validates reported checksums;
- aria2 retries only checksum-failed/incomplete FASTQs three times by default;
  use `--checksum-retries` to change the bounded retry count;
- SRA conversion promotes FASTQs only after successful completion;
- reference preparation and every processing stage are Snakemake outputs;
- temporary scientific outputs are written in staging paths before promotion;
- completed alignments are reused when peak parameters change;
- validated external BAMs and master bundles are immutable workflow inputs;
- successful stages automatically export reusable final-BAM and master-bundle
  manifests under `provenance/manifests/`;
- identical manifests/configurations retain stable semantic hashes and do not
  replace unchanged resolved YAML;
- a changed scientific parameter set should use a new `run-id`.

Independent accessions, technical lanes, biological libraries, and contexts
are separate jobs. `--cores` limits aggregate CPU usage; each rule separately
declares threads and memory. For downloads, `--file-jobs` is concurrent files
and `--connections` is segmented connections per file.

## SLURM

All site-specific launchers and profiles belong under the ignored `slurm/`
directory. Do not run downloads, alignment, peak calling, or environment
installation on a login node.

```bash
micromamba run --prefix "$PWD/.venv" \
  python src/run_pipeline.py samples.tsv \
  --workflow-profile slurm/profile \
  --jobs 50 --cores 200 --max-threads 16
```

`--jobs` caps submitted/running jobs, `--cores` caps aggregate requested CPUs,
and `--max-threads` caps one rule. Cluster hostnames, accounts, partitions, and
paths must remain in ignored files under `slurm/`.

## IGV session

Build a portable session from the final ATAC contexts and optional ChIP run:

```bash
micromamba run --prefix "$PWD/.venv" \
  python src/build_igv_session.py \
  results/drosophila-atlas.atac.dm6/ip-only/atac \
  --chip-root results/drosophila-atlas.chip_histone.dm6/ip-only \
  --output results/atlas.igv.xml --genome dm6 \
  --final-atac-only --chip-one-per-context
```

The final-only view contains pooled ATAC pileup/qpois tracks and the
replicate-supported ATAC peaks. The ChIP context option deterministically
selects the first sorted replicate (normally `rep1`) for each context.
When `atac/master/master_dhs.bed` exists, it is added automatically as the
first feature track; use `--master-bed` to select another registry explicitly.

## Verification

```bash
export MAMBA_ROOT_PREFIX="$PWD/.micromamba"
export XDG_CACHE_HOME="$PWD/.cache"

micromamba run --prefix "$PWD/.venv" pytest -q
micromamba run --prefix "$PWD/.venv" \
  snakemake --snakefile workflow/Snakefile \
  --configfile tests/fixtures/workflow_config.yaml --lint
micromamba run --prefix "$PWD/.venv" \
  snakemake --snakefile workflow/Snakefile \
  --configfile docs/workflow-dag.config.yaml --cores 16 --dry-run
```
