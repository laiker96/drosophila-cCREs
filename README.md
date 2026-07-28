# Drosophila cCREs

Reproducible ATAC-seq and ChIP-seq processing from public accessions to a
summit-aware master DHS registry and a context-resolved regulatory-element
catalog. Contexts are not hard-coded: they come from the canonical input table.

The production path has three phases separated by an explicit manual QC gate:

1. process reads through library QC and export `pending_review` BAM manifests;
2. review every library, then construct the master DHS set from accepted ATAC
   libraries only;
3. quantify accepted ATAC/H3K27ac BAMs, normalize with background TMM, fit the
   guarded H3K27ac mixtures, and write regulatory-element catalogs.

There is no external-reference or quantile-normalization branch.

## Environment

Create the repository-local orchestration environment:

```bash
export MAMBA_ROOT_PREFIX="$PWD/.micromamba"
mamba env create --prefix "$PWD/.venv" --file environment.yml
```

Snakemake rule environments are created below `.snakemake/conda` by the local
profile. Do not create named environments outside the repository.

## Input table

The public input is a CSV or TSV conforming to
`schemas/sample-sheet.schema.yaml`. Its canonical columns are:

| Column | Meaning |
|---|---|
| `accession` | Public run or experiment accession |
| `library_id` | Biological library; technical runs share an ID |
| `assay` | `atac`, `h3k27ac`, or another supported ChIP assay |
| `context` | Tissue, stage, or cell type |

Optional control and peak-caller columns are documented in the schema. The two
reviewed dm6 atlas inputs are `resources/atlas_samples_ip_only.tsv` and
`resources/atlas_samples_with_inputs.tsv`.

Each distinct ATAC `context` becomes a replicate-pooling group. Consequently,
the same code can build a master from the nine atlas contexts or any other set
of contexts present in an input table. The default consensus requires at least
two biological ATAC libraries per context.

## Logical stages and restart behavior

`--until-stage` selects a reproducible stopping boundary:

| Stage | Completed result |
|---|---|
| `trimming` | raw/trimmed FastQC, Cutadapt metrics, trimmed FASTQs |
| `alignment` | filtered final BAM/BAI, alignment QC, final-BAM manifest |
| `qc` | peaks, FRiP and assay QC, MultiQC |
| `master` | replicate-supported context peaks and master DHS bundle/manifest |
| `quantification` | raw, CPM/kb, background-TMM factors and normalized master-element signals |
| `catalog` | max-window H3K27ac mixtures, long/wide catalogs, active sets, BED/BigWig tracks, and IGV sessions |
| `report` | integrated, checksummed HTML/PDF QC report spanning inputs through the catalog |

Snakemake owns completeness. Re-running the same `project`/`run_id` resumes
from existing valid outputs and never realigns merely because a later stage is
requested. Advancing `--until-stage` does not change the scientific semantic
digest. Incomplete outputs are rebuilt because the profiles retain
`rerun-incomplete: true`.

Cross-run reuse is strict and manifest-based:

| `--from-stage` | Required artifact | Allowed stopping points |
|---|---|---|
| `accessions` | accession table/download manifest | trimming through QC |
| `final-bam` | `--final-bam-manifest` | alignment or QC with pending/accepted BAMs; master with a fully reviewed ATAC-only manifest |
| `master` | `--master-manifest` | validate/export the master bundle |
| `quantification` | master plus accepted ATAC/H3K27ac BAM manifests | quantification, catalog, or report |

The workflow automatically exports:

```text
results/<project>/<run_id>/provenance/manifests/final-bams.tsv
results/<project>/<run_id>/qc/library-review.tsv
results/<project>/<run_id>/provenance/manifests/master-dhs.tsv
```

Automatically exported new BAMs have `qc_status=pending_review`. The QC endpoint
also writes a single human-readable row per library to `qc/library-review.tsv`.
Apply explicit pass/fail decisions with `src/review_final_bam_manifest.py`.
Master DHS
construction refuses any `pending_review` ATAC row, excludes documented
rejections, and records those exclusions in the resolved configuration. Each
remaining ATAC context must still satisfy `--atac-minimum-replicates`.
Quantification likewise requires reviewed ATAC and H3K27ac manifests. Rejected
libraries require a reason. Single-end H3K27ac libraries also require a
reviewed `estimated_fragment_length_bp`.

## Phase 1: master DHS construction

First process accessions through QC. Omitting `--until-stage` has the same QC
endpoint in accession mode:

```bash
python src/run_pipeline.py resources/atlas_samples_ip_only.tsv \
  --project drosophila-atlas --run-id qc-v1 --genome dm6 \
  --until-stage qc --cores 16
```

This first run:

1. downloads immutable FASTQs with checksums and resume state;
2. runs lane-level QC/trimming/alignment;
3. merges technical runs by `library_id`, marks duplicates, and filters BAMs;
4. calls replicate peaks and produces FRiP, ATAC TSS/fragment, H3K27ac
   cross-correlation, and MultiQC outputs;
5. exports a final-BAM manifest whose new libraries are `pending_review` and a
   populated `qc/library-review.tsv` for manual review.

The review table contains context and layout, final-BAM depth, proper-pair
fraction, FRiP, peak count, assay-specific metrics, and paths to MultiQC and
the per-library plots. ATAC rows include TSS enrichment, fragment-length
median, and fractions below 150 bp and from 180--250 bp. TSS enrichment is the
maximum aggregate signal in the central +/-50 bp divided by the mean signal in
the two terminal 100-bp flanks of the +/-2-kb profile. Histone/ChIP rows include
phantompeakqualtools fragment-shift estimates, NSC, RSC, and quality tag.

Copy the generated table outside the Snakemake result namespace, then edit only
the final three review columns. Set `qc_decision` to `pass` or `fail`; a failed
row requires a reason in `notes`. The suggested phantompeakqualtools shift is
pre-populated as `estimated_fragment_length_bp` for single-end H3K27ac and
should be reviewed rather than accepted blindly.

```bash
mkdir -p data/reviewed
cp results/<qc-project>/qc-v1/qc/library-review.tsv \
  data/reviewed/atlas-atac.library-review.tsv
```

Apply the decisions without editing the generated manifest in place:

```bash
python src/review_final_bam_manifest.py \
  results/<qc-project>/qc-v1/provenance/manifests/final-bams.tsv \
  --review-table data/reviewed/atlas-atac.library-review.tsv \
  --output data/reviewed/atlas-atac.reviewed.tsv
```

The command maps `pass` to final-manifest `qc_status=accepted` and `fail` to
`qc_status=rejected`, verifies that every library has exactly one decision, and
does not alter the original manifest. The older `--decisions` option and
`accepted`/`rejected` values remain supported for compatibility.

Then construct the master from the reviewed ATAC manifest. This strict reuse
mode cannot download FASTQs, trim reads, or align again:

```bash
python src/run_pipeline.py resources/atlas_samples_ip_only.tsv \
  --from-stage final-bam \
  --final-bam-manifest data/reviewed/atlas-atac.reviewed.tsv \
  --project drosophila-atlas --run-id master-v1 --genome dm6 \
  --until-stage master --cores 16
```

This second run calls/refines peaks for accepted ATAC libraries, pools them
within each context, applies replicate support, and builds the summit-aware,
variable-width master registry. Rejected rows remain in the reviewed input
manifest and are copied into `provenance.excluded_master_libraries` in the
resolved configuration.

The master is not a simple interval union or fixed-width resize. It preserves
representative summits, reciprocal-summit clustering, the configured maximum
summit span, and context membership.

To resume the pre-review QC run, repeat it with the same project/run ID and
download manifest, for example:

```bash
python src/run_pipeline.py resources/atlas_samples_ip_only.tsv \
  --project drosophila-atlas --run-id qc-v1 --genome dm6 \
  --skip-download --manifest data/raw/download_manifest.tsv \
  --until-stage qc --cores 16
```

There is no interactive pause inside Snakemake. The completed QC job and the
review command form the explicit checkpoint; a separate master invocation is
required. A master invocation with any pending ATAC decision fails before
Snakemake starts.

## Phase 2: quantification and catalog

Use the same accession table as metadata, the exported master manifest, and
one or more reviewed final-BAM manifests:

```bash
python src/run_pipeline.py resources/atlas_samples_ip_only.tsv \
  --from-stage quantification \
  --master-manifest data/raw/drosophila-atlas/master-dhs.tsv \
  --activity-bam-manifest data/raw/drosophila-atlas/atac.accepted.tsv \
  --activity-bam-manifest data/raw/drosophila-atlas/h3k27ac.accepted.tsv \
  --report-source-root results/drosophila-atlas/master-v1 \
  --report-source-root results/drosophila-atlas-h3k27ac/qc-v1 \
  --project drosophila-atlas --run-id catalog-v1 --genome dm6 \
  --until-stage report --cores 16
```

This mode validates hashes and BAM/reference compatibility before computation.
It cannot download FASTQs or invoke alignment.

Use `--until-stage quantification` to stop after counting and normalization.
Re-run the same command with `--until-stage catalog` to reuse those outputs and
continue with the mixture/catalog stage, or with `--until-stage report` to add
the final integrated report. Repeated `--report-source-root` values freeze the
upstream manifests, JSON metrics, TSS profiles, fragment histograms,
cross-correlation results, and MultiQC locations into the report configuration.
The master-manifest result root is discovered automatically when possible.
Report inputs affect only reporting and are excluded from the scientific
semantic digest.

## Quantification method

Assay units are prepared independently:

- paired-end ATAC retains proper pairs with `0 < abs(TLEN) < 150`, applies the
  Tn5 shift, and counts one-base insertions over each exact variable-width
  master DHS;
- paired-end H3K27ac counts one fragment per proper pair;
- single-end H3K27ac extends reads using the reviewed fragment length.

The pipeline saves per-library raw counts and CPM/kb. It estimates TMM factors
separately for ATAC and H3K27ac from raw counts in fixed 10-kb autosomal
background bins. Context values are equal-weight means of normalized biological
libraries. The quantification table retains raw, CPM/kb, normalized CPM/kb,
replicate SD, and `sqrt(ATAC × H3K27ac)` values.

## H3K27ac signal and guarded mixtures

For every master summit, H3K27ac is counted in three clipped, non-overlapping
500-bp windows:

```text
left:   [summit-750, summit-250)
center: [summit-250, summit+250)
right:  [summit+250, summit+750)
```

Replicates are averaged within context before selecting the maximum normalized
window. The two-Gaussian model is fitted to positive `log10(max-window
H3K27ac)` values among DHSs open in that context.

A fit is marked supported only if all guards pass:

- at least 200 positive member DHSs;
- `BIC(one Gaussian) - BIC(two Gaussians) >= 10`;
- Ashman's D is at least 2;
- both component weights are at least 0.10;
- exactly one posterior-0.5 crossing lies between the component means;
- the fitted population has nonzero variance.

Whenever a two-component fit exists, every positive member DHS is assigned to
`low` or `high` using posterior probability 0.5, even if a guard fails. Such
rows are retained with `mixture_supported=0`,
`mixture_guardrail_warning=1`, and the exact semicolon-separated failures in
`mixture_guardrail_failures`. No unsupported call is silently promoted to a
supported one.

## Regulatory classes

Distance is measured from the master DHS summit to the nearest reference TSS:

| Class | Summit-to-TSS distance |
|---|---:|
| `promoter_associated` | ≤250 bp |
| `proximal_enhancer_like` | 251–1,000 bp |
| `distal_enhancer_like` | >1,000 bp |
| `unclassified_no_tss_on_contig` | no TSS on that contig |

An active element is an open DHS assigned to the high H3K27ac component. This
binary annotation is intended for an encyclopedia-style catalog; continuous
ATAC, H3K27ac, and combined activity values remain the appropriate inputs for
ABC scoring.

## Browser tracks and IGV sessions

The catalog stage creates a self-contained visualization bundle for every
context. The context DHS BED uses the master-DHS coordinates for rows with
`context_matrix` membership equal to one; it is therefore directly aligned
with the catalog rather than reproducing the wider upstream pooled-peak
boundaries. The active-element BED is BED9: the thick one-base interval marks
the representative summit, the score is the high-component posterior scaled
to 0--1,000, and item RGB distinguishes promoter-associated, proximal,
distal, and unclassified elements. Its name records any mixture-guardrail
warning.

For each assay and context, every accepted library is first scaled by
`1e6 / effective_library_size` using the same assay-specific background-TMM
factor as the activity table. The context BigWig is the arithmetic mean of
those basewise tracks. ATAC values represent Tn5 insertion coverage;
H3K27ac values represent inferred fragment coverage. Consequently the browser
signals and quantitative catalog share the same library selection and
normalization, while retaining their assay-specific unit semantics.

Opening `activity/catalog/igv/<context>.xml` loads five relative, portable
resources: mean ATAC, mean H3K27ac, context DHSs, active elements, and the
global master-DHS registry. Keep the `bed/`, `tracks/`, and `igv/` directories
together if the catalog is moved.

## Outputs

Master bundle:

```text
results/<project>/<run_id>/atac/master/master_dhs.bed
results/<project>/<run_id>/atac/master/master_dhs_summits.bed
results/<project>/<run_id>/atac/master/master_dhs_membership.tsv
results/<project>/<run_id>/atac/master/master_dhs_context_matrix.tsv
results/<project>/<run_id>/atac/master/master_dhs.json
```

Quantification:

```text
results/<project>/<run_id>/activity/quantification/libraries/*.signal.tsv.gz
results/<project>/<run_id>/activity/quantification/tmm_input_counts.tsv.gz
results/<project>/<run_id>/activity/quantification/normalization_factors.tsv
results/<project>/<run_id>/activity/quantification/context_signal.tsv.gz
results/<project>/<run_id>/activity/quantification/master_dhs_activity.tsv.gz
results/<project>/<run_id>/activity/quantification/activity_metrics.json
results/<project>/<run_id>/activity/quantification/activity_provenance.json
```

Catalog:

```text
results/<project>/<run_id>/activity/catalog/master_elements_long.tsv.gz
results/<project>/<run_id>/activity/catalog/master_elements_wide.tsv.gz
results/<project>/<run_id>/activity/catalog/active/<context>.active_elements.tsv.gz
results/<project>/<run_id>/activity/catalog/mixture_models.tsv
results/<project>/<run_id>/activity/catalog/regulatory_element_summary.tsv
results/<project>/<run_id>/activity/catalog/h3k27ac_mixture_distributions.svg
results/<project>/<run_id>/activity/catalog/regulatory_element_metrics.json
results/<project>/<run_id>/activity/catalog/regulatory_element_provenance.json
results/<project>/<run_id>/activity/catalog/bed/master_dhs.bed
results/<project>/<run_id>/activity/catalog/bed/<context>.dhs.bed
results/<project>/<run_id>/activity/catalog/bed/<context>.active_elements.bed
results/<project>/<run_id>/activity/catalog/bed/bed_tracks.json
results/<project>/<run_id>/activity/catalog/tracks/<context>.atac.mean.background_tmm.bw
results/<project>/<run_id>/activity/catalog/tracks/<context>.h3k27ac.mean.background_tmm.bw
results/<project>/<run_id>/activity/catalog/tracks/<context>.<assay>.mean.background_tmm.json
results/<project>/<run_id>/activity/catalog/igv/<context>.xml
```

Integrated report:

```text
results/<project>/<run_id>/activity/report/integrated_qc_report.html
results/<project>/<run_id>/activity/report/integrated_qc_report.pdf
results/<project>/<run_id>/activity/report/integrated_qc_report.json
```

The report includes frozen input and output inventories with checksums, BAM QC
and FRiP statistics, H3K27ac cross-correlation, ATAC TSS plots, master-registry
statistics, TMM factors, active-element TSS classes, mixture fits and exact
guardrail warnings. The JSON sidecar records every report input and output hash.

The long table contains one master-element/context row and reports both mixture
components. The wide table contains one row per master element with
context-prefixed membership, signal, mixture, warning, and activity columns.
Each per-context active file contains the high-component rows, including
unsupported-fit assignments with their explicit warning fields. The BED and
mean-signal sidecars are also explicit report inputs, so their paths and
checksums appear in the integrated audit trail.

## Verification

```bash
export MAMBA_ROOT_PREFIX="$PWD/.micromamba"
export XDG_CACHE_HOME="$PWD/.cache"

mamba run --prefix "$PWD/.venv" pytest -q
mamba run --prefix "$PWD/.venv" \
  snakemake --snakefile workflow/Snakefile \
  --configfile tests/fixtures/workflow_config.yaml --lint
mamba run --prefix "$PWD/.venv" \
  snakemake --snakefile workflow/Snakefile \
  --configfile tests/fixtures/workflow_config.yaml --cores 8 --dry-run
git diff --check
```
