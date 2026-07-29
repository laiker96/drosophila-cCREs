# Master DHS and regulatory-catalog implementation plan

## Objective

Build an idempotent, reproducible pipeline that:

1. constructs a variable-width, summit-aware master DHS registry from every
   ATAC context declared in the canonical sample table;
2. quantifies accepted ATAC and H3K27ac libraries without downloading or
   realigning when reusable BAM artifacts are supplied;
3. applies assay-specific background TMM normalization;
4. fits guarded two-component H3K27ac mixtures within context-member DHSs;
5. reports low and high component assignments, explicit guardrail failures,
   one wide master table, and active element sets per context;
6. renders a final integrated HTML/PDF QC report from frozen upstream and
   current manifests, JSON metrics, tables, and plots.

External-reference quantile normalization and cell-line reference cohorts are
not part of this design.

## Stage contract

| Stage | Inputs | Durable result |
|---|---|---|
| trimming | downloaded FASTQs | trimmed FASTQs and lane QC |
| alignment | trimmed FASTQs | filtered final BAM/BAI and checksummed manifest |
| QC | final BAMs | peaks, FRiP, assay QC, MultiQC, and a `pending_review` BAM manifest |
| master | fully reviewed ATAC BAM manifest | accepted-library replicate/context evidence, master DHS bundle, and manifest |
| quantification | master manifest plus accepted ATAC/H3K27ac BAM manifests | raw, CPM/kb and background-TMM activity tables |
| catalog | completed quantification | guarded mixtures, long/wide catalogs, active sets, BED/BigWig tracks and IGV sessions |
| report | completed catalog plus frozen upstream QC artifacts | integrated HTML/PDF report and checksummed JSON sidecar |

The output stopping point is not part of the scientific semantic hash. Within
one project/run namespace, requesting a later stage reuses every complete
upstream Snakemake output. Strict cross-run entry points validate checksums and
never fall back to upstream download/alignment.

Master construction is a manual QC checkpoint, not an interactive workflow
pause. Accession runs stop no later than QC. The user records one accepted or
rejected decision per ATAC library and starts a separate strict `final-bam`
run. Pending decisions are errors; rejected rows remain in the reviewed input
manifest, are excluded from peak/master inputs, and are recorded in resolved
provenance. Every retained context must still meet the configured biological-
replicate minimum.

## Master DHS construction

1. Resolve technical runs to biological `library_id` values.
2. Process lanes independently, merge and duplicate-mark by library, run QC,
   and export the final-BAM manifest for manual review.
3. Require a complete reviewed ATAC manifest and select only accepted
   biological libraries.
4. For paired-end ATAC qpois, retain proper pairs satisfying
   `0 < abs(TLEN) < 150`, apply Tn5 offsets, and represent both mates as
   one-base insertions.
5. Call and refine qpois peaks independently per biological library.
6. Pool libraries by sample-table `context` and retain peaks with configured
   replicate support.
7. Cluster context peaks through reciprocal summit containment and
   complete-linkage summit distance, with at most one peak per context.
8. Write master BED, summit BED, membership table, wide context matrix, metrics,
   and a checksummed reusable manifest.

Verification: test context discovery, minimum replicate validation, interval
ordering/non-overlap, summit containment, deterministic IDs, and the dry-run
master branch.

## Quantification

### Assay units

- ATAC: Tn5-shifted one-base insertions from proper pairs with
  `0 < abs(TLEN) < 150`.
- paired-end H3K27ac: one genomic fragment for each positive-TLEN proper pair.
- single-end H3K27ac: extended fragment using the reviewed manifest value.

### Exact master-element counts

Count each library over the variable-width master intervals and save:

- raw count;
- total usable assay units;
- CPM/kb using the original element width.

### Background TMM

1. Create deterministic 10-kb bins on configured autosomes.
2. Count each library's assay units in those bins.
3. Estimate edgeR TMM factors separately by assay.
4. Apply each factor to exact master-interval counts.
5. Aggregate biological libraries by equal-weight context mean.

Save the input matrix, metadata, factors, software receipt, normalized context
signals, master-element activity table, metrics, and checksummed provenance.

Verification: test bin boundaries, matrix identity/order, assay-specific
factors, raw/CPM/kb preservation, deterministic output, and a quantification-only
workflow dry run.

## H3K27ac summit-window signal

For each master summit construct clipped, zero-based half-open windows:

- left: `[summit-750, summit-250)`;
- center: `[summit-250, summit+250)`;
- right: `[summit+250, summit+750)`.

Count raw H3K27ac fragments in every window, apply the canonical background-TMM
factor, average biological libraries within context, and only then take the
maximum of the three context means. Resolve ties center, left, then right.

ATAC remains quantified over the exact variable-width master interval.

Verification: test chromosome clipping, zero-width windows, replicate averaging
before maximum selection, and deterministic tie handling.

## Guarded mixture model

Fit deterministic one- and two-Gaussian models to positive log10 maximum-window
H3K27ac values among DHSs with membership in that context.

A context is supported only when all conditions hold:

- at least 200 positive member DHSs;
- `BIC_1G - BIC_2G >= 10`;
- Ashman's `D >= 2`;
- both component weights are at least 0.10;
- exactly one posterior-0.5 crossing occurs between fitted means;
- the input variance is greater than `1e-12`.

If a two-component fit exists, assign positive member DHSs to `low` or `high`
at posterior probability 0.5 regardless of support. Unsupported assignments
must carry `mixture_supported=0`, `mixture_guardrail_warning=1`, and every
failed guard name. If no fit can be estimated, use `not_applicable` and report
the reason.

Verification: test separated and unimodal synthetic populations, deterministic
fits, every guard, low/high assignment from unsupported fits, and failure text.

## Regulatory classes

Use absolute summit-to-nearest-TSS distance:

- promoter-associated: at most 250 bp;
- proximal enhancer-like: 251–1,000 bp;
- distal enhancer-like: more than 1,000 bp;
- unclassified: no TSS on the element's contig.

An active element is a context-member DHS assigned to the high component.
Mixture support is an orthogonal quality flag and never silently removes an
unsupported assignment.

## Final outputs

The catalog stage writes:

1. `master_elements_long.tsv.gz`: one element/context row with membership,
   TSS class, raw/normalized signals, both mixture assignments and warnings;
2. `master_elements_wide.tsv.gz`: one master-element row with context-prefixed
   membership, signal, mixture, warning and activity columns;
3. `active/<context>.active_elements.tsv.gz`: high-component member DHSs with
   their promoter/proximal/distal class and guardrail status;
4. `mixture_models.tsv`: fitted parameters and exact support reasons;
5. mixture distribution SVG and machine-readable histogram bins;
6. per-context active-element and context-DHS BED9 tracks plus a copied master
   BED;
7. background-TMM-normalized mean ATAC and H3K27ac BigWigs, a portable
   five-track IGV session per context, and one session containing every
   context;
8. summary, metrics, and checksummed provenance files.

The continuous normalized ATAC and maximum-window H3K27ac values are retained
for ABC scoring; the binary catalog labels are annotations rather than a
replacement for continuous activity.

## Integrated QC report

The final report inventories the accession/sample sheet, reviewed BAM
manifests, master manifest, resolved upstream configurations, QC metrics,
quantification outputs, and catalog outputs. It embeds library-yield and FRiP
plots, H3K27ac cross-correlation summaries, ATAC TSS profiles, master context
statistics, TMM factors, active TSS-distance classes, and the guarded mixture
plot. Every supplied report input is checksum-validated, and an accompanying
JSON file records the complete audit trail. Reporting inputs do not affect the
scientific semantic hash.

## Final verification

Run:

```bash
export MAMBA_ROOT_PREFIX="$PWD/.micromamba"
export XDG_CACHE_HOME="$PWD/.cache"
mamba run --prefix "$PWD/.venv" pytest -q
mamba run --prefix "$PWD/.venv" snakemake --snakefile workflow/Snakefile \
  --configfile tests/fixtures/workflow_config.yaml --lint
mamba run --prefix "$PWD/.venv" snakemake --snakefile workflow/Snakefile \
  --configfile tests/fixtures/workflow_config.yaml --cores 8 --dry-run
git diff --check
```

Also dry-run `tests/fixtures/activity_config.yaml` to assert the
`quantification`, `catalog`, and `report` branches.
