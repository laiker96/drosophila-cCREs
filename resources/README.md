# Pipeline resources

This directory contains reviewed accession inputs for the dm6 atlas and
cell-line reference runs, plus the adapter sequences used by the workflow.

## Atlas inputs

- `atlas_samples_ip_only.tsv`: minimal four-column pipeline input with 23 ATAC
  accession rows and 15 H3K27ac IP rows. It contains no input-DNA libraries.
- `atlas_samples_with_inputs.tsv`: alternative combined input containing the
  same ATAC selection and 15 H3K27ac IPs, with matched inputs for the two
  GSE140539 stage-5 IPs and the selected adult-brain IP. Other H3K27ac contexts
  remain IP-only because no suitable matched input is included.

These two reviewed tables are the versioned pipeline inputs. Detailed dataset
selection and QC metadata belong to the separate atlas-analysis repository and
are not required by this accession-processing workflow.

Runtime final-BAM and master-bundle manifests contain filesystem paths and
checksums and therefore belong under the ignored `data/raw/<project>/`
namespace, not in this versioned resource directory.

Both tables use the default two-ended Tn5/MACS3-qpois ATAC branch. The
with-inputs table deliberately substitutes the matched GSE140539 stage-5
H3K27ac IPs for the IP-only PRJEB37091 stage-5 pair; it is an alternative
scientific selection, not merely the IP-only table plus control rows.
The encoded matched pairs are GSE140539 stage-5 IP/input runs
`SRR10485675/SRR10485676` and `SRR10485677/SRR10485678`, plus adult-brain
`SRR5319052/SRR5319047`.

## High-quality cell-line reference candidate

- `hq_cell_line_samples.tsv`: paired ATAC-seq and H3K27ac ChIP-seq for the
  Drosophila ML-DmD17-c3 (D17) cell line, with two replicate libraries per
  assay. The table deliberately contains the H3K27ac IPs only.

The four selected libraries come from the same NovaSeq 6000 multi-omics study,
[GSE245079](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245079).
The ATAC samples are in
[GSE245076](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245076), and
the H3K27ac samples are in
[GSE245077](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245077).
All are paired-end 50-bp libraries aligned by the submitters to dm6/BDGP6.

| Assay | Library | GEO sample | SRA run | Read pairs |
|---|---|---|---|---:|
| ATAC | replicate 1 | GSM7835773 | SRR26349810 | 25,184,931 |
| ATAC | replicate 2 | GSM7835774 | SRR26349809 | 16,606,924 |
| H3K27ac | replicate 1 | GSM7835791 | SRR26350407 | 13,579,699 |
| H3K27ac | replicate 2 | GSM7835792 | SRR26350406 | 14,913,574 |

This is a modern, within-study comparator because both assays use the exact
same named cell line and the study provides replicate libraries. It is not
selected automatically as the downstream ABC normalization reference. The
table is directly accepted by `src/run_pipeline.py`; H3K27ac runs IP-only. The
study's available ChIP input is intentionally excluded because the planned
activity table quantifies IP fragments directly over the ATAC-derived master
DHS set and does not perform input subtraction or input-controlled peak
calling.

“High quality” is provisional until the raw reads pass the local pipeline QC.
Before using this reference for quantile normalization, require acceptable
mapping and duplication rates, ATAC TSS enrichment/FRiP and replicate
agreement, and H3K27ac FRiP and replicate agreement. The table does not make
the processing workflow perform ABC scoring or quantile normalization.
GEO explicitly describes the ChIP replicates as coming from different flasks,
but does not state whether the two ATAC libraries came from independent
cultures. The table assigns the ATAC samples distinct `library_id` values, so
the pipeline treats them as biological replicates; confirm that interpretation
before relying on replicate support as strictly biological.

The leading metadata-level matched candidate is untreated S2 at T0 from
[GSE95689](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE95689), which
has two ATAC-seq and two conventional H3K27ac replicates from the same cell
line, study, state, and time point. Its raw paired-end reads must be remapped
from the submitter's dm3 analysis to dm6 and pass the same local QC before it
can be selected.

Untreated S2 and Kc ATAC replicates from
[GSE119708](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119708) can
also be compared with the S2-DRSC and Kc167 H3K27ac libraries in GSE245077.
Those are cross-study sensitivity references: S2 must not be assumed identical
to S2-DRSC, nor Kc to Kc167. Candidate selection and QC are scientific run
metadata; they must be frozen explicitly and must never trigger a silent
fallback between cohorts.

## Adapter sequences

The 158 sequence records in `adapters/adapters.fa` come from
`BBMap_39.96.tar.gz:bbmap/resources/adapters.fa` in the parent workspace.

- BBMap archive SHA-256: `e173bdd0d3ca047f378c71dad568a148596c1690bf36abca93e918569c9fb382`
- Upstream extracted FASTA SHA-256: `85abe9d3e40dc37c968f7e4c1227e05976a4ed0583d1dd442d375aa7516f13a9`
- Repository FASTA SHA-256: `74e19a3b2b09f8fa84bf2e59877f025591bfe64fb50fcb8ee6169631ccb58468`
- Records: 158

The complete file contains TruSeq, Nextera, and other sequences. It must not be
passed wholesale to Cutadapt as though all 158 entries were equivalent 3-prime
adapters.

The reviewed `adapters/nextera.fa` and `adapters/truseq.fa` subsets contain the
read-through sequences associated with the named workflow presets. The full
BBMap collection is retained for provenance and custom adapter review.
