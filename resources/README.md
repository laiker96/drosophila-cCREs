# Pipeline resources

This directory contains reviewed accession and contact-source inputs for the
dm6 atlas plus the adapter sequences used by the workflow.

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

The workflow writes generated final-BAM, master-bundle, and stage-checkpoint
manifests below `results/<project>/<run_id>/provenance/manifests/`. Human-edited
review tables and the accepted/rejected final-BAM manifests derived from them
belong below the ignored `data/reviewed/` namespace. Neither kind belongs in
this versioned resource directory.

Both tables use the default two-ended Tn5/MACS3-qpois ATAC branch. The
with-inputs table deliberately substitutes the matched GSE140539 stage-5
H3K27ac IPs for the IP-only PRJEB37091 stage-5 pair; it is an alternative
scientific selection, not merely the IP-only table plus control rows.
The encoded matched pairs are GSE140539 stage-5 IP/input runs
`SRR10485675/SRR10485676` and `SRR10485677/SRR10485678`, plus adult-brain
`SRR5319052/SRR5319047`.

## Contact sources

`atlas_contact_sources.tsv` is the versioned source manifest for the contact
link stage. It contains the primary GEO download URL, format, replicate label,
target atlas context, biological match description, and caveat for each matrix.
The seven observed contexts are `ab`, `e5`, `e11`, `ead`, `lb`, `o`, and `wid`;
`e13` and `hid` are intentionally absent because the workflow models them by
distance. The table conforms to `schemas/contact-sources.schema.yaml`.

The source repositories do not publish a checksum for every supplementary
matrix, so blank checksum cells are explicit rather than invented. The
workflow checks any checksum present, preserves aria2 resume state, and records
the computed SHA-256 of every downloaded source and normalized output. Raw
matrices are downloaded below the ignored `data/raw/contacts/` namespace and
are never committed.

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
