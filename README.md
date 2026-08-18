# sv-evidence-extraction

Extracts PE (paired-end), SR (split-read), and RD (read depth) structural
variant evidence from [GATK-SV](https://github.com/broadinstitute/gatk-sv)
`GatherBatchEvidence` outputs (`merged_PE`, `merged_SR`, `merged_bincov`,
`median_cov`), for a given individual and genomic window.

Built to run standalone, from the command line, or as a
[WDL](sv_evidence_extraction.wdl) task in a [Terra](https://terra.bio)
workflow, pulling its inputs from a workspace's DATA/TABLES.

## How it works

PE and SR evidence files are per-batch, bgzipped and tabix-indexed, with
every sample in the batch in one file. RD (bincov) files are similarly
per-batch, but wide-format (one column per sample). This tool:

- Groups requested samples by batch so a batch's evidence file is only
  queried once per region, regardless of how many of its samples (e.g. a
  trio) are involved.
- Opens each batch file exactly once per run via `pysam.TabixFile`
  (falling back to shelling out to the `tabix` CLI if the runtime's
  htslib build lacks GCS support), and reuses that handle across every
  region touching that batch.
- Pads PE/SR queries as small windows around each of an SV's two
  breakpoints (since breakpoint support clusters there, not across the
  whole event), while padding RD queries across the SV's full span
  (since read depth needs the whole interval).

See the docstrings in [`sv_evidence_extraction/core.py`](sv_evidence_extraction/core.py)
for the padding/windowing details.

## Install

```bash
pip install -r requirements.txt
```

## Usage

Two subcommands, sharing the same extraction code:

```bash
# One region, one or a few samples -- an ad hoc spot check.
python -m sv_evidence_extraction.cli query \
    --evidence-paths-tsv evidence_paths.tsv \
    --sample-batch-map-tsv sample_batch_map.tsv \
    --region chr1:161016017-161022143 \
    --sample-ids sample_child,sample_father \
    --out-prefix results/chr1_161016017

# Many independently-scoped regions (e.g. one row per de novo candidate,
# each naming its own child/mother/father) -- the efficient bulk path for
# populating a Terra sample_set row's PE/SR/RD evidence tables.
python -m sv_evidence_extraction.cli build-tables \
    --evidence-paths-tsv evidence_paths.tsv \
    --sample-batch-map-tsv sample_batch_map.tsv \
    --regions-tsv regions.tsv \
    --out-prefix results/my_cohort
```

Both write `<out-prefix>.{pe,sr,rd}.{tsv,parquet}`.

### Inputs

- `--evidence-paths-tsv`: one row per GATK-SV batch, columns
  `entity:sample_set_id`, `median_cov`, `merged_PE`, `merged_SR`,
  `merged_bincov`.
- `--sample-batch-map-tsv`: two columns, no header -- `batch_id`,
  `sample_id`.
- `--regions-tsv` (build-tables mode): columns `chrom`, `start`, `end`,
  `sample_id` (comma-separated for multiple individuals), and optionally
  `name`.

## Running as a Terra workflow

[`sv_evidence_extraction.wdl`](sv_evidence_extraction.wdl) wraps both
subcommands as WDL tasks, switched by a `mode` input (`"query"` or
`"build_tables"`). Build and push the image from [`Dockerfile`](Dockerfile)
to a registry your Terra workspace can pull from, point the WDL's `docker`
input at it, and import the workflow into your workspace (directly, or via
Dockstore for version tracking).

Note: `evidence_paths_tsv`/`sample_batch_map_tsv` are the only WDL `File`
inputs Cromwell localizes -- the `merged_PE`/`merged_SR`/`merged_bincov`
URIs inside `evidence_paths_tsv` stay as plain `gs://` strings read
directly by htslib at runtime, since they're large per-batch files meant
for random access, not bulk download. This means the container needs live
GCS access at runtime, which a Terra/Cromwell GCE worker provides for free
via its attached service account.
