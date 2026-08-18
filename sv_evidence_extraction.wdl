version 1.0

## Extract PE, SR, and RD structural-variant evidence from GATK-SV batch
## evidence files (merged_PE / merged_SR / merged_bincov / median_cov).
##
## Both modes below share one Python CLI (sv_evidence_extraction.cli):
##   - "query": one ad hoc region + one or more samples -- a spot check.
##   - "build_tables": a regions TSV with many independently-scoped rows
##     (each naming its own sample_id(s)) -- the efficient bulk path for
##     populating a sample_set row's PE/SR/RD evidence-table attributes.
##
## IMPORTANT: evidence_paths_tsv and sample_batch_map_tsv are the only
## WDL File inputs Cromwell localizes. The merged_PE/merged_SR/merged_bincov
## /median_cov URIs *inside* evidence_paths_tsv are plain gs:// strings
## read directly by pysam/htslib at runtime -- they are deliberately NOT
## declared as File inputs, since they're large per-batch files meant to
## be randomly accessed via tabix, not downloaded wholesale. This means
## the container needs live GCS access at runtime; on a Terra/Cromwell
## GCE worker this comes for free from the VM's attached service account
## (confirmed: htslib's built-in GCS backend handles this with no
## GCS_OAUTH_TOKEN/subprocess workaround required).

workflow SVEvidenceExtraction {
  input {
    String mode  # "query" or "build_tables"

    File evidence_paths_tsv
    File sample_batch_map_tsv
    String output_prefix

    Float pad_pct = 0.30
    Int pad_floor = 1000
    Int pad_ceiling_pe_sr = 5000

    # query-mode inputs
    String? region
    String? sample_ids
    String? region_name

    # build_tables-mode input
    File? regions_tsv

    String docker = "us-central1-docker.pkg.dev/PROJECT_ID/sv-evidence-extraction/sv-evidence-extraction:latest"
    Int disk_gb = 20
    Int mem_gb = 8
  }

  if (mode == "query") {
    call QueryEvidence {
      input:
        evidence_paths_tsv    = evidence_paths_tsv,
        sample_batch_map_tsv  = sample_batch_map_tsv,
        region                = select_first([region]),
        sample_ids            = select_first([sample_ids]),
        region_name           = region_name,
        output_prefix         = output_prefix,
        pad_pct               = pad_pct,
        pad_floor             = pad_floor,
        pad_ceiling_pe_sr     = pad_ceiling_pe_sr,
        docker                = docker,
        disk_gb               = disk_gb,
        mem_gb                = mem_gb,
    }
  }

  if (mode == "build_tables") {
    call BuildEvidenceTables {
      input:
        evidence_paths_tsv    = evidence_paths_tsv,
        sample_batch_map_tsv  = sample_batch_map_tsv,
        regions_tsv           = select_first([regions_tsv]),
        output_prefix         = output_prefix,
        pad_pct               = pad_pct,
        pad_floor             = pad_floor,
        pad_ceiling_pe_sr     = pad_ceiling_pe_sr,
        docker                = docker,
        disk_gb               = disk_gb,
        mem_gb                = mem_gb,
    }
  }

  output {
    File? pe_tsv      = select_first([QueryEvidence.pe_tsv, BuildEvidenceTables.pe_tsv])
    File? pe_parquet  = select_first([QueryEvidence.pe_parquet, BuildEvidenceTables.pe_parquet])
    File? sr_tsv      = select_first([QueryEvidence.sr_tsv, BuildEvidenceTables.sr_tsv])
    File? sr_parquet  = select_first([QueryEvidence.sr_parquet, BuildEvidenceTables.sr_parquet])
    File? rd_tsv      = select_first([QueryEvidence.rd_tsv, BuildEvidenceTables.rd_tsv])
    File? rd_parquet  = select_first([QueryEvidence.rd_parquet, BuildEvidenceTables.rd_parquet])
  }
}

task QueryEvidence {
  input {
    File evidence_paths_tsv
    File sample_batch_map_tsv
    String region
    String sample_ids
    String? region_name
    String output_prefix
    Float pad_pct
    Int pad_floor
    Int pad_ceiling_pe_sr
    String docker
    Int disk_gb
    Int mem_gb
  }

  command <<<
    set -euo pipefail
    python3 -m sv_evidence_extraction.cli query \
      --evidence-paths-tsv ~{evidence_paths_tsv} \
      --sample-batch-map-tsv ~{sample_batch_map_tsv} \
      --region ~{region} \
      --sample-ids ~{sample_ids} \
      ~{"--name " + region_name} \
      --pad-pct ~{pad_pct} \
      --pad-floor ~{pad_floor} \
      --pad-ceiling-pe-sr ~{pad_ceiling_pe_sr} \
      --out-prefix ~{output_prefix}
  >>>

  output {
    File pe_tsv     = "~{output_prefix}.pe.tsv"
    File pe_parquet = "~{output_prefix}.pe.parquet"
    File sr_tsv     = "~{output_prefix}.sr.tsv"
    File sr_parquet = "~{output_prefix}.sr.parquet"
    File rd_tsv     = "~{output_prefix}.rd.tsv"
    File rd_parquet = "~{output_prefix}.rd.parquet"
  }

  runtime {
    docker: docker
    memory: "~{mem_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
    cpu: 2
  }
}

task BuildEvidenceTables {
  input {
    File evidence_paths_tsv
    File sample_batch_map_tsv
    File regions_tsv
    String output_prefix
    Float pad_pct
    Int pad_floor
    Int pad_ceiling_pe_sr
    String docker
    Int disk_gb
    Int mem_gb
  }

  command <<<
    set -euo pipefail
    python3 -m sv_evidence_extraction.cli build-tables \
      --evidence-paths-tsv ~{evidence_paths_tsv} \
      --sample-batch-map-tsv ~{sample_batch_map_tsv} \
      --regions-tsv ~{regions_tsv} \
      --pad-pct ~{pad_pct} \
      --pad-floor ~{pad_floor} \
      --pad-ceiling-pe-sr ~{pad_ceiling_pe_sr} \
      --out-prefix ~{output_prefix}
  >>>

  output {
    File pe_tsv     = "~{output_prefix}.pe.tsv"
    File pe_parquet = "~{output_prefix}.pe.parquet"
    File sr_tsv     = "~{output_prefix}.sr.tsv"
    File sr_parquet = "~{output_prefix}.sr.parquet"
    File rd_tsv     = "~{output_prefix}.rd.tsv"
    File rd_parquet = "~{output_prefix}.rd.parquet"
  }

  runtime {
    docker: docker
    memory: "~{mem_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
    cpu: 4
  }
}
