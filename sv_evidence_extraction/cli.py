"""Command-line entry point for extracting PE/SR/RD structural variant evidence.

Two subcommands, sharing the same underlying extraction code:

query
    One region, one or a few samples. For ad hoc/interactive spot
    checks -- e.g. the WDL-task equivalent of "pull evidence for this
    one candidate call".

build-tables
    A regions TSV with many rows, each independently scoped to its own
    sample_id(s) (e.g. one row per de novo candidate, each naming that
    candidate's child/mother/father). This is the efficient path for
    building the three Terra sample_set-level evidence tables for a
    few hundred individuals at once, since samples are grouped by
    GATK-SV batch internally so a batch's evidence file is only opened
    once no matter how many regions touch it.

Example
-------
    python -m sv_evidence_extraction.cli query \\
        --evidence-paths-tsv evidence_paths.tsv \\
        --sample-batch-map-tsv sample_batch_map.tsv \\
        --ped-file pedigree.ped \\
        --region chr1:161016017-161022143 \\
        --sample-ids sample_child,sample_father \\
        --out-prefix results/chr1_161016017

    python -m sv_evidence_extraction.cli build-tables \\
        --evidence-paths-tsv evidence_paths.tsv \\
        --sample-batch-map-tsv sample_batch_map.tsv \\
        --regions-tsv denovo_candidates.tsv \\
        --out-prefix results/asd_cohort_denovo
"""
import argparse
import sys

from .core import RegionRequest, build_evidence_tables, load_regions_table
from .io_utils import load_evidence_index, write_evidence_tables
from .pedigree import load_pedigree


def _add_common_args(parser):
    parser.add_argument(
        "--evidence-paths-tsv", required=True,
        help="Terra sample_set evidence-paths table "
             "(entity:sample_set_id, median_cov, merged_PE, merged_SR, merged_bincov).",
    )
    parser.add_argument(
        "--sample-batch-map-tsv", required=True,
        help="Two-column (batch_id, sample_id) sample->batch map, no header.",
    )
    parser.add_argument(
        "--ped-file", default=None,
        help="Optional 6-column PED file. When given, adds a \"relationship\" "
             "(child/father/mother/unknown) column to every output table, labeled "
             "per region from that region's own sample_ids.",
    )
    parser.add_argument(
        "--out-prefix", required=True,
        help="Output path prefix; writes <prefix>.{pe,sr,rd}.{tsv,parquet}.",
    )
    parser.add_argument(
        "--pad-pct", type=float, default=0.30,
        help="Fraction of event size to pad on each side (default 0.30).",
    )
    parser.add_argument(
        "--pad-floor", type=int, default=1000,
        help="Minimum padding in bp, for both RD's full-span window and PE/SR's "
             "per-breakpoint windows (default 1000).",
    )
    parser.add_argument(
        "--pad-ceiling-pe-sr", type=int, default=5000,
        help="Maximum per-breakpoint padding in bp for PE/SR windows (default 5000). "
             "RD's full-span padding has no ceiling, since it needs to cover the whole event.",
    )


def _run(regions, args):
    evidence_index = load_evidence_index(args.evidence_paths_tsv, args.sample_batch_map_tsv)
    df_ped = load_pedigree(args.ped_file) if args.ped_file else None
    tables = build_evidence_tables(
        regions, evidence_index,
        pad_pct=args.pad_pct, pad_floor=args.pad_floor, pad_ceiling=args.pad_ceiling_pe_sr,
        df_ped=df_ped,
    )
    written = write_evidence_tables(tables, args.out_prefix)
    for evidence_class, paths in written.items():
        print(f"{evidence_class.upper()}: {len(tables[evidence_class])} rows -> "
              f"{paths['tsv']}, {paths['parquet']}", file=sys.stderr)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    p_query = subparsers.add_parser("query", help="Extract evidence for one region and one or more samples.")
    _add_common_args(p_query)
    p_query.add_argument("--region", required=True, help="chrom:start-end, 1-based inclusive.")
    p_query.add_argument("--sample-ids", required=True, help="Comma-separated sample IDs, e.g. child,mother,father.")
    p_query.add_argument("--name", default=None, help="Label for this region (default: derived from the region string).")

    p_build = subparsers.add_parser("build-tables", help="Extract evidence for many independently-scoped regions.")
    _add_common_args(p_build)
    p_build.add_argument(
        "--regions-tsv", required=True,
        help="TSV with columns chrom, start, end, sample_id[, name]. "
             "sample_id may be a comma-separated list per row.",
    )

    args = parser.parse_args(argv)

    if args.command == "query":
        chrom, coords = args.region.split(":")
        start_str, end_str = coords.split("-")
        start, end = int(start_str), int(end_str)
        sample_ids = [s.strip() for s in args.sample_ids.split(",") if s.strip()]
        name = args.name or f"{chrom}_{start}_{end}"
        regions = [RegionRequest(name=name, chrom=chrom, start=start, end=end, sample_ids=sample_ids)]
    else:
        regions = load_regions_table(args.regions_tsv)

    _run(regions, args)


if __name__ == "__main__":
    main()
