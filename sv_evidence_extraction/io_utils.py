"""I/O helpers for the Terra-table inputs and dual-format (TSV + Parquet) outputs.

Kept separate from core.py so the extraction logic in core.py has no
opinions about where inputs come from or how outputs are laid out on
disk -- useful for testing core.py against synthetic fixtures without
touching real Terra tables.
"""
import pandas as pd

from .core import EvidenceIndex


def load_evidence_index(evidence_paths_tsv, sample_batch_map_tsv):
    """Build an `EvidenceIndex` from the two Terra-table TSVs this tool needs.

    Parameters
    ----------
    evidence_paths_tsv : str
        Path (local or gs://) to the sample_set-level evidence-paths
        table -- one row per batch, columns entity:sample_set_id,
        median_cov, merged_PE, merged_SR, merged_bincov.
    sample_batch_map_tsv : str
        Path to the sample -> batch map: two columns, no header
        (batch_id, sample_id).

    Returns
    -------
    EvidenceIndex
    """
    df_evidence = pd.read_csv(evidence_paths_tsv, sep="\t")
    df_batch = pd.read_csv(sample_batch_map_tsv, sep="\t", names=["batch_id", "sample_id"])
    return EvidenceIndex.from_tables(df_evidence, df_batch)


def write_evidence_tables(tables, out_prefix):
    """Write each evidence DataFrame to both TSV and Parquet.

    Files are named `{out_prefix}.{evidence_class}.{tsv,parquet}`, e.g.
    `my_run.pe.tsv` / `my_run.pe.parquet`. TSV is convenient for
    downstream cluster tools that expect plain text; Parquet is far
    cheaper to pull down and load locally, especially for the RD table.

    Parameters
    ----------
    tables : dict of str -> pandas.DataFrame
        As returned by `core.build_evidence_tables` (keys "pe", "sr", "rd").
    out_prefix : str
        Local path prefix for output files.

    Returns
    -------
    dict of str -> dict of str -> str
        {evidence_class: {"tsv": path, "parquet": path}}
    """
    written = {}
    for evidence_class, df in tables.items():
        tsv_path = f"{out_prefix}.{evidence_class}.tsv"
        parquet_path = f"{out_prefix}.{evidence_class}.parquet"
        df.to_csv(tsv_path, sep="\t", index=False)
        df.to_parquet(parquet_path, index=False)
        written[evidence_class] = {"tsv": tsv_path, "parquet": parquet_path}
    return written
