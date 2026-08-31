"""Pedigree-based family-role labeling, shared by the pipeline and local notebooks.

This lives in the shipped package (unlike notebooks/pedigree_utils.py's
input-building helpers) because `core.build_evidence_tables` applies it
at extraction time when a pedigree is supplied, so a colleague can be
handed an output table with the relationship already built in rather
than needing to run a local notebook step to get it.
"""
import pandas as pd

PED_HEADER = ["FamID", "IndividualID", "FatherID", "MotherID", "Gender", "Affected"]

# PLINK-style PED files use "0" (a literal string, not a blank/NaN cell)
# as the sentinel for "no parent recorded" -- confirmed against a real
# cohort's pedigree file, which had zero blank FatherID/MotherID cells
# but ~23k rows with FatherID == "0". A plain pd.notna() check would
# silently treat "0" as a real sample ID.
MISSING_PARENT_SENTINELS = {"0", 0}


def load_pedigree(ped_file_uri):
    """Load a standard 6-column PED file into a DataFrame.

    Parameters
    ----------
    ped_file_uri : str
        Local path or gs:// URI to the PED file.

    Returns
    -------
    pandas.DataFrame
        Columns: FamID, IndividualID, FatherID, MotherID, Gender, Affected.
    """
    return pd.read_csv(ped_file_uri, sep="\t", names=PED_HEADER, comment="#")


def _has_parent(value):
    """True if `value` is a real parent ID, not missing/the "0" sentinel."""
    return pd.notna(value) and value not in MISSING_PARENT_SENTINELS


def label_family_roles(sample_ids, df_ped):
    """Label each sample_id in a family group as "child", "father", or "mother".

    Roles are inferred from the pedigree relationships *among the given
    sample_ids themselves*, rather than assuming a fixed input order:
    whichever sample appears as another's FatherID is "father", whichever
    appears as another's MotherID is "mother", and a sample whose own
    FatherID/MotherID resolves to one of the others is "child".

    Deliberately scoped to just the sample_ids passed in, not the whole
    pedigree: a build-tables run touches many unrelated families in one
    process, and a sample's role can genuinely differ by context (e.g. a
    parent in one family who is a childless founder in the query for a
    different one). Callers must call this once per family group -- see
    `core.build_evidence_tables`, which does so per region -- rather than
    on the union of sample_ids across an entire run.

    Parameters
    ----------
    sample_ids : iterable of str
        The sample_ids present in one family group.
    df_ped : pandas.DataFrame
        As returned by `load_pedigree`.

    Returns
    -------
    dict of str -> str
        sample_id -> one of "child", "father", "mother", "unknown" (for
        a sample with no resolvable relationship to the others, e.g. a
        lone founder or an unrelated sample in the mix).
    """
    sample_ids = set(sample_ids)
    ped_rows = df_ped.loc[df_ped["IndividualID"].isin(sample_ids)].set_index("IndividualID")

    roles = {sample_id: "unknown" for sample_id in sample_ids}
    for sample_id in sample_ids:
        if sample_id not in ped_rows.index:
            continue
        father_id = ped_rows.loc[sample_id, "FatherID"]
        mother_id = ped_rows.loc[sample_id, "MotherID"]
        if father_id in sample_ids:
            roles[sample_id] = "child"
            roles[father_id] = "father"
        if mother_id in sample_ids:
            roles[sample_id] = "child"
            roles[mother_id] = "mother"
    return roles


def add_relationship_column(df, role_map):
    """Insert a "relationship" column right after "sample_id", mapped from `role_map`.

    Parameters
    ----------
    df : pandas.DataFrame
        Must have a "sample_id" column; may be empty.
    role_map : dict of str -> str
        sample_id -> role, as returned by `label_family_roles`.

    Returns
    -------
    pandas.DataFrame
        Copy of `df` with the new column. Samples not present in
        `role_map` get "unknown" rather than NaN.
    """
    out = df.copy()
    if out.empty:
        out["relationship"] = pd.Series(dtype=str)
        return out
    out.insert(
        out.columns.get_loc("sample_id") + 1,
        "relationship",
        out["sample_id"].map(role_map).fillna("unknown"),
    )
    return out
