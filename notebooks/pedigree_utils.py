"""Pedigree helpers shared by the local analysis notebooks.

Lives outside sv_evidence_extraction/ deliberately: this is local
analysis tooling for building inputs and labeling results, not part of
the extraction pipeline that ships in the Docker image. Import from
here in any notebook that needs pedigree lookups, rather than pasting a
fresh copy -- see sv_review_utils.py in the sibling ASC_denovo_SV repo
for what happens when that discipline slips (six drifted copies of the
same function across six notebooks).
"""
import pandas as pd

PED_HEADER = ["FamID", "IndividualID", "FatherID", "MotherID", "Gender", "Affected"]

# PLINK-style PED files use "0" (a literal string, not a blank/NaN cell)
# as the sentinel for "no parent recorded" -- confirmed against this
# cohort's actual pedigree file, which has zero blank FatherID/MotherID
# cells but ~23k rows with FatherID == "0". A plain pd.notna() check
# would silently treat "0" as a real sample ID and try to pull evidence
# for a nonexistent sample "0".
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


def resolve_family_sample_ids(child_id, df_ped, include="both"):
    """Resolve a child's sample_ids list (itself plus requested parent(s)) from a PED table.

    Parameters
    ----------
    child_id : str
        IndividualID of the child, as it appears in df_ped["IndividualID"].
    df_ped : pandas.DataFrame
        As returned by `load_pedigree`.
    include : {"both", "father", "mother", "none"}, default "both"
        Which parent(s) to include alongside the child. "none" returns
        just the child (e.g. for a call where only the proband's own
        evidence is of interest).

    Returns
    -------
    list of str
        [child_id] plus any resolved, non-missing parent IDs, in that order.

    Raises
    ------
    ValueError
        If child_id isn't found in df_ped, since a silently-empty result
        would otherwise look like "no evidence found" downstream rather
        than "this sample isn't even in the pedigree".
    """
    matches = df_ped.loc[df_ped["IndividualID"] == child_id]
    if matches.empty:
        raise ValueError(f"'{child_id}' not found in the pedigree table.")
    row = matches.iloc[0]

    sample_ids = [child_id]
    if include in ("both", "father") and _has_parent(row["FatherID"]):
        sample_ids.append(row["FatherID"])
    if include in ("both", "mother") and _has_parent(row["MotherID"]):
        sample_ids.append(row["MotherID"])
    return sample_ids


def label_family_roles(sample_ids, df_ped):
    """Label each sample_id in a family group as "child", "father", or "mother".

    Roles are inferred from the pedigree relationships *among the given
    sample_ids themselves*, rather than assuming a fixed input order:
    whichever sample appears as another's FatherID is "father", whichever
    appears as another's MotherID is "mother", and a sample whose own
    FatherID/MotherID resolves to one of the others is "child". This
    works regardless of which parent(s) `resolve_family_sample_ids` was
    asked to include, and doesn't depend on any naming convention in the
    evidence tables themselves.

    Parameters
    ----------
    sample_ids : iterable of str
        The sample_ids present in one result set (e.g. from an evidence
        table's "sample_id" column).
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


def label_evidence_relationships(evidence, df_ped):
    """Add a "relationship" column (child/father/mother/unknown) to every table in an evidence dict.

    Parameters
    ----------
    evidence : dict of str -> pandas.DataFrame
        As returned by `examine_evidence.ipynb`'s `load_evidence_set` --
        each DataFrame must have a "sample_id" column.
    df_ped : pandas.DataFrame
        As returned by `load_pedigree`.

    Returns
    -------
    dict of str -> pandas.DataFrame
        Same keys, each with a "relationship" column inserted right
        after "sample_id".
    """
    all_sample_ids = set()
    for df in evidence.values():
        if not df.empty:
            all_sample_ids.update(df["sample_id"].unique())

    sample_roles = label_family_roles(all_sample_ids, df_ped)

    labeled = {}
    for evidence_class, df in evidence.items():
        out = df.copy()
        if out.empty:
            out["relationship"] = pd.Series(dtype=str)
        else:
            out.insert(
                out.columns.get_loc("sample_id") + 1,
                "relationship",
                out["sample_id"].map(sample_roles).fillna("unknown"),
            )
        labeled[evidence_class] = out
    return labeled
