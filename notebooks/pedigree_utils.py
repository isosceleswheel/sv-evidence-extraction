"""Pedigree helpers for the local analysis notebooks.

`load_pedigree`, `label_family_roles`, and `add_relationship_column` now
live in sv_evidence_extraction/pedigree.py (part of the shipped package,
since core.build_evidence_tables applies them at extraction time too) --
re-exported here so notebook imports don't need to change. This module
keeps only the input-building/local-analysis helpers that have no reason
to ship in the Docker image: `resolve_family_sample_ids` and
`label_evidence_relationships`.

Import from here in any notebook that needs pedigree lookups, rather
than pasting a fresh copy -- see sv_review_utils.py in the sibling
ASC_denovo_SV repo for what happens when that discipline slips (six
drifted copies of the same function across six notebooks).
"""
import sys
from pathlib import Path

# sv_evidence_extraction lives at the repo root, one level up from
# notebooks/ -- add it via this file's own location rather than the
# caller's cwd, so this works regardless of where a notebook is launched from.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from sv_evidence_extraction.pedigree import (  # noqa: E402
    PED_HEADER,
    MISSING_PARENT_SENTINELS,
    load_pedigree,
    _has_parent,
    label_family_roles,
    add_relationship_column,
)


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


def label_evidence_relationships(evidence, df_ped):
    """Add a "relationship" column (child/father/mother/unknown) to every table in an evidence dict.

    Unlike `core.build_evidence_tables` (which labels each region's own
    family group separately, since a build-tables run can mix many
    unrelated families), this labels from the union of sample_ids across
    all three already-loaded tables -- fine here since `load_evidence_set`
    loads one region's result set at a time.

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

    role_map = label_family_roles(all_sample_ids, df_ped)
    return {evidence_class: add_relationship_column(df, role_map) for evidence_class, df in evidence.items()}
