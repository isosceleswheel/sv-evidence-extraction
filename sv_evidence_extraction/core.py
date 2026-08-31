"""Core evidence-extraction logic for PE, SR, and RD structural variant evidence.

This module has no Terra/GCS-specific I/O in it beyond the ability to open
a tabix-indexed file by URL (local path or gs://) -- see io_utils.py for
loading the Terra evidence-paths / sample-batch-map tables that feed
`EvidenceIndex`, and cli.py for the command-line wrapper.

Design notes
------------
- PE and SR evidence files are per-batch (not per-sample): one bgzipped,
  tabix-indexed file holds every sample in a GATK-SV batch, with the
  sample ID as the last column. RD (bincov) files are similarly
  per-batch, but wide-format (one column per sample). Because of this,
  every extraction here is keyed by *batch*, not by sample: samples
  requested for a region are grouped by batch first so a batch's file
  is only queried once per region regardless of how many of its samples
  are involved (e.g. a trio that rebatched together).
- Opening a remote tabix file has real latency (index fetch + auth
  handshake). `TabixHandleCache` keeps one open handle per URL for the
  life of a whole run, so a batch file touched by many regions across a
  build-tables run is only opened once total, not once per region.
- PE/SR support a variant's two breakpoints, not its interior. Padding
  the *entire* span of a multi-Mb SV would pull a large amount of
  irrelevant evidence from the middle and cost a much bigger tabix scan
  for no benefit, so `breakpoint_windows()` returns small windows around
  each breakpoint (merged into one if the event is small enough that
  they'd overlap). RD genuinely needs the full padded span -- that's the
  depth signal across the whole event -- so `pad_window()` is used for it
  instead.
"""
import os
import sys
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

try:
    import pysam
    _HAVE_PYSAM = True
except ImportError:  # pragma: no cover - exercised only where pysam truly isn't installed
    pysam = None
    _HAVE_PYSAM = False

import subprocess


def _refresh_gcs_oauth_token():
    """Set GCS_OAUTH_TOKEN to a fresh access token, for htslib's GCS backend.

    htslib's built-in GCS support (used by both pysam and the `tabix` CLI)
    does not discover Application Default Credentials on its own the way
    Python's google-auth/google-cloud-storage libraries do -- confirmed by
    hand against a real private bucket, where pysam.TabixFile() failed to
    open a file the same ADC credentials could read fine via
    google.cloud.storage. It only works if GCS_OAUTH_TOKEN is set
    explicitly, so this shells out to `gcloud auth print-access-token`
    (matching the token-minting approach already proven in the original
    analysis notebook this tool replaced) before every new file is opened.
    """
    result = subprocess.run(
        ["gcloud", "auth", "print-access-token"], capture_output=True, text=True,
    )
    if result.returncode != 0:
        print(
            f"WARNING: could not obtain a GCS access token ({result.stderr.strip()}); "
            "access to private gs:// evidence files will likely fail.",
            file=sys.stderr,
        )
        return
    os.environ["GCS_OAUTH_TOKEN"] = result.stdout.strip()


# ======================================================================
# WINDOWING
# ======================================================================

def pad_window(start, end, pad_pct=0.30, pad_floor=1000):
    """Symmetrically pad a [start, end] interval for full-span evidence queries (RD).

    Parameters
    ----------
    start, end : int
        1-based inclusive event coordinates.
    pad_pct : float, default 0.30
        Fraction of event size to pad on each side.
    pad_floor : int, default 1000
        Minimum padding in bp -- without this, tiny SVs (e.g. a 50bp
        insertion) would get almost no padding at all.

    Returns
    -------
    (int, int)
        Padded (start, end); start is clipped to a minimum of 1.
    """
    event_size = max(end - start, 0)
    pad = max(int(event_size * pad_pct), pad_floor)
    return max(1, start - pad), end + pad


def breakpoint_windows(start, end, pad_pct=0.30, pad_floor=1000, pad_ceiling=5000):
    """Compute one or two padded windows around an SV's breakpoints, for PE/SR queries.

    The per-breakpoint pad is capped at `pad_ceiling` regardless of event
    size: breakpoint-supporting reads cluster tightly around the true
    breakpoint, so a multi-Mb SV doesn't warrant a multi-Mb pad the way it
    would for `pad_window`'s full-span RD query.

    If the two padded breakpoint windows overlap (small event, or pad
    large relative to event size), they're merged into a single window
    rather than returned as two overlapping queries against the same
    tabix file.

    Parameters
    ----------
    start, end : int
        1-based inclusive event coordinates.
    pad_pct : float, default 0.30
        Fraction of event size used to compute the pad, before clamping.
    pad_floor : int, default 1000
        Minimum pad in bp.
    pad_ceiling : int, default 5000
        Maximum pad in bp, regardless of event size.

    Returns
    -------
    list of (int, int)
        A single (start, end) tuple if the breakpoint windows merge,
        otherwise two tuples -- the start-breakpoint window followed by
        the end-breakpoint window.
    """
    event_size = max(end - start, 0)
    pad = max(pad_floor, min(int(event_size * pad_pct), pad_ceiling))

    left_window = (max(1, start - pad), start + pad)
    right_window = (max(1, end - pad), end + pad)

    if left_window[1] >= right_window[0]:
        return [(left_window[0], right_window[1])]
    return [left_window, right_window]


# ======================================================================
# SAMPLE / BATCH RESOLUTION
# ======================================================================

@dataclass
class EvidenceIndex:
    """Lookup tables mapping sample IDs to their per-batch evidence file URLs.

    Built once per run from a Terra evidence-paths table (one row per
    GATK-SV batch, with PE/SR/RD/median_cov URIs) and a sample -> batch
    map, then reused for every region in that run.

    Attributes
    ----------
    sample_to_batch : dict
        sample_id -> batch_id.
    pe_url, sr_url, rd_url, cov_url : dict
        batch_id -> GCS/local URI of that batch's merged_PE, merged_SR,
        merged_bincov, and median_cov files respectively.
    """
    sample_to_batch: dict
    pe_url: dict
    sr_url: dict
    rd_url: dict
    cov_url: dict

    @classmethod
    def from_tables(cls, df_evidence, df_batch,
                     batch_entity_col="entity:sample_set_id",
                     pe_col="merged_PE", sr_col="merged_SR",
                     rd_col="merged_bincov", cov_col="median_cov"):
        """Build an EvidenceIndex from the evidence-paths and sample-batch-map DataFrames.

        Parameters
        ----------
        df_evidence : pandas.DataFrame
            One row per batch, as loaded from the Terra evidence-paths
            table (e.g. evidence_paths.tsv).
        df_batch : pandas.DataFrame
            Two columns: "batch_id" and "sample_id".
        batch_entity_col, pe_col, sr_col, rd_col, cov_col : str
            Column names in `df_evidence`, overridable in case the Terra
            table schema drifts from the current convention.
        """
        sample_to_batch = dict(zip(df_batch["sample_id"], df_batch["batch_id"]))
        return cls(
            sample_to_batch=sample_to_batch,
            pe_url=dict(zip(df_evidence[batch_entity_col], df_evidence[pe_col])),
            sr_url=dict(zip(df_evidence[batch_entity_col], df_evidence[sr_col])),
            rd_url=dict(zip(df_evidence[batch_entity_col], df_evidence[rd_col])),
            cov_url=dict(zip(df_evidence[batch_entity_col], df_evidence[cov_col])),
        )

    def group_by_batch(self, sample_ids):
        """Group sample IDs by their batch, dropping (with a warning) any with no known batch.

        Returns
        -------
        dict
            batch_id -> list of sample_ids from `sample_ids` in that batch.
        """
        groups = {}
        for sample_id in sample_ids:
            batch_id = self.sample_to_batch.get(sample_id)
            if batch_id is None:
                print(f"WARNING: no batch found for sample '{sample_id}'; skipping.", file=sys.stderr)
                continue
            groups.setdefault(batch_id, []).append(sample_id)
        return groups


# ======================================================================
# TABIX ACCESS (pysam-native, with subprocess fallback)
# ======================================================================

class TabixSource:
    """A cached, reusable handle to one tabix-indexed evidence file.

    Wraps `pysam.TabixFile` when it can open the URL directly (opened
    once, then reused for many `.fetch()` calls against different
    regions -- this is what makes batch queries across many regions
    cheap). Before opening a gs:// URL, refreshes GCS_OAUTH_TOKEN, since
    htslib's GCS backend needs it explicitly rather than discovering
    Application Default Credentials on its own. Falls back to shelling
    out to the `tabix` CLI per-query if pysam/htslib still can't open
    the URL, so the tool still works, just at the notebook's original
    speed.
    """

    def __init__(self, url):
        self.url = url
        self._handle = None
        self._use_subprocess = False
        self._open()

    def _open(self):
        if self.url.startswith("gs://"):
            _refresh_gcs_oauth_token()

        if _HAVE_PYSAM:
            try:
                self._handle = pysam.TabixFile(self.url)
                return
            except Exception as exc:
                print(
                    f"NOTE: pysam could not open '{self.url}' directly ({exc}); "
                    "falling back to `tabix` subprocess calls for this file.",
                    file=sys.stderr,
                )
        self._use_subprocess = True

    def fetch_lines(self, chrom, start, end):
        """Yield raw tab-delimited data lines overlapping chrom:start-end (1-based inclusive)."""
        if not self._use_subprocess:
            try:
                # pysam fetch uses 0-based half-open coordinates regardless
                # of the underlying file's convention.
                yield from self._handle.fetch(chrom, max(start - 1, 0), end)
                return
            except Exception as exc:
                print(
                    f"NOTE: pysam fetch failed on '{self.url}' ({exc}); "
                    "falling back to `tabix` subprocess for this file.",
                    file=sys.stderr,
                )
                self._use_subprocess = True

        region = f"{chrom}:{start}-{end}"
        cmd = f'GCS_OAUTH_TOKEN="$(gcloud auth print-access-token)" tabix "{self.url}" "{region}"'
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"WARNING: tabix subprocess failed for {self.url} {region}: {result.stderr.strip()}", file=sys.stderr)
            return
        for line in result.stdout.splitlines():
            if line:
                yield line

    def header_line(self):
        """Return the tab-delimited header line (e.g. bincov matrix column names), or None if there isn't one."""
        if not self._use_subprocess:
            try:
                header_lines = list(self._handle.header)
                if header_lines:
                    return header_lines[-1].lstrip("#")
                return None
            except Exception:
                pass

        cmd = f'GCS_OAUTH_TOKEN="$(gcloud auth print-access-token)" tabix -H "{self.url}"'
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in result.stdout.splitlines() if l]
        return lines[-1].lstrip("#") if lines else None

    def close(self):
        if self._handle is not None and hasattr(self._handle, "close"):
            self._handle.close()


class TabixHandleCache:
    """Keeps one open TabixSource per URL for the lifetime of a run.

    Passed through a whole `build_evidence_tables` run so a given
    batch's PE/SR/RD file is opened exactly once, no matter how many
    regions in the run need it.
    """

    def __init__(self):
        self._sources = {}

    def get(self, url):
        if url not in self._sources:
            self._sources[url] = TabixSource(url)
        return self._sources[url]

    def close_all(self):
        for source in self._sources.values():
            source.close()
        self._sources.clear()


# ======================================================================
# EVIDENCE EXTRACTION
# ======================================================================

PE_COLUMNS = ["chrom1", "pos1", "strand1", "chrom2", "pos2", "strand2", "sample_id"]
SR_COLUMNS = ["chrom", "pos", "orientation", "count", "sample_id"]
RD_COLUMNS = ["chrom", "start", "end", "sample_id", "read_depth", "median_cov"]


def extract_pe(handle_cache, url, windows, sample_ids):
    """Extract PE (discordant paired-end) evidence for a set of samples from one batch's PE file.

    Parameters
    ----------
    handle_cache : TabixHandleCache
    url : str
        URI of this batch's merged_PE file.
    windows : list of (chrom, start, end)
        One or two windows to query -- typically the output of
        `breakpoint_windows()`.
    sample_ids : list of str
        Only rows whose sample_id is in this set are kept (the file
        contains every sample in the batch).

    Returns
    -------
    pandas.DataFrame
        Columns: PE_COLUMNS, deduplicated across overlapping windows.
    """
    source = handle_cache.get(url)
    sample_set = set(sample_ids)
    rows = []
    for chrom, start, end in windows:
        for line in source.fetch_lines(chrom, start, end):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != len(PE_COLUMNS):
                continue
            if fields[-1] not in sample_set:
                continue
            rows.append(fields)

    df = pd.DataFrame(rows, columns=PE_COLUMNS)
    if not df.empty:
        df[["pos1", "pos2"]] = df[["pos1", "pos2"]].astype(int)
        df = df.drop_duplicates().reset_index(drop=True)
    return df


def extract_sr(handle_cache, url, windows, sample_ids):
    """Extract SR (split-read) evidence for a set of samples. Same shape as `extract_pe`."""
    source = handle_cache.get(url)
    sample_set = set(sample_ids)
    rows = []
    for chrom, start, end in windows:
        for line in source.fetch_lines(chrom, start, end):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != len(SR_COLUMNS):
                continue
            if fields[-1] not in sample_set:
                continue
            rows.append(fields)

    df = pd.DataFrame(rows, columns=SR_COLUMNS)
    if not df.empty:
        df[["pos", "count"]] = df[["pos", "count"]].astype(int)
        df = df.drop_duplicates().reset_index(drop=True)
    return df


def _load_median_cov(cov_url, cov_cache):
    """Load (and cache) a batch's whole median_cov table -- small enough to read in full, not tabix-queried."""
    if cov_cache is not None and cov_url in cov_cache:
        return cov_cache[cov_url]
    df = pd.read_csv(cov_url, sep="\t")
    if cov_cache is not None:
        cov_cache[cov_url] = df
    return df


def extract_rd(handle_cache, rd_url, window, sample_ids, cov_url=None, cov_cache=None):
    """Extract RD (bincov) evidence for a set of samples over one padded window, long-format, with median_cov merged in.

    Parameters
    ----------
    handle_cache : TabixHandleCache
    rd_url : str
        URI of this batch's merged_bincov file.
    window : (chrom, start, end)
        The full padded event span (unlike PE/SR, RD needs the whole
        interval -- that's the depth signal across the event).
    sample_ids : list of str
    cov_url : str, optional
        URI of this batch's median_cov file, for the normalization column.
    cov_cache : dict, optional
        Cache of already-loaded median_cov DataFrames keyed by cov_url,
        so a batch's (small) median_cov file isn't reloaded for every
        region that touches that batch.

    Returns
    -------
    pandas.DataFrame
        Long-format columns: RD_COLUMNS.
    """
    chrom, start, end = window
    source = handle_cache.get(rd_url)
    header = source.header_line()
    if header is None:
        return pd.DataFrame(columns=RD_COLUMNS)

    columns = header.split("\t")
    rows = [line.split("\t") for line in source.fetch_lines(chrom, start, end)]
    if not rows:
        return pd.DataFrame(columns=RD_COLUMNS)

    df = pd.DataFrame(rows, columns=columns)
    chrom_col, start_col, end_col = columns[:3]
    valid_samples = [s for s in sample_ids if s in df.columns]
    if not valid_samples:
        return pd.DataFrame(columns=RD_COLUMNS)

    df = df[[chrom_col, start_col, end_col] + valid_samples].copy()
    df[valid_samples] = df[valid_samples].astype(float)
    df[[start_col, end_col]] = df[[start_col, end_col]].astype(int)

    df_long = df.melt(
        id_vars=[chrom_col, start_col, end_col],
        value_vars=valid_samples,
        var_name="sample_id",
        value_name="read_depth",
    ).rename(columns={chrom_col: "chrom", start_col: "start", end_col: "end"})

    if cov_url:
        df_cov = _load_median_cov(cov_url, cov_cache)
        cov_lookup = df_cov.iloc[0].to_dict() if not df_cov.empty else {}
        df_long["median_cov"] = df_long["sample_id"].map(cov_lookup)
    else:
        df_long["median_cov"] = np.nan

    return df_long[RD_COLUMNS]


# ======================================================================
# REGION REQUESTS + TOP-LEVEL ORCHESTRATION
# ======================================================================

@dataclass
class RegionRequest:
    """One event needing PE/SR/RD evidence pulled for a specific set of samples.

    Attributes
    ----------
    name : str
        Identifier for this event, carried through to every output row
        so results can be joined back to the source call.
    chrom : str
    start, end : int
        1-based inclusive core event coordinates (padding is applied
        internally by `build_evidence_tables`, not here).
    sample_ids : list of str
        The individuals this event's evidence should be pulled for
        (e.g. a de novo candidate's child, mother, and father).
    """
    name: str
    chrom: str
    start: int
    end: int
    sample_ids: list = field(default_factory=list)


def load_regions_table(path):
    """Load a regions TSV into a list of `RegionRequest`.

    Expected columns
    -----------------
    chrom, start, end : genomic coordinates, 1-based inclusive.
    sample_id : one sample ID, or several comma-separated (e.g. a de
        novo candidate's "child,mother,father").
    name : optional; auto-generated from chrom/start/end if absent.
    """
    df = pd.read_csv(path, sep="\t")
    required = {"chrom", "start", "end", "sample_id"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Regions table is missing required columns: {sorted(missing)}")

    requests = []
    for _, row in df.iterrows():
        sample_ids = [s.strip() for s in str(row["sample_id"]).split(",") if s.strip()]
        if "name" in df.columns and pd.notna(row.get("name")):
            name = row["name"]
        else:
            name = f"{row['chrom']}_{row['start']}_{row['end']}"
        requests.append(RegionRequest(
            name=name, chrom=row["chrom"], start=int(row["start"]), end=int(row["end"]),
            sample_ids=sample_ids,
        ))
    return requests


def build_evidence_tables(regions, evidence_index, pad_pct=0.30, pad_floor=1000, pad_ceiling=5000):
    """Extract PE, SR, and RD evidence for every region request.

    Samples in each region are grouped by batch so a batch's PE/SR/RD
    file is queried once per region regardless of how many of its
    samples are requested, and tabix handles are cached across the
    *entire* run so a batch file shared by many regions is opened once
    total.

    Parameters
    ----------
    regions : list of RegionRequest
    evidence_index : EvidenceIndex
    pad_pct, pad_floor, pad_ceiling : see `pad_window` / `breakpoint_windows`.

    Returns
    -------
    dict of str -> pandas.DataFrame
        Keys "pe", "sr", "rd", each with a leading "name" column
        identifying which region request the rows came from.
    """
    handle_cache = TabixHandleCache()
    cov_cache = {}
    pe_frames, sr_frames, rd_frames = [], [], []

    try:
        for region in regions:
            pe_sr_windows = [
                (region.chrom, w_start, w_end)
                for w_start, w_end in breakpoint_windows(region.start, region.end, pad_pct, pad_floor, pad_ceiling)
            ]
            rd_start, rd_end = pad_window(region.start, region.end, pad_pct, pad_floor)
            rd_window = (region.chrom, rd_start, rd_end)

            batches = evidence_index.group_by_batch(region.sample_ids)
            for batch_id, batch_samples in batches.items():
                pe_url = evidence_index.pe_url.get(batch_id)
                sr_url = evidence_index.sr_url.get(batch_id)
                rd_url = evidence_index.rd_url.get(batch_id)
                cov_url = evidence_index.cov_url.get(batch_id)

                if pe_url:
                    df_pe = extract_pe(handle_cache, pe_url, pe_sr_windows, batch_samples)
                    if not df_pe.empty:
                        df_pe.insert(0, "name", region.name)
                        pe_frames.append(df_pe)

                if sr_url:
                    df_sr = extract_sr(handle_cache, sr_url, pe_sr_windows, batch_samples)
                    if not df_sr.empty:
                        df_sr.insert(0, "name", region.name)
                        sr_frames.append(df_sr)

                if rd_url:
                    df_rd = extract_rd(handle_cache, rd_url, rd_window, batch_samples, cov_url=cov_url, cov_cache=cov_cache)
                    if not df_rd.empty:
                        df_rd.insert(0, "name", region.name)
                        rd_frames.append(df_rd)
    finally:
        handle_cache.close_all()

    return {
        "pe": pd.concat(pe_frames, ignore_index=True) if pe_frames else pd.DataFrame(columns=["name"] + PE_COLUMNS),
        "sr": pd.concat(sr_frames, ignore_index=True) if sr_frames else pd.DataFrame(columns=["name"] + SR_COLUMNS),
        "rd": pd.concat(rd_frames, ignore_index=True) if rd_frames else pd.DataFrame(columns=["name"] + RD_COLUMNS),
    }
