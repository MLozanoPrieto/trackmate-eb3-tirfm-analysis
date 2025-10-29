#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
trackmate-eb3-tirfm-analysis
Version: v1.0.0
DOI: 10.5281/zenodo.XXXXXXX.
Repository: https://github.com/MLozanoPrieto/trackmate-eb3-tirfm-analysis
License: see LICENSE file in this repository

Description
-----------
Post-processing pipeline for EB3 TIRFM tracking data exported from TrackMate.
This script:
  1. Aggregates per-track and per-cell metrics from TrackMate CSVs
     (edges.csv / tracks.csv) and writes summary tables (1_METRICS).
  2. Optionally runs nonparametric statistics (Mann–Whitney / Van Elteren-style
     stratified test) and FDR correction across conditions (2_STATS).
  3. Optionally generates publication-quality plots by condition and by
     experiment, including significance annotations if stats are available
     (3_PLOTS).

Typical usage
-------------
Basic metrics only (writes 1_METRICS under the input folder):

    python trackmate_eb3_tirfm_analysis.py \
        --input examples/dummy_dataset/input

Metrics + stats (also writes 2_STATS):

    python trackmate_eb3_tirfm_analysis.py \
        --input examples/dummy_dataset/input \
        --stats

Metrics + stats + plots (also writes 3_PLOTS):

    python trackmate_eb3_tirfm_analysis.py \
        --input examples/dummy_dataset/input \
        --stats \
        --plots

Custom output folder (all results + logs go to OUTPUT_DIR instead of INPUT_DIR):

    python trackmate_eb3_tirfm_analysis.py \
        --input examples/dummy_dataset/input \
        --output OUTPUT_DIR \
        --stats \
        --plots

Citation
--------
If you use this pipeline, please cite both the software release and the Methods chapter:

Software:
    Lozano-Prieto M. TrackMate EB3 TIRFM Analysis Pipeline (v1.0.0).
    Zenodo. DOI: 10.5281/zenodo.XXXXXXX.

Methods chapter:
    Gómez-Morón Á, Lozano-Prieto M, Martín Cófreces NB.
    Studying actin and tubulin cytoskeleton dynamics at planar immunological
    synapses in primary human T lymphocytes through TIRF Microscopy.
    Methods in Molecular Biology. 2025.

Author(s)
---------
Marta Lozano-Prieto
"""


import os
import re
import sys
import shutil
import argparse
from pathlib import Path
from typing import Dict, Tuple, Optional, List
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, rankdata
import seaborn as sns
import matplotlib.pyplot as plt
from math import sqrt, erfc
from datetime import datetime


# =============================================================================
#                              LOGGING (TEE)
# =============================================================================

class _TeeStream:
    def __init__(self, original_stream, file_obj):
        self._orig = original_stream
        self._file = file_obj

    def write(self, data):
        self._orig.write(data)
        self._file.write(data)
        self._file.flush()

    def flush(self):
        self._orig.flush()
        self._file.flush()


def _setup_logging(base_output: Path) -> Path:
    """
    Create logs folder under base_output and tee stdout/stderr to a timestamped log file.
    Returns the log file path.
    """
    logs_dir = base_output / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_path = logs_dir / f"trackmate_pipeline_{ts}.log"
    f = open(log_path, "w", encoding="utf-8")
    sys.stdout = _TeeStream(sys.stdout, f)
    sys.stderr = _TeeStream(sys.stderr, f)
    print(f"[INFO] Logging to: {log_path}")
    return log_path


# =============================================================================
#                                  METRICS
# =============================================================================

# =========================
# UTILITIES
# =========================

def normalize_col_name(name: str) -> str:
    """
    Normalize column names: uppercase, replace spaces/symbols with '_',
    and remove units in parentheses.
    """
    if name is None:
        return ""
    name = re.sub(r"\(.*?\)", "", str(name))
    name = name.upper()
    name = re.sub(r"[^A-Z0-9]+", "_", name)
    name = re.sub(r"_+", "_", name).strip("_")
    return name

def has_key_column(cols, key_candidates) -> bool:
    norm_cols = [normalize_col_name(c) for c in cols]
    norm_keys = [normalize_col_name(k) for k in key_candidates]
    return any(k in norm_cols for k in norm_keys)

def detect_file_kind(fname: str) -> Optional[str]:
    """
    Return 'edges' or 'tracks' based on filename suffix.
    Ignores any '*all-spots*.csv' if present.
    """
    fl = fname.lower()
    if re.search(r"-all[-_ ]?spots\.csv$", fl):
        return None
    if fl.endswith("-edges.csv"):
        return "edges"
    if fl.endswith("-tracks.csv"):
        return "tracks"
    return None

def basename_cell_from_filename(fname: str) -> str:
    """
    Extract the common base (cell/video) part before -edges/-tracks.
    """
    fl = fname.lower()
    for suf in ["-edges.csv", "-tracks.csv"]:
        if fl.endswith(suf):
            return fname[: -len(suf)]
    return os.path.splitext(fname)[0]

def parse_metadata(parent: Path, file_path: Path) -> Tuple[str, str, str]:
    """
    From the path relative to `parent`, return (Experiment, Condition, Plate).
    Expected: parent / Experiment / Condition / [Plate] / file.csv
    """
    rel_parts = file_path.parent.relative_to(parent).parts
    exp = rel_parts[0] if len(rel_parts) >= 1 else ""
    cond = rel_parts[1] if len(rel_parts) >= 2 else ""
    plate = rel_parts[2] if len(rel_parts) >= 3 else ""
    return exp, cond, plate

def read_trackmate_csv_with_header_guess(csv_path: Path, key_candidates) -> pd.DataFrame:
    """
    Read a TrackMate CSV guessing the true header (multi-line possible).
    - Try different header offsets (first 25 lines)
    - Auto-detect delimiter
    - Try common encodings
    - Verify presence of key columns
    """
    encodings_to_try = ["utf-8", "utf-8-sig", "latin-1"]
    for enc in encodings_to_try:
        try:
            for hdr in range(0, 25):
                try:
                    df = pd.read_csv(
                        csv_path,
                        engine="python",
                        sep=None,
                        header=hdr,
                        dtype=str,
                        encoding=enc,
                        on_bad_lines="skip"
                    )
                except Exception:
                    continue
                if df.empty and len(df.columns) == 0:
                    continue
                df = df.dropna(how="all")
                if has_key_column(df.columns, key_candidates):
                    print(f"[INFO] Parsed CSV '{csv_path.name}' with encoding={enc}, header={hdr}")
                    return df
        except UnicodeDecodeError:
            continue
    df = pd.read_csv(csv_path, engine="python", sep=None, header=0, dtype=str, encoding="utf-8", on_bad_lines="skip")
    print(f"[WARN] Fallback parse CSV '{csv_path.name}' with utf-8, header=0")
    return df

def to_numeric(series: pd.Series) -> pd.Series:
    """Convert to numeric (NaN if not convertible)."""
    return pd.to_numeric(series, errors="coerce")

def pick_column(df: pd.DataFrame, candidates) -> Optional[str]:
    """
    Return the first column whose normalized name matches a candidate (EXACT only).
    No partial matches.
    """
    norm = {normalize_col_name(c): c for c in df.columns}
    for cand in candidates:
        nc = normalize_col_name(cand)
        if nc in norm:
            return norm[nc]
    return None

KEY_COLS = {
    "edges": ["TRACK_ID", "SPEED", "DIRECTIONAL_CHANGE_RATE"],
    "tracks": [
        "TRACK_MEDIAN_SPEED", "TRACK_DURATION", "TOTAL_DISTANCE_TRAVELED",
        "TRACK_DISPLACEMENT", "CONFINEMENT_RATIO",
        "LINEARITY_OF_FORWARD_PROGRESSION", "MEAN_STRAIGHT_LINE_SPEED",
        "MEAN_DIRECTIONAL_CHANGE_RATE"
    ],
}

EDGES_CANDS = {
    "TRACK_ID": ["TRACK_ID"],
    "SPEED": ["SPEED"],
    "DCR": ["DIRECTIONAL_CHANGE_RATE"],
}
TRACKS_CANDS = {
    "TRACK_MEDIAN_SPEED": ["TRACK_MEDIAN_SPEED"],
    "TRACK_DURATION": ["TRACK_DURATION"],
    "TOTAL_DISTANCE_TRAVELED": ["TOTAL_DISTANCE_TRAVELED"],
    "TRACK_DISPLACEMENT": ["TRACK_DISPLACEMENT"],
    "CONFINEMENT_RATIO": ["CONFINEMENT_RATIO"],
    "LINEARITY_OF_FORWARD_PROGRESSION": ["LINEARITY_OF_FORWARD_PROGRESSION"],
    "MEAN_STRAIGHT_LINE_SPEED": ["MEAN_STRAIGHT_LINE_SPEED"],
    "MEAN_DIRECTIONAL_CHANGE_RATE": ["MEAN_DIRECTIONAL_CHANGE_RATE"],
}

# =========================
# SCAN & AGGREGATION
# =========================

def collect_cell_files(parent_dir: Path) -> Dict[Tuple[str, str, str, str], Dict[str, Path]]:
    """
    Walk the directory tree and group file paths per cell:
    key = (Experiment, Condition, Plate, CellBase)
    value = dict with possible keys 'edges', 'tracks'
    """
    mapping: Dict[Tuple[str, str, str, str], Dict[str, Path]] = {}
    for root, _, files in os.walk(parent_dir):
        for f in files:
            if not f.lower().endswith(".csv"):
                continue
            kind = detect_file_kind(f)
            if kind is None:
                continue
            fpath = Path(root) / f
            exp, cond, plate = parse_metadata(parent_dir, fpath)
            cell_base = basename_cell_from_filename(f)
            key = (exp, cond, plate, cell_base)
            mapping.setdefault(key, {})
            mapping[key][kind] = fpath
    return mapping

def per_track_stats_from_edges(df_e: pd.DataFrame) -> Optional[pd.DataFrame]:
    """
    From an 'edges' DataFrame, return a per-track table with:
      - SPEED: median, p10, p90, IQR
      - DIRECTIONAL_CHANGE_RATE: median, p10, p90, IQR
    Cleans TRACK_ID by enforcing numeric.
    Requires TRACK_ID and at least one of SPEED or DIRECTIONAL_CHANGE_RATE.
    """
    col_tid = pick_column(df_e, EDGES_CANDS["TRACK_ID"])
    if not col_tid:
        return None

    col_speed = pick_column(df_e, EDGES_CANDS["SPEED"])
    col_dcr   = pick_column(df_e, EDGES_CANDS["DCR"])
    if not col_speed and not col_dcr:
        return None

    tid_num = pd.to_numeric(df_e[col_tid], errors="coerce")
    mask_valid = tid_num.notna()
    if not mask_valid.any():
        return None

    tmp = pd.DataFrame({col_tid: tid_num[mask_valid].astype(int)})

    if col_speed:
        tmp["_SPEED_"] = to_numeric(df_e.loc[mask_valid, col_speed])
    if col_dcr:
        tmp["_DCR_"] = to_numeric(df_e.loc[mask_valid, col_dcr])

    gb = tmp.groupby(col_tid)

    # Robust aggregation: join by TRACK_ID instead of concatenating by position
    out = pd.DataFrame(index=gb.size().index)  # index = TRACK_ID

    if col_speed:
        agg_s = gb["_SPEED_"].agg(
            SPEED_median="median",
            SPEED_p10=lambda s: s.quantile(0.10),
            SPEED_p90=lambda s: s.quantile(0.90),
            SPEED_IQR=lambda s: s.quantile(0.75) - s.quantile(0.25),
        )
        out = out.join(agg_s, how="left")

    if col_dcr:
        agg_d = gb["_DCR_"].agg(
            DIRECTIONAL_CHANGE_RATE_median="median",
            DIRECTIONAL_CHANGE_RATE_p10=lambda s: s.quantile(0.10),
            DIRECTIONAL_CHANGE_RATE_p90=lambda s: s.quantile(0.90),
            DIRECTIONAL_CHANGE_RATE_IQR=lambda s: s.quantile(0.75) - s.quantile(0.25),
        )
        out = out.join(agg_d, how="left")

    out = out.reset_index().rename(columns={"index": col_tid})
    return out if not out.empty else None

def compute_metrics_for_cell(files_dict: Dict[str, Path]) -> Tuple[Dict[str, float], Optional[pd.DataFrame]]:
    """
    Read required CSVs and compute per-cell metrics (final table) and,
    if edges are present, return the intermediate per-track table.

    Per-cell metrics from tracks.csv:
      - median TRACK_MEDIAN_SPEED
      - median TRACK_DURATION
      - p90    TRACK_DURATION
      - catastrophe_rate = 1 / mean(TRACK_DURATION)
      - medians: TOTAL_DISTANCE_TRAVELED, TRACK_DISPLACEMENT, CONFINEMENT_RATIO,
                 LINEARITY_OF_FORWARD_PROGRESSION, MEAN_STRAIGHT_LINE_SPEED,
                 MEAN_DIRECTIONAL_CHANGE_RATE
      - p90    MEAN_DIRECTIONAL_CHANGE_RATE

    Per-cell metrics derived from edges (via per-track table):
      - median of per-track (median/p10/p90/IQR) for SPEED
      - median of per-track (median/p10/p90/IQR) for DIRECTIONAL_CHANGE_RATE
    """
    out = {
        "N_TRACKS": np.nan,
        "TRACK_MEDIAN_SPEED_median": np.nan,
        "TRACK_DURATION_median": np.nan,
        "TRACK_DURATION_p90": np.nan,
        "CATASTROPHE_RATE": np.nan,
        "TOTAL_DISTANCE_TRAVELED_median": np.nan,
        "TRACK_DISPLACEMENT_median": np.nan,
        "CONFINEMENT_RATIO_median": np.nan,
        "LINEARITY_OF_FORWARD_PROGRESSION_median": np.nan,
        "MEAN_STRAIGHT_LINE_SPEED_median": np.nan,
        "MEAN_DIRECTIONAL_CHANGE_RATE_median": np.nan,
        "MEAN_DIRECTIONAL_CHANGE_RATE_p90": np.nan,
        "SPEED_TRACK_MEDIAN_median": np.nan,
        "SPEED_TRACK_P10_median": np.nan,
        "SPEED_TRACK_P90_median": np.nan,
        "SPEED_TRACK_IQR_median": np.nan,
        "DIRECTIONAL_CHANGE_RATE_TRACK_MEDIAN_median": np.nan,
        "DIRECTIONAL_CHANGE_RATE_TRACK_P10_median": np.nan,
        "DIRECTIONAL_CHANGE_RATE_TRACK_P90_median": np.nan,
        "DIRECTIONAL_CHANGE_RATE_TRACK_IQR_median": np.nan,
    }

    per_track_df = None

    if "tracks" in files_dict:
        df_t = read_trackmate_csv_with_header_guess(files_dict["tracks"], KEY_COLS["tracks"])

        col_tid = pick_column(df_t, ["TRACK_ID", "TRACK_INDEX"])
        if col_tid:
            n_tracks = pd.to_numeric(df_t[col_tid], errors="coerce").dropna().astype(int).nunique()
        else:
            # Fallback robusto: cuenta filas no vacías si no hay columna de ID detectable
            n_tracks = df_t.dropna(how="all").shape[0]
        out["N_TRACKS"] = float(n_tracks)

        def get(colkey: str) -> Optional[pd.Series]:
            col = pick_column(df_t, TRACKS_CANDS[colkey])
            return to_numeric(df_t[col]) if col else None

        s_track_median_speed = get("TRACK_MEDIAN_SPEED")
        s_track_duration = get("TRACK_DURATION")
        s_total_dist = get("TOTAL_DISTANCE_TRAVELED")
        s_track_disp = get("TRACK_DISPLACEMENT")
        s_conf_ratio = get("CONFINEMENT_RATIO")
        s_lofp = get("LINEARITY_OF_FORWARD_PROGRESSION")
        s_mean_sls = get("MEAN_STRAIGHT_LINE_SPEED")
        s_mean_dcr = get("MEAN_DIRECTIONAL_CHANGE_RATE")

        if s_track_median_speed is not None:
            out["TRACK_MEDIAN_SPEED_median"] = float(s_track_median_speed.dropna().median()) if not s_track_median_speed.dropna().empty else np.nan

        if s_track_duration is not None:
            s_clean = s_track_duration.dropna()
            if not s_clean.empty:
                out["TRACK_DURATION_median"] = float(s_clean.median())
                out["TRACK_DURATION_p90"]    = float(s_clean.quantile(0.90))
                m = float(s_clean.mean())
                out["CATASTROPHE_RATE"] = (1.0 / m) if m > 0 else np.nan

        if s_total_dist is not None:
            out["TOTAL_DISTANCE_TRAVELED_median"] = float(s_total_dist.dropna().median()) if not s_total_dist.dropna().empty else np.nan

        if s_track_disp is not None:
            out["TRACK_DISPLACEMENT_median"] = float(s_track_disp.dropna().median()) if not s_track_disp.dropna().empty else np.nan

        if s_conf_ratio is not None:
            out["CONFINEMENT_RATIO_median"] = float(s_conf_ratio.dropna().median()) if not s_conf_ratio.dropna().empty else np.nan

        if s_lofp is not None:
            out["LINEARITY_OF_FORWARD_PROGRESSION_median"] = float(s_lofp.dropna().median()) if not s_lofp.dropna().empty else np.nan

        if s_mean_sls is not None:
            out["MEAN_STRAIGHT_LINE_SPEED_median"] = float(s_mean_sls.dropna().median()) if not s_mean_sls.dropna().empty else np.nan

        if s_mean_dcr is not None:
            s_clean = s_mean_dcr.dropna()
            if not s_clean.empty:
                out["MEAN_DIRECTIONAL_CHANGE_RATE_median"] = float(s_clean.median())
                out["MEAN_DIRECTIONAL_CHANGE_RATE_p90"]    = float(s_clean.quantile(0.90))

    if "edges" in files_dict:
        df_e = read_trackmate_csv_with_header_guess(files_dict["edges"], KEY_COLS["edges"])
        per_track_df = per_track_stats_from_edges(df_e)

        if per_track_df is not None and not per_track_df.empty:
            for src_col, dst_col in [
                ("SPEED_median", "SPEED_TRACK_MEDIAN_median"),
                ("SPEED_p10", "SPEED_TRACK_P10_median"),
                ("SPEED_p90", "SPEED_TRACK_P90_median"),
                ("SPEED_IQR", "SPEED_TRACK_IQR_median"),
            ]:
                if src_col in per_track_df.columns:
                    s = to_numeric(per_track_df[src_col]).dropna()
                    out[dst_col] = float(s.median()) if not s.empty else np.nan

            for src_col, dst_col in [
                ("DIRECTIONAL_CHANGE_RATE_median", "DIRECTIONAL_CHANGE_RATE_TRACK_MEDIAN_median"),
                ("DIRECTIONAL_CHANGE_RATE_p10", "DIRECTIONAL_CHANGE_RATE_TRACK_P10_median"),
                ("DIRECTIONAL_CHANGE_RATE_p90", "DIRECTIONAL_CHANGE_RATE_TRACK_P90_median"),
                ("DIRECTIONAL_CHANGE_RATE_IQR", "DIRECTIONAL_CHANGE_RATE_TRACK_IQR_median"),
            ]:
                if src_col in per_track_df.columns:
                    s = to_numeric(per_track_df[src_col]).dropna()
                    out[dst_col] = float(s.median()) if not s.empty else np.nan

    return out, per_track_df

def build_cell_metrics_table(parent_dir: Path) -> pd.DataFrame:
    """
    Scan folders and subfolders, identify edges/tracks files (ignore any all-spots CSVs if present),
    compute the requested per-cell metrics, and return a DataFrame with columns:
      Experiment, Condition, Plate, Cell, [metrics]
    Also save:
      - '1_METRICS/cell_metrics.csv' (per-cell table)
      - '1_METRICS/track_metrics_from_edges.csv' (per-track intermediate table)
    """
    parent = parent_dir.resolve()
    cell_map = collect_cell_files(parent)

    rows: List[dict] = []
    track_tables: List[pd.DataFrame] = []

    metric_cols = [
        "N_TRACKS",
        "TRACK_MEDIAN_SPEED_median",
        "TRACK_DURATION_median",
        "TRACK_DURATION_p90",
        "CATASTROPHE_RATE",
        "TOTAL_DISTANCE_TRAVELED_median",
        "TRACK_DISPLACEMENT_median",
        "CONFINEMENT_RATIO_median",
        "LINEARITY_OF_FORWARD_PROGRESSION_median",
        "MEAN_STRAIGHT_LINE_SPEED_median",
        "MEAN_DIRECTIONAL_CHANGE_RATE_median",
        "MEAN_DIRECTIONAL_CHANGE_RATE_p90",
        "SPEED_TRACK_MEDIAN_median",
        "SPEED_TRACK_P10_median",
        "SPEED_TRACK_P90_median",
        "SPEED_TRACK_IQR_median",
        "DIRECTIONAL_CHANGE_RATE_TRACK_MEDIAN_median",
        "DIRECTIONAL_CHANGE_RATE_TRACK_P10_median",
        "DIRECTIONAL_CHANGE_RATE_TRACK_P90_median",
        "DIRECTIONAL_CHANGE_RATE_TRACK_IQR_median",
    ]

    for (exp, cond, plate, cell), files_dict in sorted(cell_map.items()):
        metrics, per_track_df = compute_metrics_for_cell(files_dict)

        metric_vals = [metrics.get(m) for m in metric_cols]
        def is_nan_or_zero(x):
            return (x is None) or (isinstance(x, float) and (np.isnan(x) or x == 0.0)) or (x == 0)
        if all(is_nan_or_zero(v) for v in metric_vals):
            print(
                f"[WARN] Cell omitted due to empty metrics (NaN/0): "
                f"{exp} / {cond} / {plate} / {cell} "
                f"| files: edges={'edges' in files_dict}, tracks={'tracks' in files_dict}"
            )
            continue

        row = {
            "Experiment": exp,
            "Condition": cond,
            "Plate": plate,
            "Cell": cell,
        }
        row.update(metrics)
        rows.append(row)

        if per_track_df is not None and not per_track_df.empty:
            per_track_df = per_track_df.copy()
            per_track_df.insert(0, "Experiment", exp)
            per_track_df.insert(1, "Condition", cond)
            per_track_df.insert(2, "Plate", plate)
            per_track_df.insert(3, "Cell", cell)
            track_tables.append(per_track_df)

    cols_order = ["Experiment", "Condition", "Plate", "Cell"] + metric_cols
    df = pd.DataFrame(rows)

    if not df.empty:
        df = df[[c for c in cols_order if c in df.columns] + [c for c in df.columns if c not in cols_order]]

        # Final content-only filter: drop rows where all metric values are NaN or 0
        metric_present = [c for c in metric_cols if c in df.columns]
        if metric_present:
            mask_all_empty = (df[metric_present].isna() | (df[metric_present] == 0)).all(axis=1)
            df = df[~mask_all_empty].copy()

    out_dir = parent / "1_METRICS"
    out_dir.mkdir(parents=True, exist_ok=True)

    out_path_cells = out_dir / "cell_metrics.csv"
    df.to_csv(out_path_cells, index=False)
    print(f"[OK] Per-cell metrics saved to: {out_path_cells}")
    print(f"[INFO] {len(df)} cells included (after filtering empties).")

    if track_tables:
        df_tracks = pd.concat(track_tables, axis=0, ignore_index=True)
    else:
        df_tracks = pd.DataFrame(columns=[
            "Experiment", "Condition", "Plate", "Cell",
            "TRACK_ID",
            "SPEED_median", "SPEED_p10", "SPEED_p90", "SPEED_IQR",
            "DIRECTIONAL_CHANGE_RATE_median", "DIRECTIONAL_CHANGE_RATE_p10",
            "DIRECTIONAL_CHANGE_RATE_p90", "DIRECTIONAL_CHANGE_RATE_IQR",
        ])

    out_path_tracks = out_dir / "track_metrics_from_edges.csv"
    df_tracks.to_csv(out_path_tracks, index=False)
    print(f"[OK] Per-track metrics (from edges) saved to: {out_path_tracks}")

    return df


def run_step1(input_dir: Path, output_base: Path):
    """Run script 1 and, if needed, move 1_METRICS from input_dir to output_base."""
    _ = build_cell_metrics_table(input_dir)

    # If output base differs from input, relocate 1_METRICS to output base
    src = input_dir.resolve() / "1_METRICS"
    dst = output_base.resolve() / "1_METRICS"
    if src != dst:
        dst.parent.mkdir(parents=True, exist_ok=True)
        if dst.exists():
            shutil.rmtree(dst)
        shutil.move(str(src), str(dst))
        print(f"[OK] Moved 1_METRICS to: {dst}")


# =============================================================================
#                               STATS
# =============================================================================

# =========================
# CONFIG (parametrized at runtime)
# =========================
CSV_PATH = Path()

# =========================
# UTILITIES
# =========================
META_COLS = ["Experiment", "Condition", "Plate", "Cell"]

def fdr_bh(pvals: pd.Series) -> pd.Series:
    """
    Benjamini–Hochberg FDR for a series of p-values (may contain NaN).
    Returns a Series aligned to the original index.
    """
    p = pvals.copy()
    mask = p.notna()
    m = mask.sum()
    if m == 0:
        return pd.Series(np.nan, index=p.index)

    p_sorted = p[mask].sort_values()
    ranks = np.arange(1, m + 1, dtype=float)
    q = (p_sorted.values * m / ranks)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    q_series = pd.Series(q, index=p_sorted.index)
    out = pd.Series(np.nan, index=p.index)
    out.loc[q_series.index] = q_series
    return out

def to_numeric_safe(df: pd.DataFrame, exclude: list) -> pd.DataFrame:
    """
    Convert all columns to numeric except those in 'exclude'.
    Non-numeric strings become NaN.
    """
    df2 = df.copy()
    for c in df2.columns:
        if c in exclude:
            continue
        df2[c] = pd.to_numeric(df2[c], errors="coerce")
    return df2

def compute_effect_sizes(g1: np.ndarray, g2: np.ndarray) -> dict:
    """
    Compute CLES and Cliff's delta from U (via rank sums).
    - CLES = U1 / (n1*n2)
    - Cliff's delta = 2*CLES - 1
    Uses average ranks to handle ties.
    """
    n1, n2 = len(g1), len(g2)
    combined = np.concatenate([g1, g2])
    ranks = rankdata(combined, method="average")
    R1 = ranks[:n1].sum()
    U1 = R1 - n1 * (n1 + 1) / 2.0
    cles = U1 / (n1 * n2) if n1 > 0 and n2 > 0 else np.nan
    cliffs = 2 * cles - 1 if np.isfinite(cles) else np.nan
    return {"U1": U1, "CLES": cles, "CliffsDelta": cliffs}

def _var_u_with_ties(combined: np.ndarray, n1: int, n2: int) -> float:
    """
    Variance of Mann–Whitney U with tie correction (normal approximation).
    Var(U) = n1*n2/12 * [ N + 1 - sum(t^3 - t)/(N*(N-1)) ]
    where t are tie group sizes in the pooled sample.
    """
    N = n1 + n2
    if N <= 1:
        return 0.0
    vc = pd.Series(combined).value_counts()
    tie_term = 0.0
    if len(vc) > 0:
        tie_term = float(((vc**3 - vc).sum()) / (N * (N - 1)))
    var_u = (n1 * n2 / 12.0) * ((N + 1.0) - tie_term)
    return float(var_u)

def stratified_mannwhitney_van_elteren(df: pd.DataFrame, metric: str, cond1: str, cond2: str) -> Tuple[float, int]:
    """
    Van Elteren-like stratified MW across Experiment strata.
    Returns (two-sided p-value, n_strata_used). Uses normal approximation with tie correction.
    Only strata where both conditions have data are used.
    """
    Z_num = 0.0
    Z_den = 0.0
    n_used = 0

    for exp in sorted(df["Experiment"].dropna().unique().tolist()):
        sub = df.loc[df["Experiment"] == exp, ["Condition", metric]].dropna()
        g1 = sub.loc[sub["Condition"] == cond1, metric].values
        g2 = sub.loc[sub["Condition"] == cond2, metric].values
        n1, n2 = len(g1), len(g2)
        if n1 == 0 or n2 == 0:
            continue

        combined = np.concatenate([g1, g2])
        ranks = rankdata(combined, method="average")
        R1 = ranks[:n1].sum()
        U1 = R1 - n1 * (n1 + 1) / 2.0
        mean_u = n1 * n2 / 2.0
        var_u = _var_u_with_ties(combined, n1, n2)
        if var_u <= 0:
            continue

        Z_num += (U1 - mean_u)
        Z_den += var_u
        n_used += 1

    if n_used == 0 or Z_den <= 0:
        return (np.nan, 0)

    Z = Z_num / sqrt(Z_den)
    # two-sided p-value via normal tail
    p = erfc(abs(Z) / sqrt(2.0))
    return (float(p), n_used)

def stats_main():
    # Load and basic checks
    if not CSV_PATH.exists():
        raise FileNotFoundError(f"CSV not found: {CSV_PATH}")

    df = pd.read_csv(CSV_PATH)

    for col in META_COLS:
        if col not in df.columns:
            raise ValueError(f"Missing required column '{col}' in the CSV.")
    df["Condition"] = df["Condition"].astype(str).str.strip()
    df = to_numeric_safe(df, exclude=META_COLS)

    # Detect conditions
    conds = df["Condition"].dropna().unique().tolist()
    if len(conds) != 2:
        print(f"[WARN] Expected exactly 2 conditions, found: {conds}. Skipping stats generation.")
        return

    cond1, cond2 = conds[0], conds[1]

    # Numeric metric columns
    metric_cols = [c for c in df.columns if c not in META_COLS and pd.api.types.is_numeric_dtype(df[c])]

    results = []

    # How many experiments overall?
    exps_all = sorted(df["Experiment"].dropna().unique().tolist())
    print(f"[INFO] Experiments detected for stats: {exps_all}")

    for col in metric_cols:
        g1_all = df.loc[df["Condition"] == cond1, col].dropna().values
        g2_all = df.loc[df["Condition"] == cond2, col].dropna().values

        n1, n2 = len(g1_all), len(g2_all)

        # Descriptives
        mean1 = float(np.mean(g1_all)) if n1 > 0 else np.nan
        mean2 = float(np.mean(g2_all)) if n2 > 0 else np.nan
        med1  = float(np.median(g1_all)) if n1 > 0 else np.nan
        med2  = float(np.median(g2_all)) if n2 > 0 else np.nan

        # Pooled Mann–Whitney
        if n1 >= 1 and n2 >= 1:
            try:
                test = mannwhitneyu(g1_all, g2_all, alternative="two-sided", method="auto")
                p_pooled = float(test.pvalue)
                es = compute_effect_sizes(g1_all, g2_all)
                U1 = float(es["U1"])
                cles = float(es["CLES"])
                cliffs = float(es["CliffsDelta"])
            except Exception:
                p_pooled = np.nan
                U1 = np.nan
                cles = np.nan
                cliffs = np.nan
        else:
            p_pooled = np.nan
            U1 = np.nan
            cles = np.nan
            cliffs = np.nan

        # Stratified Mann–Whitney (Van Elteren-like) by Experiment
        p_strat, n_strata_used = stratified_mannwhitney_van_elteren(df, col, cond1, cond2)

        # Decide which p-value to report as primary
        if n_strata_used >= 2 and np.isfinite(p_strat):
            p_primary = p_strat
            test_used = "MW_stratified_vanElteren"
        else:
            p_primary = p_pooled
            test_used = "MW_pooled"

        print(f"[STATS] {col}: test_used={test_used} | p_primary={p_primary} | "
              f"p_pooled={p_pooled} | p_strat={p_strat} | strata_used={n_strata_used}")

        # Direction (by median)
        direction = (
            f"{cond1}>{cond2}" if np.isfinite(med1) and np.isfinite(med2) and med1 > med2
            else (f"{cond2}>{cond1}" if np.isfinite(med1) and np.isfinite(med2) and med2 > med1 else "tie/NA")
        )

        results.append({
            "metric": col,
            "cond1": cond1,
            "cond2": cond2,
            "n_cond1": n1,
            "n_cond2": n2,
            "mean_cond1": mean1,
            "mean_cond2": mean2,
            "median_cond1": med1,
            "median_cond2": med2,
            "direction_by_median": direction,
            "U1": U1,
            "p_value_two_sided": p_primary,          # primary p (stratified if available, else pooled)
            "p_value_pooled": p_pooled,
            "p_value_stratified": p_strat,
            "n_strata_used": n_strata_used,
            "test_used": test_used,
            "CLES": cles,
            "Cliffs_delta": cliffs
        })

    res_df = pd.DataFrame(results)
    res_df["q_value_fdr_bh"] = fdr_bh(res_df["p_value_two_sided"])
    res_df["significant_q<0.05"] = res_df["q_value_fdr_bh"] < 0.05
    res_df = res_df.sort_values(["q_value_fdr_bh", "p_value_two_sided"], na_position="last")

    # Save
    out_dir = CSV_PATH.parents[1] / "2_STATS"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / "mannwhitney_all_metrics.csv"
    res_df.to_csv(out_file, index=False)

    print(f"Done. Saved to: {out_file}")
    print(f"Total metrics analyzed: {len(metric_cols)}")
    print("Columns in output:")
    print(", ".join(res_df.columns))


def run_step2(base_dir: Path):
    """Parameterize CSV_PATH and run stats_main()."""
    global CSV_PATH
    CSV_PATH = base_dir / "1_METRICS" / "cell_metrics.csv"
    print(f"[INFO] Running stats on: {CSV_PATH}")
    stats_main()


# =============================================================================
#                               PLOTS
# =============================================================================

# =========================
# CONFIG (parametrized at runtime)
# =========================
BASE_DIR = Path()

# One row per cell (metadata + metrics)
METRICS_CSV_PATH = Path()

# Global Mann–Whitney results per metric
STATS_CSV_PATH = Path()

# Output folders
PLOTS_DIR = Path()
OUT_DIR_COND = Path()
OUT_DIR_COND_UNCORR = Path()
OUT_DIR_COND_FDR = Path()
OUT_DIR_EXP = Path()

# Style
FIGSIZE = (4, 4.8)
DPI = 600
PALETTE_NAME = "Set2"
DOT_ALPHA = 0.65
DOT_SIZE = 3.6
JITTER = 0.18
BOX_LINEWIDTH = 1.4
MEDIAN_LINEWIDTH = 2.0
LABEL_FONTSIZE = 11.5
TICK_FONTSIZE = 10.5
STAR_FONTSIZE = 16

# Flag set at runtime depending on stats presence
HAS_STATS = False


# =========================
# UTILS
# =========================
def ensure_dirs():
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)
    OUT_DIR_COND.mkdir(parents=True, exist_ok=True)
    if HAS_STATS:
        OUT_DIR_COND_UNCORR.mkdir(parents=True, exist_ok=True)
        OUT_DIR_COND_FDR.mkdir(parents=True, exist_ok=True)
    OUT_DIR_EXP.mkdir(parents=True, exist_ok=True)

def sanitize_filename(name: str) -> str:
    safe = re.sub(r"[^\w\-_.]+", "_", str(name))
    safe = re.sub(r"_+", "_", safe).strip("_")
    return safe or "metric"

def is_numeric_series(s: pd.Series) -> bool:
    return pd.api.types.is_numeric_dtype(s)

def load_data():
    df = pd.read_csv(METRICS_CSV_PATH)

    stats = None
    try:
        if STATS_CSV_PATH.exists():
            tmp = pd.read_csv(STATS_CSV_PATH)
            required_cols = {"metric", "p_value_two_sided", "q_value_fdr_bh"}
            miss = required_cols - set(tmp.columns)
            if miss:
                print(f"[WARN] Missing columns in {STATS_CSV_PATH.name}: {sorted(miss)}. "
                      "Plots will be generated without significance annotations.")
            else:
                stats = tmp.set_index("metric")
        else:
            print(f"[WARN] {STATS_CSV_PATH} not found. "
                  "Plots will be generated without significance annotations.")
    except Exception as e:
        print(f"[WARN] Could not read statistics from {STATS_CSV_PATH} ({e}). "
              "Plots will be generated without significance annotations.")
        stats = None

    return df, stats

def get_metric_lists(df: pd.DataFrame):
    meta_cols = ["Experiment", "Condition", "Plate", "Cell"]
    for col in meta_cols:
        if col not in df.columns:
            raise ValueError(f"Required column not found in metrics CSV: '{col}'")
    metric_cols = [c for c in df.columns if c not in meta_cols and is_numeric_series(df[c])]
    if not metric_cols:
        raise ValueError("No numeric metrics were found in the CSV.")
    return meta_cols, metric_cols

def setup_style():
    sns.set_theme(context="notebook", style="white")
    plt.rcParams.update({
        "axes.labelsize": LABEL_FONTSIZE,
        "xtick.labelsize": TICK_FONTSIZE,
        "ytick.labelsize": TICK_FONTSIZE,
        "figure.dpi": DPI,
        "savefig.dpi": DPI,
        "axes.grid": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

def get_sig_label(p: float | None) -> str | None:
    """
    If p is None/NaN do not annotate; otherwise '*' if p < 0.05, else 'ns'.
    """
    if p is None or pd.isna(p):
        return None
    return "*" if float(p) < 0.05 else "ns"

def draw_sig_bracket(ax: plt.Axes, x1: float, x2: float, y: float, h: float, label: str):
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.2, c="black", clip_on=False)
    ax.text((x1 + x2) * 0.5, y + h + (h * 0.2), label,
            ha="center", va="bottom", fontsize=STAR_FONTSIZE, fontweight="bold")

def annotate_significance_boxplot(ax: plt.Axes, data: pd.DataFrame, x_col: str, y_col: str, label: str | None):
    """
    Draws a bracket only if a significance label is provided.
    """
    if label is None:
        return
    groups = [g for g in data[x_col].dropna().unique()]
    if len(groups) < 2:
        return
    groups = sorted(groups)[:2]

    y_vals = data[y_col].dropna()
    if y_vals.empty:
        return
    y_min, y_max = float(y_vals.min()), float(y_vals.max())
    span = (y_max - y_min) if (y_max > y_min) else (abs(y_max) + 1.0)

    x1, x2 = 0, 1
    y = y_max + span * 0.08
    h = span * 0.045
    ax.set_ylim(ax.get_ylim()[0], max(ax.get_ylim()[1], y_max + span * 0.2))

    draw_sig_bracket(ax, x1, x2, y, h, label)

def counts_per_group(data: pd.DataFrame, x_col: str, y_col: str, order: list[str]) -> dict[str, int]:
    counts = {}
    for cat in order:
        n = int(data.loc[data[x_col] == cat, y_col].dropna().shape[0])
        counts[cat] = n
    return counts

def set_xticklabels_with_counts(ax: plt.Axes, order: list[str], counts: dict[str, int]):
    positions = np.arange(len(order))
    labels = [f"{cat}\n(n={counts.get(cat, 0)})" for cat in order]
    ax.set_xticks(positions)
    ax.set_xticklabels(labels)


# =========================
# PLOTS
# =========================
def _single_condition_plot(data: pd.DataFrame, metric: str, cats: list[str], palette):
    fig, ax = plt.subplots(figsize=FIGSIZE)
    sns.boxplot(
        data=data, x="Condition", y=metric, ax=ax, palette=palette,
        linewidth=BOX_LINEWIDTH, showfliers=False, order=cats
    )
    sns.stripplot(
        data=data, x="Condition", y=metric, ax=ax, color="black",
        alpha=DOT_ALPHA, size=DOT_SIZE, jitter=JITTER, order=cats
    )
    ax.set_xlabel("")
    ax.set_ylabel(metric)
    ax.grid(False)

    cts = counts_per_group(data, "Condition", metric, cats)
    set_xticklabels_with_counts(ax, cats, cts)

    sns.despine(ax=ax)
    for art in ax.artists:
        art.set_edgecolor("black")
    for line in ax.lines:
        xdata, ydata = getattr(line, "get_xdata", lambda: [])(), getattr(line, "get_ydata", lambda: [])()
        if len(xdata) == 2 and ydata[0] == ydata[1]:
            line.set_linewidth(MEDIAN_LINEWIDTH)
    return fig, ax

def plot_by_condition(df: pd.DataFrame, metric: str, p_uncorr, q_fdr):
    data = df[["Condition", metric]].dropna()
    if data.empty:
        print(f"[WARN] '{metric}': no data for conditions (skipped).")
        return

    cats = sorted(data["Condition"].unique().tolist())
    palette = sns.color_palette(PALETTE_NAME, n_colors=len(cats))

    # If no statistics available, produce a single plot without creating UNCORRECTED/FDR subfolders
    if not HAS_STATS:
        fig, ax = _single_condition_plot(data, metric, cats, palette)
        fig.tight_layout()
        fpath = OUT_DIR_COND / f"{sanitize_filename(metric)}__by_condition.png"
        fig.savefig(fpath, dpi=DPI)
        plt.close(fig)
        return

    # ===== Uncorrected =====
    fig, ax = _single_condition_plot(data, metric, cats, palette)
    annotate_significance_boxplot(ax, data, "Condition", metric, get_sig_label(p_uncorr))
    fig.tight_layout()
    fpath = OUT_DIR_COND_UNCORR / f"{sanitize_filename(metric)}__by_condition_p.png"
    fig.savefig(fpath, dpi=DPI)
    plt.close(fig)

    # ===== FDR-corrected =====
    fig, ax = _single_condition_plot(data, metric, cats, palette)
    annotate_significance_boxplot(ax, data, "Condition", metric, get_sig_label(q_fdr))
    fig.tight_layout()
    fpath = OUT_DIR_COND_FDR / f"{sanitize_filename(metric)}__by_condition_q_fdr.png"
    fig.savefig(fpath, dpi=DPI)
    plt.close(fig)

def plot_by_experiment(df: pd.DataFrame, metric: str):
    exps = sorted(df["Experiment"].dropna().unique().tolist())
    if not exps:
        print(f"[WARN] No values found in 'Experiment'. Skipping per-experiment plots for {metric}.")
        return

    for exp in exps:
        sub = df.loc[df["Experiment"] == exp, ["Condition", metric]].dropna()
        if sub.empty:
            print(f"[WARN] '{metric}' in experiment '{exp}': no data (skipped).")
            continue

        cats = sorted(sub["Condition"].unique().tolist())
        palette = sns.color_palette(PALETTE_NAME, n_colors=len(cats))

        fig, ax = plt.subplots(figsize=FIGSIZE)

        # Always draw a boxplot, even with a single condition
        sns.boxplot(
            data=sub, x="Condition", y=metric, ax=ax, palette=palette,
            linewidth=BOX_LINEWIDTH, showfliers=False, order=cats
        )

        sns.stripplot(
            data=sub, x="Condition", y=metric, ax=ax, color="black",
            alpha=DOT_ALPHA, size=DOT_SIZE, jitter=JITTER, order=cats
        )

        ax.set_xlabel("")
        ax.set_ylabel(metric)
        ax.grid(False)

        cts = counts_per_group(sub, "Condition", metric, cats)
        set_xticklabels_with_counts(ax, cats, cts)

        sns.despine(ax=ax)
        for line in ax.lines:
            xdata, ydata = getattr(line, "get_xdata", lambda: [])(), getattr(line, "get_ydata", lambda: [])()
            if len(xdata) == 2 and ydata[0] == ydata[1]:
                line.set_linewidth(MEDIAN_LINEWIDTH)

        fig.tight_layout()
        fpath = OUT_DIR_EXP / f"{sanitize_filename(metric)}__experiment_{sanitize_filename(exp)}.png"
        fig.savefig(fpath, dpi=DPI)
        plt.close(fig)


# =========================
# MAIN (plots)
# =========================
def plots_main():
    setup_style()

    print("Loading data…")
    df, stats = load_data()
    _, metric_cols = get_metric_lists(df)

    global HAS_STATS
    if stats is None or stats.empty:
        HAS_STATS = False
        p_map, q_map = {}, {}
    else:
        HAS_STATS = True
        p_map = stats["p_value_two_sided"].to_dict()
        q_map = stats["q_value_fdr_bh"].to_dict()

    ensure_dirs()

    print(f"{len(metric_cols)} numeric metrics detected.")
    for i, metric in enumerate(metric_cols, start=1):
        p_uncorr = p_map.get(metric, np.nan)
        q_fdr = q_map.get(metric, np.nan)
        print(f"[{i}/{len(metric_cols)}] {metric} | p={p_uncorr} | q(FDR)={q_fdr}")

        plot_by_condition(df, metric, p_uncorr, q_fdr)
        plot_by_experiment(df, metric)

    print("\nDone.")
    if HAS_STATS:
        print(f"BY_CONDITION (uncorrected): {OUT_DIR_COND_UNCORR}")
        print(f"BY_CONDITION (FDR BH):      {OUT_DIR_COND_FDR}")
    else:
        print(f"BY_CONDITION:               {OUT_DIR_COND}")
    print(f"BY_EXPERIMENT:              {OUT_DIR_EXP}")

def run_step3(base_dir: Path):
    """Parameterize BASE_DIR and derived paths, then run plots_main()."""
    global BASE_DIR, METRICS_CSV_PATH, STATS_CSV_PATH
    global PLOTS_DIR, OUT_DIR_COND, OUT_DIR_COND_UNCORR, OUT_DIR_COND_FDR, OUT_DIR_EXP

    BASE_DIR = base_dir
    METRICS_CSV_PATH = BASE_DIR / "1_METRICS" / "cell_metrics.csv"
    STATS_CSV_PATH   = BASE_DIR / "2_STATS"   / "mannwhitney_all_metrics.csv"

    PLOTS_DIR = BASE_DIR / "3_PLOTS"
    OUT_DIR_COND = PLOTS_DIR / "BY_CONDITION"
    OUT_DIR_COND_UNCORR = OUT_DIR_COND / "UNCORRECTED"
    OUT_DIR_COND_FDR    = OUT_DIR_COND / "FDR_CORRECTED"
    OUT_DIR_EXP = PLOTS_DIR / "BY_EXPERIMENT"

    print(f"[INFO] Plotting from: {METRICS_CSV_PATH}")
    if STATS_CSV_PATH.exists():
        print(f"[INFO] Using stats from: {STATS_CSV_PATH}")
    else:
        print(f"[INFO] No stats file found at: {STATS_CSV_PATH} (plots without significance)")

    plots_main()


# =============================================================================
#                                 CLI ORCHESTRATOR
# =============================================================================

def main_cli():
    parser = argparse.ArgumentParser(
        description="Unified pipeline: (1) metrics -> (2) stats (optional) -> (3) plots (optional)."
    )
    parser.add_argument(
        "-i", "--input", required=True, type=Path,
        help="Input root directory (was PARENT_DIR)."
    )
    parser.add_argument(
        "-o", "--output", required=False, type=Path,
        help="Base output directory for 1_METRICS, 2_STATS, 3_PLOTS, logs/. If omitted, outputs go under the input directory."
    )
    parser.add_argument(
        "-s", "--stats", action="store_true",
        help="Compute Mann–Whitney statistics (creates 2_STATS)."
    )
    parser.add_argument(
        "-p", "--plots", action="store_true",
        help="Generate plots (creates 3_PLOTS). Does not require --stats."
    )

    args = parser.parse_args()

    input_dir = args.input.resolve()
    if not input_dir.exists():
        print(f"[ERROR] Input directory not found: {input_dir}", file=sys.stderr)
        sys.exit(1)

    base_output = args.output.resolve() if args.output else input_dir

    # Setup log file tee
    _setup_logging(base_output)

    print(f"[STEP 1] Metrics from: {input_dir}")
    run_step1(input_dir, base_output)

    if args.stats:
        print(f"\n[STEP 2] Statistics in: {base_output / '2_STATS'}")
        run_step2(base_output)

    if args.plots:
        print(f"\n[STEP 3] Plots in: {base_output / '3_PLOTS'}")
        run_step3(base_output)


if __name__ == "__main__":
    main_cli()
