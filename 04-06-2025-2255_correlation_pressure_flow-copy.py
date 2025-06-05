#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

"""
Two plots per patch + two pooled plots (“ALL”)
Correct SD handling, with correlation between flow & pressure in curve plots
"""


"""
04-06-2025-2255_correlation_pressure_flow.py

Overview:
This script processes raw flow‐rate and pressure CSV recordings for each “patch” (experimental sample),
bins them by frequency, computes per‐patch and pooled statistics, and generates dual‐axis curves and
scatter plots illustrating how flow and pressure vary with frequency, as well as their mutual correlation.
It is used in Section 4.6 of the thesis to demonstrate the relationship between flow rate and pressure for
each patch and for the combined dataset.

Purpose:
  • Recursively scan a root directory (ROOT) for subfolders named like “<testID>_<start>Hz_<end>Hz_<duration>s”.
  • Within each such folder, detect whether it belongs to a “pressure_tests” branch or a “flow_rate_tests” branch:
      – If in “pressure_tests”, load the CSV’s pressure column (mbar) as the measurement.
      – Otherwise, load the CSV’s flow column (µL/min).
  • For each CSV file:
      1. Read into a pandas DataFrame, drop any “Air detected” rows ±2 seconds around each flagged timestamp.
      2. Convert the timestamp column to pandas datetime, compute elapsed seconds relative to the start.
      3. Compute a “Freq” column: each elapsed second t maps to a frequency in [start_freq, end_freq] via
         Freq = start + (end − start) × ((t mod sweep_duration) / sweep_duration).
  • Group the cleaned data by patch ID (extracted from folder names or accompanying .txt files). Maintain
    two dictionaries: flows[patchID] = list of per‐run DataFrames (flow vs. frequency); presses[patchID] = list
    of per‐run DataFrames (pressure vs. frequency).
  • For each patchID (iterating in sorted numeric order):
      1. Build “pooled curve” DataFrames:
         – curve_flow   = pooled_mean_sd(flowL, CURVE_BIN_HZ): For a list of flow DataFrames, bin by rounding
                         Frequency to the nearest multiple of CURVE_BIN_HZ (e.g., 3 Hz), then compute
                         {mean, std} across all runs in that patch for each bin.
         – curve_press  = pooled_mean_sd(pressL, CURVE_BIN_HZ) similarly for pressure DataFrames.
      2. If either curve_flow or curve_press is non‐empty:
         a. Identify common frequency bins where both flow and pressure have data.
         b. Compute Pearson r_cp and Spearman ρ_cp between the mean flow and mean pressure arrays at those bins.
         c. Plot a dual‐axis figure:
            • Left y-axis: flow (µL/min) vs. frequency (Hz) with error bars = std (if curve_flow exists).
            • Right y-axis: pressure (mbar) vs. frequency (Hz) with error bars = std (if curve_press exists).
            • Title: “Patch <pid> | Flow Rate & Pressure vs Frequency (Hz) | (<CURVE_BIN_HZ> Hz Binned)
                       \nPearson r=<r_cp> | Spearman ρ=<ρ_cp>” (only if both flow and pressure exist).
            • Legends: “Flow (µL/min)” on left, “Pressure (mbar)” on right.
            • X‐axis limits: [FREQ_MIN_HZ, FREQ_MAX_HZ].
            • Save to “<ROOT>/analysis/patch-<pid>/patch-<pid>_combined_freq.png”.
      3. If both flowL and pressL are non‐empty:
         a. Build “pooled scatter” DataFrames:
            – scat_flow  = pooled_mean_sd(flowL, SCATTER_BIN_HZ): mean±std at each frequency bin (SCATTER_BIN_HZ).
            – scat_press = pooled_mean_sd(pressL, SCATTER_BIN_HZ).
            – merged = scat_flow.join(scat_press, lsuffix="_flow", rsuffix="_press").
         b. If merged is non‐empty:
            – Compute Pearson r and Spearman ρ between merged["mean_flow"] and merged["mean_press"].
            – Plot a scatter with xerr = std_flow, yerr = std_press:
              • ax.errorbar(mean_flow, mean_press, xerr=std_flow, yerr=std_press,
                            fmt='o', capsize=ERROR_CAPSIZE, elinewidth=ERROR_ELINEWIDTH,
                            alpha=ERROR_ALPHA, ecolor='grey', markersize=POINT_SIZE)
              • X‐label: “Mean flow (µL/min)”, Y‐label: “Mean pressure (mbar)”
              • Title: “Patch <pid> | Flow vs Pressure | <SCATTER_BIN_HZ> Hz Binned
                         \nPearson r=<r> | Spearman ρ=<ρ>”
              • Grid on.
              • Save to “<ROOT>/analysis/patch-<pid>/patch-<pid>_flow_vs_press.png”.
      4. Accumulate all per‐patch flow Frames into all_flow_frames and all per‐patch pressure Frames into all_press_frames.

  • After processing every patchID:
      1. Build “grand-pooled” curves:
         – curve_all_flow  = pooled_mean_sd(all_flow_frames, CURVE_BIN_HZ)
         – curve_all_press = pooled_mean_sd(all_press_frames, CURVE_BIN_HZ)
         – common = intersection of their bin indices.
         – If common is non‐empty, compute Pearson r_ga and Spearman ρ_ga between grand‐mean flow & grand‐mean pressure.
         – Plot dual‐axis “ALL combined freq”:
            • Left axis: grand pooled flow vs. frequency with error bars.
            • Right axis: grand pooled pressure vs. frequency with error bars.
            • Title: “Grand average – Flow & Pressure vs Frequency (Hz) (BIN <CURVE_BIN_HZ> Hz)
                       \nPearson r=<r_ga>, Spearman ρ=<ρ_ga>”
            • Save to “<ROOT>/analysis/ALL/ALL_combined_freq.png”.
      2. Build “grand-pooled” scatter:
         – scat_all_flow  = pooled_mean_sd(all_flow_frames, SCATTER_BIN_HZ)
         – scat_all_press = pooled_mean_sd(all_press_frames, SCATTER_BIN_HZ)
         – mergedG = scat_all_flow.join(scat_all_press, lsuffix="_flow", rsuffix="_press")
         – If mergedG is non‐empty:
            • Compute Pearson rG, Spearman ρG between mergedG["mean_flow"] and mergedG["mean_press"].
            • Plot error‐bar scatter: xerr=std_flow, yerr=std_press
              – X‐label: “Mean flow (µL/min)”, Y‐label: “Mean pressure (mbar)”
              – Title: “Grand average – Flow vs Pressure (BIN <SCATTER_BIN_HZ> Hz)
                         \nPearson r=<rG>, Spearman ρ=<ρG>”
              – Grid on.
              – Save to “<ROOT>/analysis/ALL/ALL_flow_vs_press.png”.

  • At the end, if any files or folders were “skipped” (because parsing or loading failed), print them to console.

How It Works (Step-by-Step):

1. Configuration & Constants
   – `ROOT`: Root directory containing both “flow_rate_tests” and “pressure_tests” subdirectories.
   – `CURVE_BIN_HZ`: Bin width (Hz) for curve‐generation (default 3 Hz).
   – `SCATTER_BIN_HZ`: Bin width (Hz) for scatter‐plots (default 3 Hz).
   – `FREQ_MIN_HZ`, `FREQ_MAX_HZ`: Frequency limits for plots (e.g., 0 to 1500 Hz).
   – `POINT_SIZE`: Marker size for scatter/line plots.
   – `ERROR_CAPSIZE`, `ERROR_ELINEWIDTH`, `ERROR_ALPHA`: Error‐bar styling parameters.
   – `TIME_OFFSET`: Seconds to subtract from each elapsed time (if logger started late).
   – `SAVE_DPI`: Dots‐per‐inch for saved PNGs (e.g., 300).
   – `VERBOSE`: If True, print every successfully loaded CSV to console.
   – Regex patterns:
     • `PATCH_TXT_RE = re.compile(r"[Pp]atch\s*0*?(\d+)")`: Extract numeric patch ID from filenames like “patch1” or “Patch 01”.
     • `SWEEP_RE = re.compile(r"^\d+_[0-9]+Hz_[0-9]+Hz_[0-9]+s$")`: Identify subfolders containing sweep data.
     • `SWEEP_INFO_RE = re.compile(r"(\d+)Hz_(\d+)Hz_(\d+)s")`: Extract start_freq, end_freq, duration (in seconds) from folder name.

2. Helper Functions
   – `sweep_info(name)`: Given a folder name matching SWEEP_INFO_RE, returns `(s_f, e_f, d)` as floats.
   – `load_csv(path, is_pressure)`:
       • Read CSV at `path` with pandas.read_csv(...).
       • If only one column detected, re‐read without delimiter argument.
       • Detect “Air detected” column (optional).
       • If `is_pressure=True`, find a column matching `IPS\s*\(`, else find a column containing “Flowboard” or “Flow Unit”.
       • Rename the DataFrame’s columns to `[“Time”, “Value”]` (and optionally `[“Time”, “Value”, “Air”]` if air column exists).
       • Convert “Time” to datetime. If “Air” exists and any row has Air=True, remove ±2 seconds around those timestamps.
       • Return a DataFrame with columns `[“Time”, “Value”]`. Errors raised if no Value column or if all rows removed.

   – `add_freq(df, sweep)`:
       • Given a DataFrame `df` with columns `[“Time”, “Value”]` and a folder name `sweep` (e.g., “10Hz_500Hz_60s”):
           1. Extract `(s_f, e_f, d)` via `sweep_info(sweep)`.
           2. Compute `Elapsed_s = (Time − first_time).total_seconds()`.
           3. Compute `Freq = s_f + (e_f − s_f) × ((Elapsed_s − TIME_OFFSET) mod d) / d`.
           4. Return a new DataFrame with columns `["Freq", "Value"]`.
       • Effectively, each timestamp is mapped to a corresponding instantaneous frequency in the sweep.

   – `bin_means_only(df, step)`:
       • Round `df["Freq"]` to nearest multiple of `step`.
       • Filter bins in [FREQ_MIN_HZ, FREQ_MAX_HZ].
       • Return a Series indexed by bin value, containing the mean of `Value` within each bin.

   – `pooled_mean_sd(frames, step)`:
       • Input: `frames` = list of DataFrames, each with columns `["Freq","Value"]` (one per run).
       • For each DataFrame f in `frames`, compute `bin_means_only(f, step)` (Series).
       • Concatenate those Series along columns to form a DataFrame “table” with rows = frequency bins, columns per run.
       • Return a new DataFrame with index = frequency bins, columns:
           { "mean": row‐wise mean across runs, "std": row‐wise standard deviation (ddof=1) across runs }.

   – `save(fig, path)`:
       • Ensure parent directories of `path` exist.  
       • Save `fig` at `path` with DPI=SAVE_DPI.  
       • Call `plt.show()` (to display in interactive sessions).

3. Scanning Data
   – Loop `for sweep in ROOT.rglob("*")`: recursively walk through all directories under ROOT.
   – If `sweep` is a directory and its name matches SWEEP_RE (e.g., “10_100Hz_1Hz_60s”), proceed:
       • Attempt to detect `patchID` by searching for a “*.txt” file in `sweep` whose stem matches PATCH_TXT_RE.
       • If no patch‐ID found (no matching txt), append `(relative_path, "no Patch txt")` to `skipped`, continue.
       • Determine `is_press = True if "pressure_tests" in sweep.parts else False`.
       • Set `dest = presses if is_press else flows`.
       • For each CSV file in `sweep.glob("*.csv")`:
           1. Try `tidy = load_csv(csv, is_press)`. If it fails, record in `skipped` with error message.
           2. If successful, call `df_freq = add_freq(tidy, sweep.name)`.
           3. Append `df_freq` to `dest[patchID]` and, if VERBOSE, print a success message to console.

4. Per-Patch Plotting
   – Initialize two empty lists: `all_flow_frames = []`, `all_press_frames = []`.
   – For each unique `patchID` (sorted numerically by int):
       • `flowL = flows.get(patchID, [])` (list of DataFrames of flow vs. frequency).
       • `pressL = presses.get(patchID, [])` (list of DataFrames of pressure vs. frequency).
       • `out = ROOT / "analysis" / f"patch-{patchID}"` → output directory for this patch.
       • Build per-patch “curve” DataFrames:
           – `curve_flow  = pooled_mean_sd(flowL, CURVE_BIN_HZ)`.
           – `curve_press = pooled_mean_sd(pressL, CURVE_BIN_HZ)`.
       • If `curve_flow` or `curve_press` is non‐empty:
           a. Compute `common = intersection(curve_flow.index, curve_press.index)`.
           b. If `common` non‐empty:
              – `cf = curve_flow.loc[common, "mean"]` (Series of mean flow at common bins).
              – `cp = curve_press.loc[common, "mean"]` (Series of mean pressure).
              – `r_cp, p_cp = pearsonr(cf, cp)`; `rho_cp, _ = spearmanr(cf, cp)`.
           c. Create figure `fig, axF = plt.subplots(figsize=(8,5))`, assign `axF` to flow axis.
           d. Title string:  
              ```
              f"Patch {patchID} | Flow Rate & Pressure vs Frequency (Hz) |  ({CURVE_BIN_HZ} Hz Binned)"
              + (f"\nPearson r={r_cp:.3f} |  Spearman ρ={rho_cp:.3f}" if common non‐empty else "")
              ```
           e. If `curve_flow` non‐empty:
              – `axF.errorbar(curve_flow.index, curve_flow["mean"], curve_flow["std"], fmt='-o', capsize=ERROR_CAPSIZE, elinewidth=ERROR_ELINEWIDTH, alpha=ERROR_ALPHA, markersize=POINT_SIZE, label="Flow (µL/min)")`.
              – `axF.set_ylabel("Flow (µL/min)")`.
           f. If `curve_press` non‐empty:
              – If `curve_flow` non‐empty: `axP = axF.twinx()`; else `axP = axF` (single‐axis).  
              – `axP.errorbar(curve_press.index, curve_press["mean"], curve_press["std"], fmt='-s', capsize=ERROR_CAPSIZE, elinewidth=ERROR_ELINEWIDTH, alpha=ERROR_ALPHA, markersize=POINT_SIZE, color="C1", label="Pressure (mbar)")`.  
              – `axP.set_ylabel("Pressure (mbar)")`.  
           g. `axF.set(xlabel="Frequency (Hz)", xlim=(FREQ_MIN_HZ, FREQ_MAX_HZ))`, `axF.grid(True)`, `axF.legend(loc="upper left", fontsize="small")`.  
              If both flow & pressure present: `axP.legend(loc="upper right", fontsize="small")`.  
           h. `plt.tight_layout()`; `save(fig, out / f"patch-{patchID}_combined_freq.png")`, `plt.close(fig)`.

       • If `flowL` and `pressL` are both non‐empty:
           a. `scat_flow  = pooled_mean_sd(flowL, SCATTER_BIN_HZ)`.
           b. `scat_press = pooled_mean_sd(pressL, SCATTER_BIN_HZ)`.
           c. `merged = scat_flow.join(scat_press, lsuffix="_flow", rsuffix="_press")`.
           d. If `merged` non‐empty:
               – `r, pval = pearsonr(merged["mean_flow"], merged["mean_press"])`.
               – `rho, _ = spearmanr(merged["mean_flow"], merged["mean_press"])`.
               – Create `fig, ax = plt.subplots(figsize=(6,5))`.
               – `ax.errorbar(merged["mean_flow"], merged["mean_press"], xerr=merged["std_flow"], yerr=merged["std_press"], fmt="o", capsize=ERROR_CAPSIZE, elinewidth=ERROR_ELINEWIDTH, alpha=ERROR_ALPHA, ecolor="grey", markersize=POINT_SIZE)`.
               – `ax.set(xlabel="Mean flow (µL/min)", ylabel="Mean pressure (mbar)")`.
               – `ax.set_title(f"Patch {patchID} |  Flow vs Pressure | {SCATTER_BIN_HZ} Hz Binned\nPearson r={r:.3f} |  Spearman ρ={rho:.3f}")`.
               – `ax.grid(True)`, `plt.tight_layout()`.
               – `save(fig, out / f"patch-{patchID}_flow_vs_press.png")`, `plt.close(fig)`.

       • Append all individual flow DataFrames in `flowL` to `all_flow_frames` and all `pressL` to `all_press_frames`.

5. Grand-Pooled Plots (across all patches)
   – `all_dir = ROOT / "analysis" / "ALL"`.
   – `curve_all_flow  = pooled_mean_sd(all_flow_frames, CURVE_BIN_HZ)`.
   – `curve_all_press = pooled_mean_sd(all_press_frames, CURVE_BIN_HZ)`.
   – If either `curve_all_flow` or `curve_all_press` is non‐empty:
       a. Compute `common = intersection(curve_all_flow.index, curve_all_press.index)`.
       b. If `common` non‐empty:
          – `caf = curve_all_flow.loc[common, "mean"]`; `cap = curve_all_press.loc[common, "mean"]`.
          – `r_ga, p_ga = pearsonr(caf, cap)`; `rho_ga, _ = spearmanr(caf, cap)`.
       c. Create `fig, axF = plt.subplots(figsize=(8,5))`. Title:  
          ```
          f"Grand average – Flow & Pressure vs Frequency (Hz) (BIN {CURVE_BIN_HZ} Hz)"
          + (f"\nPearson r={r_ga:.3f}, Spearman ρ={rho_ga:.3f}" if common non‐empty else "")
          ```
       d. If `curve_all_flow` non‐empty:
          – `axF.errorbar(curve_all_flow.index, curve_all_flow["mean"], curve_all_flow["std"], fmt='-o', capsize=ERROR_CAPSIZE, elinewidth=ERROR_ELINEWIDTH, alpha=ERROR_ALPHA, markersize=POINT_SIZE, label="Flow (µL/min)")`; `axF.set_ylabel("Flow (µL/min)")`.
       e. If `curve_all_press` non‐empty:
          – If flow plotted: `axP = axF.twinx()`; else `axP = axF`.
          – `axP.errorbar(curve_all_press.index, curve_all_press["mean"], curve_all_press["std"], fmt='-s', capsize=ERROR_CAPSIZE, elinewidth=ERROR_ELINEWIDTH, alpha=ERROR_ALPHA, markersize=POINT_SIZE, color="C1", label="Pressure (mbar)")`; `axP.set_ylabel("Pressure (mbar)")`.
       f. `axF.set(xlabel="Frequency (Hz)", xlim=(FREQ_MIN_HZ, FREQ_MAX_HZ))`, `axF.grid(True)`, `axF.legend(loc="upper left", fontsize="small")`. If both plotted: `axP.legend(loc="upper right", fontsize="small")`.
       g. `plt.tight_layout()`, `save(fig, all_dir / "ALL_combined_freq.png")`, `plt.close(fig)`.

   – Create Grand-Pooled Scatter:
       a. `scat_all_flow  = pooled_mean_sd(all_flow_frames, SCATTER_BIN_HZ)`.
       b. `scat_all_press = pooled_mean_sd(all_press_frames, SCATTER_BIN_HZ)`.
       c. `mergedG = scat_all_flow.join(scat_all_press, lsuffix="_flow", rsuffix="_press")`.
       d. If `mergedG` non‐empty:
          – `rG, pG = pearsonr(mergedG["mean_flow"], mergedG["mean_press"])`.
          – `rhoG, _ = spearmanr(mergedG["mean_flow"], mergedG["mean_press"])`.
          – Create `fig, ax = plt.subplots(figsize=(6,5))`.
          – `ax.errorbar(mergedG["mean_flow"], mergedG["mean_press"], xerr=mergedG["std_flow"], yerr=mergedG["std_press"], fmt="o", capsize=ERROR_CAPSIZE, elinewidth=ERROR_ELINEWIDTH, alpha=ERROR_ALPHA, ecolor="grey", markersize=POINT_SIZE)`.
          – `ax.set(xlabel="Mean flow (µL/min)", ylabel="Mean pressure (mbar)")`.
          – `ax.set_title(f"Grand average – Flow vs Pressure (BIN {SCATTER_BIN_HZ} Hz)\nPearson r={rG:.3f}, Spearman ρ={rhoG:.3f}")`.
          – `ax.grid(True)`, `plt.tight_layout()`.
          – `save(fig, all_dir / "ALL_flow_vs_press.png")`, `plt.close(fig)`.

6. Diagnostics & Skipped Files
   – If `skipped` list is non‐empty:
       • Print a header “Skipped items:” then each `(path, warning)`.

Usage:
 1. Set `ROOT` to the parent folder containing two subdirectories: “flow_rate_tests” and “pressure_tests”.  
 2. Adjust bin widths `CURVE_BIN_HZ` and `SCATTER_BIN_HZ` to desired frequency resolution.  
 3. Run `python 04-06-2025-2255_correlation_pressure_flow.py`.  
 4. Inspect console:  
    – “✓ path/to/csv → patch <pid>” for every loaded CSV.  
    – If any CSV fails, it will be listed under “Skipped items.”  
 5. Output: PNG files in `<ROOT>/analysis/patch-<pid>/` for each patch and `<ROOT>/analysis/ALL/` for pooled plots.

Outputs Produced:
  • For each patch `<pid>`:  
      – `patch-<pid>_combined_freq.png`: Dual‐axis curve: Flow vs. Frequency and Pressure vs. Frequency, with binned means & std.  
      – `patch-<pid>_flow_vs_press.png`: Scatter plot of mean flow vs. mean pressure (binned), with error bars & correlation in title.  
  • For all patches combined:  
      – `ALL_combined_freq.png`: Dual‐axis curve for grand‐pooled data.  
      – `ALL_flow_vs_press.png`: Scatter for grand‐pooled mean flow vs. mean pressure.  
  • No CSV/TXT side‐files are produced by this script. All output is PNGs.


"""






from pathlib import Path
import re, numpy as np, pandas as pd, matplotlib.pyplot as plt
from collections import defaultdict
from scipy import stats

# ───────── USER SETTINGS ────────────────────────────────────────────────────
ROOT            = Path(
    "/Users/oskarjohnbruunsmith/Desktop/BP/Experimental_data/"
    "123_patch_test_pressure_test_02-05-2025"
)
CURVE_BIN_HZ    = 3
SCATTER_BIN_HZ  = 3
FREQ_MIN_HZ     = 0
FREQ_MAX_HZ     = 1500
POINT_SIZE      = 3           # marker area (pts²)
# Use thinner error‐bar style:
ERROR_CAPSIZE   = 2           # capsize (pts)
ERROR_ELINEWIDTH = 0.4        # elinewidth (pts)
ERROR_ALPHA     = 0.7         # transparency
TIME_OFFSET     = 0.0
SAVE_DPI        = 300
VERBOSE         = True
# ────────────────────────────────────────────────────────────────────────────

PATCH_TXT_RE  = re.compile(r"[Pp]atch\s*0*?(\d+)")
SWEEP_RE      = re.compile(r"^\d+_[0-9]+Hz_[0-9]+Hz_[0-9]+s$")
SWEEP_INFO_RE = re.compile(r"(\d+)Hz_(\d+)Hz_(\d+)s")

flows, presses = defaultdict(list), defaultdict(list)
skipped        = []

# ═════════ helpers ══════════════════════════════════════════════════════════
def sweep_info(name):
    s,e,d = map(float, SWEEP_INFO_RE.search(name).groups())
    return s,e,d

def load_csv(path,is_pressure):
    raw = pd.read_csv(path, sep=";", engine="python")
    if raw.shape[1]==1:
        raw = pd.read_csv(path)
    air = next((c for c in raw.columns if "Air detected" in c), None)
    if is_pressure:
        col = next((c for c in raw.columns if re.search(r"IPS\s*\(", c)), None)
    else:
        col = next((c for c in raw.columns if "Flowboard" in c or "Flow Unit" in c), None)
    if col is None:
        raise ValueError("no Flow/Pressure column")
    df = raw[[raw.columns[0], col] + ([air] if air else [])].copy()
    df.columns = ["Time", "Value"] + (["Air"] if air else [])
    df["Time"] = pd.to_datetime(df["Time"])
    if "Air" in df and df["Air"].any():
        t0 = df["Time"].iloc[0]
        ts = (df["Time"] - t0).dt.total_seconds().to_numpy()
        bad = np.zeros(len(df), bool)
        for t in ts[df["Air"].to_numpy(bool)]:
            bad |= np.abs(ts - t) <= 2
        df = df[~bad]
    if df.empty:
        raise ValueError("all rows removed")
    return df[["Time", "Value"]]

def add_freq(df, sweep):
    s,e,d = sweep_info(sweep)
    span = e - s
    t = (df["Time"] - df["Time"].iloc[0]).dt.total_seconds()
    df = df.copy()
    df["Freq"] = s + span * (((t - TIME_OFFSET) % d) / d)
    return df[["Freq", "Value"]]

def bin_means_only(df, step):
    bins = np.round(df["Freq"]/step) * step
    df = df.assign(Bin=bins)
    df = df[df["Bin"].between(FREQ_MIN_HZ, FREQ_MAX_HZ)]
    return df.groupby("Bin")["Value"].mean()

def pooled_mean_sd(frames, step):
    if not frames:
        return pd.DataFrame()
    table = pd.concat([bin_means_only(f, step) for f in frames], axis=1)
    return pd.DataFrame({
        "mean": table.mean(axis=1),
        "std":  table.std(axis=1, ddof=1)
    })

def save(fig, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=SAVE_DPI)
    plt.show()

# ═════════ scan data ════════════════════════════════════════════════════════
for sweep in ROOT.rglob("*"):
    if not (sweep.is_dir() and SWEEP_RE.match(sweep.name)):
        continue
    pid = None
    for txt in sweep.glob("*.txt"):
        m = PATCH_TXT_RE.search(txt.stem)
        if m:
            pid = m.group(1)
            break
    if pid is None:
        skipped.append((sweep.relative_to(ROOT), "no Patch txt"))
        continue
    is_press = "pressure_tests" in sweep.parts
    dest = presses if is_press else flows
    for csv in sweep.glob("*.csv"):
        try:
            tidy = load_csv(csv, is_press)
            dest[pid].append(add_freq(tidy, sweep.name))
            if VERBOSE:
                print(f"✓ {csv.relative_to(ROOT)} → patch {pid}")
        except Exception as e:
            skipped.append((csv.relative_to(ROOT), str(e)))

# ═════════ per-patch plots ══════════════════════════════════════════════════
all_flow_frames = []
all_press_frames = []
for pid in sorted(set(flows) | set(presses), key=int):
    flowL, pressL = flows.get(pid, []), presses.get(pid, [])
    out = ROOT / "analysis" / f"patch-{pid}"

    # pooled curve
    curve_flow  = pooled_mean_sd(flowL, CURVE_BIN_HZ)
    curve_press = pooled_mean_sd(pressL, CURVE_BIN_HZ)

    if not curve_flow.empty or not curve_press.empty:
        # correlation between flow and pressure means
        common = curve_flow.index.intersection(curve_press.index)
        if not common.empty:
            cf = curve_flow.loc[common, "mean"]
            cp = curve_press.loc[common, "mean"]
            r_cp, p_cp = stats.pearsonr(cf, cp)
            rho_cp, _ = stats.spearmanr(cf, cp)
        fig, axF = plt.subplots(figsize=(8,5))
        title = f"Patch {pid} | Flow Rate & Pressure vs Frequency (Hz) |  ({CURVE_BIN_HZ} Hz Binned)"
        if not common.empty:
            title += f"\nPearson r={r_cp:.3f} |  Spearman ρ={rho_cp:.3f}"
        fig.suptitle(title)

        if not curve_flow.empty:
            axF.errorbar(
                curve_flow.index, curve_flow["mean"], curve_flow["std"],
                fmt="-o",
                capsize=ERROR_CAPSIZE,
                elinewidth=ERROR_ELINEWIDTH,
                alpha=ERROR_ALPHA,
                markersize=POINT_SIZE,
                label="Flow (µL/min)"
            )
            axF.set_ylabel("Flow (µL/min)")
        if not curve_press.empty:
            axP = axF.twinx() if not curve_flow.empty else axF
            axP.errorbar(
                curve_press.index, curve_press["mean"], curve_press["std"],
                fmt="-s",
                capsize=ERROR_CAPSIZE,
                elinewidth=ERROR_ELINEWIDTH,
                alpha=ERROR_ALPHA,
                markersize=POINT_SIZE,
                color="C1",
                label="Pressure (mbar)"
            )
            axP.set_ylabel("Pressure (mbar)")

        axF.set(xlabel="Frequency (Hz)", xlim=(FREQ_MIN_HZ, FREQ_MAX_HZ))
        axF.grid(True)
        axF.legend(loc="upper left", fontsize="small")
        if not curve_flow.empty and not curve_press.empty:
            axP.legend(loc="upper right", fontsize="small")
        plt.tight_layout()
        save(fig, out / f"patch-{pid}_combined_freq.png")

    # scatter remains
    if flowL and pressL:
        scat_flow  = pooled_mean_sd(flowL, SCATTER_BIN_HZ)
        scat_press = pooled_mean_sd(pressL, SCATTER_BIN_HZ)
        merged = scat_flow.join(scat_press, lsuffix="_flow", rsuffix="_press")
        if not merged.empty:
            r, pval = stats.pearsonr(merged["mean_flow"], merged["mean_press"])
            rho, _ = stats.spearmanr(merged["mean_flow"], merged["mean_press"])
            fig, ax = plt.subplots(figsize=(6,5))
            ax.errorbar(
                merged["mean_flow"], merged["mean_press"],
                xerr=merged["std_flow"], yerr=merged["std_press"],
                fmt="o",
                capsize=ERROR_CAPSIZE,
                elinewidth=ERROR_ELINEWIDTH,
                alpha=ERROR_ALPHA,
                ecolor="grey",
                markersize=POINT_SIZE
            )
            ax.set(xlabel="Mean flow (µL/min)", ylabel="Mean pressure (mbar)")
            ax.set_title(
                f"Patch {pid} |  Flow vs Pressure | {SCATTER_BIN_HZ} Hz Binned\n"
                f"Pearson r={r:.3f} |  Spearman ρ={rho:.3f}"
            )
            ax.grid(True)
            plt.tight_layout()
            save(fig, out / f"patch-{pid}_flow_vs_press.png")

    all_flow_frames.extend(flowL)
    all_press_frames.extend(pressL)

# ═════════ ALL plots – pooled curve stats ════════════════════════════════════
all_dir = ROOT / "analysis" / "ALL"
curve_all_flow  = pooled_mean_sd(all_flow_frames, CURVE_BIN_HZ)
curve_all_press = pooled_mean_sd(all_press_frames, CURVE_BIN_HZ)

if not curve_all_flow.empty or not curve_all_press.empty:
    common = curve_all_flow.index.intersection(curve_all_press.index)
    if not common.empty:
        caf = curve_all_flow.loc[common, "mean"]
        cap = curve_all_press.loc[common, "mean"]
        r_ga, p_ga = stats.pearsonr(caf, cap)
        rho_ga, _ = stats.spearmanr(caf, cap)
    fig, axF = plt.subplots(figsize=(8,5))
    title = f"Grand average – Flow & Pressure vs Frequency (Hz) (BIN {CURVE_BIN_HZ} Hz)"
    if not common.empty:
        title += f"\nPearson r={r_ga:.3f}, Spearman ρ={rho_ga:.3f}"
    fig.suptitle(title)
    if not curve_all_flow.empty:
        axF.errorbar(
            curve_all_flow.index, curve_all_flow["mean"], curve_all_flow["std"],
            fmt="-o",
            capsize=ERROR_CAPSIZE,
            elinewidth=ERROR_ELINEWIDTH,
            alpha=ERROR_ALPHA,
            markersize=POINT_SIZE,
            label="Flow (µL/min)"
        )
        axF.set_ylabel("Flow (µL/min)")
    if not curve_all_press.empty:
        axP = axF.twinx() if not curve_all_flow.empty else axF
        axP.errorbar(
            curve_all_press.index, curve_all_press["mean"], curve_all_press["std"],
            fmt="-s",
            capsize=ERROR_CAPSIZE,
            elinewidth=ERROR_ELINEWIDTH,
            alpha=ERROR_ALPHA,
            markersize=POINT_SIZE,
            color="C1",
            label="Pressure (mbar)"
        )
        axP.set_ylabel("Pressure (mbar)")
    axF.set(xlabel="Frequency (Hz)", xlim=(FREQ_MIN_HZ, FREQ_MAX_HZ))
    axF.grid(True)
    axF.legend(loc="upper left", fontsize="small")
    if not curve_all_flow.empty and not curve_all_press.empty:
        axP.legend(loc="upper right", fontsize="small")
    plt.tight_layout()
    save(fig, all_dir / "ALL_combined_freq.png")

# grand-pooled scatter unchanged
scat_all_flow  = pooled_mean_sd(all_flow_frames, SCATTER_BIN_HZ)
scat_all_press = pooled_mean_sd(all_press_frames, SCATTER_BIN_HZ)
mergedG = scat_all_flow.join(scat_all_press, lsuffix="_flow", rsuffix="_press")
if not mergedG.empty:
    rG, pG = stats.pearsonr(mergedG["mean_flow"], mergedG["mean_press"])
    rhoG, _ = stats.spearmanr(mergedG["mean_flow"], mergedG["mean_press"])
    fig, ax = plt.subplots(figsize=(6,5))
    ax.errorbar(
        mergedG["mean_flow"], mergedG["mean_press"],
        xerr=mergedG["std_flow"], yerr=mergedG["std_press"],
        fmt="o",
        capsize=ERROR_CAPSIZE,
        elinewidth=ERROR_ELINEWIDTH,
        alpha=ERROR_ALPHA,
        ecolor="grey",
        markersize=POINT_SIZE
    )
    ax.set(xlabel="Mean flow (µL/min)", ylabel="Mean pressure (mbar)")
    ax.set_title(
        f"Grand average – Flow vs Pressure (BIN {SCATTER_BIN_HZ} Hz)\n"
        f"Pearson r={rG:.3f}, Spearman ρ={rhoG:.3f}"
    )
    ax.grid(True)
    plt.tight_layout()
    save(fig, all_dir / "ALL_flow_vs_press.png")

# diagnostic
if skipped:
    print("\nSkipped items:")
    for f,w in skipped:
        print(f"  {f} -> {w}")
