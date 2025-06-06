#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 18 20:18:15 2025

@author: oskarjohnbruunsmith
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Per-CSV pump-power analysis with per-sweep overlay (display-only)

  • Finds sweep length from a matching test-folder sitting in one of two roots
  • Resets time to 0 at the start of every sweep so curves overlay
  • Optional BIN_STEP-second smoothing
  • Shows one plot window per CSV and prints R² + mean power
"""


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
──────────────────────────────────────────────────────────────────────────────
Pump-Power “I_RMS vs Time-within-Sweep” Viewer
──────────────────────────────────────────────────────────────────────────────
S C O P E
---------
* This script is for **interactive inspection** of *raw* current-sensor
  captures (RMS current) recorded during pressure / flow-rate sweeps.
* You point it at a folder full of `*.csv` files that look like
    ──  108.csv
        109.csv
        …                              (any file-name is allowed, but it
                                        must **start with the numeric
                                        test-ID** of the sweep it belongs
                                        to; e.g. “108-noise.csv” also
                                        works).
* For *every single* CSV it finds, it opens a Matplotlib window that
  overlays **all sweeps** contained in that file on the same axes so you
  can see drift, warm-up, outliers, etc., at a glance.

D A T A   L A Y O U T   /   A S S U M P T I O N S
--------------------------------------------------
1.  **Current files**  
    Lives in one directory, configured via `I_RMS_DIR`.  
    Each CSV must have at least two columns

        timestamp ,  I_RMS(A)
          ISO8601    floating-point amps

2.  **Sweep folders**  
    The file name’s leading digits (“108” in *108.csv*) identify the
    *test number*.  
    Somewhere else (search path list `SWEEP_ROOTS`) there is a *folder*
    whose name encodes the sweep range and duration:

        <test_num>_<start>Hz_<end>Hz_<N>s
        e.g. 108_1Hz_1500Hz_100s

    The trailing “…_100s” tells us that **one sweep lasts 100 s**.
    Nothing inside that directory is read—only the folder name matters.

3.  **Time alignment**  
    The script aligns every raw sample to **“time within the current
    sweep (0 … N-s)”** by

      * turning the CSV’s absolute `timestamp` into elapsed-seconds,
      * dividing by *N* to label each row with a sweep-index,
      * subtracting whole multiples of *N* so each sweep starts at 0 s.

C O N F I G U R A B L E   P A R A M E T E R S
---------------------------------------------
* `BIN_STEP`        – if > 0, points are binned into that many seconds
                      (mean within each bin) to de-speckle the overlay;
                      0 disables binning entirely.
* `SUPPLY_VOLTAGE`  – multiplied by the file’s mean I_RMS to print an
                      average power consumption line.
* `SWEEP_ROOTS`     – list of directories to search for the matching
                      “108_…_N s” folders.
* `I_RMS_DIR`       – where the current CSVs live.

P L O T   A N N O T A T I O N S
-------------------------------
* Title shows:  file-name · sweep length · BIN size · R² of a *linear*
  fit of I_RMS vs t (within sweep) · average power.
* Legend shows one coloured line per sweep (or scatter cloud if
  BIN_STEP = 0).  Colours cycle through Matplotlib’s tab10 palette: Sweep 1 =
  blue, Sweep 2 = orange, etc.
* X-axis: 0 → N s (one full sweep).
* Y-axis: raw or binned I_RMS (Amps).

S T A T I S T I C S   P R I N T E D   T O   C O N S O L E
---------------------------------------------------------
* R² (coefficient of determination) of **all raw samples** vs their
  time-within-sweep → quick drift check (flat line ⇒ R² ≈ 0).
* Average I_RMS   (across *all* rows in that CSV).
* Average Power   (mean I_RMS × supply voltage).
* Number of sweeps detected, number of rows, sweep-length found.

A L G O R I T H M   S U M M A R Y
---------------------------------
1. **Load CSV** → strip & normalise column names.
2. **Locate sweep folder** for that test-ID → extract sweep length *N*.
3. **Convert timestamps** to elapsed-seconds; derive
       sweep_idx  = floor(t / N)
       t_in_sweep = t − sweep_idx · N
4. **(optional) BINNING**  
       t_bin = floor(t_in_sweep / BIN_STEP) · BIN_STEP
       overlay the sweep by plotting mean(I_RMS) for each bin.
5. **Statistics**  
       Linear regression of I_RMS vs t_in_sweep → R²  
       mean(I_RMS) → power.
6. **Plot** all sweeps on the same axes (colour-coded).
7. **Print** summary block to stdout.
8. Repeat for every CSV within `I_RMS_DIR` (depth-first).

C O M M O N   M O D S
---------------------
* Analyse voltage sweeps? Adjust `SUPPLY_VOLTAGE` and y-labels.
* Look at Flow-board power instead? Change the column selector for the
  measurement column (currently searches for “I_RMS(A)”).
* Batch-save PNGs instead of showing windows?  
  Replace the `plt.show()` call with `fig.savefig(...)`.

D E P E N D E N C I E S
-----------------------
Python ≥ 3.8, `numpy`, `pandas`, `matplotlib`, `scipy`.

──────────────────────────────────────────────────────────────────────────────
"""


from pathlib import Path
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

# ────────── configure paths / constants ──────────────────────────────────────
I_RMS_DIR     = Path("/Users/oskarjohnbruunsmith/Desktop/BP/Experimental_data/123_patch_test_pressure_test_02-05-2025/I_RMS_data")

# folders that *do* contain “…Hz_…Hz_<N>s” sub-dirs
# ── Replace the two entries below with the EXACT directories you provided:
SWEEP_ROOTS   = [
    Path(
        "/Users/oskarjohnbruunsmith/Desktop/BP/Experimental_data/"
        "123_patch_test_pressure_test_02-05-2025/flow_rate_tests"
    ),
    Path(
        "/Users/oskarjohnbruunsmith/Desktop/BP/Experimental_data/"
        "123_patch_test_pressure_test_02-05-2025/pressure_tests"
    ),
]

BIN_STEP      = 0        # 0 = no binning, 1 = 1 s, 5 = 5 s, …
SUPPLY_VOLTAGE = 2.5     # V
# ─────────────────────────────────────────────────────────────────────────────


def find_sweep_seconds(test_num: str) -> float:
    """
    Search the SWEEP_ROOTS for a directory whose *name* looks like
    "<test_num>_<start>Hz_<end>Hz_<N>s" and return N (seconds) as float.
    """
    pat = re.compile(fr"^{test_num}_[0-9]+Hz_[0-9]+Hz_(\d+)s$", re.IGNORECASE)

    for root in SWEEP_ROOTS:
        if not root.is_dir():
            continue
        for p in root.iterdir():
            m = pat.match(p.name)
            if m:
                return float(m.group(1))

    raise ValueError(f"Sweep-length folder for test {test_num} not found "
                     f"in any SWEEP_ROOTS ({', '.join(str(r) for r in SWEEP_ROOTS)})")


def _bin(values, step):
    return values if step <= 0 else np.floor(values / step) * step


def process_csv(csv_path: Path):
    # ---------- load ---------------------------------------------------------
    df = pd.read_csv(csv_path)
    df.columns = df.columns.str.strip()
    if {"timestamp", "I_RMS(A)"} - set(df.columns):
        raise ValueError(f"{csv_path}: missing 'timestamp' or 'I_RMS(A)' columns")

    # ---------- get sweep length --------------------------------------------
    test_num = re.match(r"(\d+)", csv_path.stem).group(1)   # e.g. '108' from '108.csv'
    sweep_len = find_sweep_seconds(test_num)

    # ---------- time alignment ----------------------------------------------
    df["timestamp"] = pd.to_datetime(df["timestamp"])
    df = df.sort_values("timestamp", ignore_index=True)
    t0 = df["timestamp"].iloc[0]
    df["elapsed_s"] = (df["timestamp"] - t0).dt.total_seconds()

    df["sweep_idx"]  = np.floor(df["elapsed_s"] / sweep_len).astype(int)
    df["t_in_sweep"] = df["elapsed_s"] - df["sweep_idx"] * sweep_len
    df["t_bin"]      = _bin(df["t_in_sweep"], BIN_STEP)

    # ---------- statistics on RAW points ------------------------------------
    slope, intercept, r_val, *_ = stats.linregress(df["t_in_sweep"], df["I_RMS(A)"])
    r_squared = r_val ** 2
    mean_I     = df["I_RMS(A)"].mean()
    mean_power = mean_I * SUPPLY_VOLTAGE

    # ---------- plot ---------------------------------------------------------
    cmap = plt.get_cmap("tab10")
    fig, ax = plt.subplots(figsize=(10, 5))
    
    for s_idx, g in df.groupby("sweep_idx"):
        colour = cmap(s_idx % 10)
        if BIN_STEP > 0:
            g_bin = (
                g.groupby("t_bin")["I_RMS(A)"]
                 .mean()
                 .reset_index()
                 .rename(columns={"t_bin": "t"})
            )
            ax.plot(
                g_bin["t"],
                g_bin["I_RMS(A)"],
                "-o",
                ms=4,
                color=colour,
                label=f"Sweep {s_idx + 1}"
            )
        else:
            ax.scatter(
                g["t_in_sweep"],
                g["I_RMS(A)"],
                s=6,
                alpha=0.5,
                color=colour,
                label=f"Sweep {s_idx + 1}"
            )
    
    # Determine the maximum I_RMS value across all raw samples
    y_max = df["I_RMS(A)"].max() * 1.3
    y_min = df["I_RMS(A)"].min() * 0.7
    
    ax.set(
        xlabel="Time within sweep (s)",
        ylabel="I$_{RMS}$ (A)",
        title=(
            f"{csv_path.name}  |  sweep = {sweep_len:.0f}s  |  "
            f"BIN = {BIN_STEP}s  |  R$^2$ = {r_squared:.4f}  |  "
            f"Mean P = {mean_power*1e3:.2f} mW"
        ),
        xlim=(0, sweep_len),
    )
    ax.set_ylim(y_min, y_max)
    
    ax.grid(True)
    ax.legend(frameon=False, fontsize="small", ncol=2)
    plt.tight_layout()
    plt.show()


    # ---------- console summary ---------------------------------------------
    print("──────────────────────────────────────────────────────────────")
    print(f"File              : {csv_path.name}")
    print(f"Sweep length      : {sweep_len:.1f} s")
    print(f"Number of sweeps  : {df['sweep_idx'].nunique()}")
    print(f"Samples           : {len(df):,}")
    print(f"R² (I_RMS vs tₛ)  : {r_squared:.6f}")
    print(f"Average I_RMS     : {mean_I:.6f} A")
    print(f"Average power @ {SUPPLY_VOLTAGE:.2f} V : {mean_power:.6f} W")
    print("──────────────────────────────────────────────────────────────\n")


def main():
    csv_files = sorted(p for p in I_RMS_DIR.rglob("*.csv") if p.is_file())
    if not csv_files:
        raise FileNotFoundError(f"No CSV files found under {I_RMS_DIR}")

    for csv_path in csv_files:
        process_csv(csv_path)


if __name__ == "__main__":
    main()
