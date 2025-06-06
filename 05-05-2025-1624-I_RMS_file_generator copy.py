#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  5 16:13:26 2025

@author: oskarjohnbruunsmith
"""

#!/usr/bin/env python3
"""
Batch‑compute RMS current from log files.

CSV columns
-----------
timestamp   ISO wall‑clock time, millisecond precision
I_RMS(A)    1‑s RMS current, 50 % overlap, rounded to 5 decimals
"""

from pathlib import Path
import re
from datetime import datetime, timedelta
import numpy as np
import pandas as pd

# -------- parameters you may tweak ------------------------------------
WINDOW_SEC     = 1.0     # RMS window length (s)
OVERLAP        = 0.50    # 50 % overlap
CURRENT_SCALE  = 0.5     # << multiply every raw current by this factor
# ----------------------------------------------------------------------

BRACKET_RE   = re.compile(r"\[([^\]]+)\]")
SPLIT_WS_RE  = re.compile(r"[ \t]+")

def load_current_file(path: Path):
    lines = path.read_text(encoding="utf‑8", errors="ignore").splitlines()
    hdr = next((i for i,l in enumerate(lines)
                if "Timestamp" in l and "Current" in l), None)
    if hdr is None:
        print(f"  ⚠  {path.name}: header not found – skipped")
        return None, None

    anchor_dt = anchor_ms = None
    times, currents = [], []

    for raw in lines[hdr + 1:]:
        m = BRACKET_RE.search(raw)
        if not m:
            continue
        try:
            wall_dt = datetime.strptime(m.group(1), "%Y-%m-%d %H:%M:%S.%f")
        except ValueError:
            continue

        data_part = raw.split("]", 1)[1].strip()
        toks = SPLIT_WS_RE.split(data_part)
        if len(toks) < 2:
            continue
        try:
            ms_ctr = float(toks[0])
            cur_A  = float(toks[-1]) * CURRENT_SCALE    # << scaling
        except ValueError:
            continue

        if anchor_dt is None:
            anchor_dt, anchor_ms = wall_dt, ms_ctr

        abs_dt = anchor_dt + timedelta(milliseconds=(ms_ctr - anchor_ms))
        times.append(abs_dt)
        currents.append(cur_A)

    if len(times) < 2:
        print(f"  ⚠  {path.name}: <2 usable rows – skipped")
        return None, None

    return times, np.asarray(currents, dtype=float)


def compute_rms(times, currents,
                window_sec=WINDOW_SEC, overlap=OVERLAP):
    t0 = times[0]
    t_sec = np.array([(t - t0).total_seconds() for t in times])
    span = t_sec[-1] - t_sec[0]
    if span == 0:
        return [], []

    fs_est = len(t_sec) / span
    win  = max(2, int(fs_est * window_sec))
    step = max(1, int(win * (1 - overlap)))

    rms_vals, rms_times = [], []
    for s in range(0, len(currents) - win, step):
        seg = currents[s:s + win]
        rms_vals.append(np.sqrt(np.mean(seg ** 2)))
        rms_times.append(times[s + win // 2])

    return rms_times, rms_vals


def main():
    base = Path("/Users/oskarjohnbruunsmith/Desktop/BP/Experimental_data/"
                "123_patch_test_pressure_test_02-05-2025")
    src  = base / "OneDrive_1_05-05-2025"
    out  = base / "I_RMS_data"
    out.mkdir(exist_ok=True)

    for txt in sorted(src.glob("*.txt")):
        times, currents = load_current_file(txt)
        if times is None:
            continue

        rms_times, rms_vals = compute_rms(times, currents)
        if not rms_vals:
            print(f"  ⚠  {txt.name}: window too long – skipped")
            continue

        pd.DataFrame({
            "timestamp": [t.isoformat(sep=" ", timespec="milliseconds")
                          for t in rms_times],
            "I_RMS(A)":  [f"{v:.5f}" for v in rms_vals]
        }).to_csv(out / f"{txt.stem}.csv", index=False)

        print(f"  ✓  {txt.stem}.csv  ({len(rms_vals)} points)")


if __name__ == "__main__":
    main()
