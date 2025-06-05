#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


FULL SWEEP‑ANALYSIS SUITE
=========================
• Flow‑rate (µL/min)  *or*  Pressure (mbar) — auto‑detected per CSV.  
• Patch‑comparison   *or*  8‑configuration comparison — auto‑detected.  
• Adjustable frequency bin width via FREQ_BIN_HZ.  
• Generates overview, all‑tests, average‑per‑label, parameter‑effects
  and global‑effects plots, each with & without error bars.  
"""

# ─────────────────────── USER SETTINGS ──────────────────────────────────────
PARENT_FOLDER = r"/Users/oskarjohnbruunsmith/Desktop/BP/Experimental_data/123_patch_test_pressure_test_02-05-2025/pressure_tests"
TIME_OFFSET   = 0.0        # seconds offset, if the logger starts late
FREQ_BIN_HZ   = 5         # 1 → every Hz, 2 → every 2 Hz, 10 → every 10 Hz …
POINT_SIZE    = 6          # global marker-size for all plots

# Individual families (set 0 to skip)
PLOT_8_CONFIGS            = 1    # ignored when comparing patches
PLOT_ALL_TESTS_PER_CONFIG = 1
PLOT_AVG_PER_LABEL        = 1
PLOT_PARAMETER_EFFECTS    = 1    # configs‑only
PLOT_GLOBAL_EFFECTS       = 1    # configs‑only
# ----------------------------------------------------------------------------

import os, re, glob, math, warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
try:
    from power_consumption_code import consumption_piecewise
except ModuleNotFoundError:
    consumption_piecewise = None

warnings.filterwarnings("ignore", category=FutureWarning)

TITLE_SUFFIX = "| 2.5V"

# ═════════════════════ utilities ════════════════════════════════════════════
def ensure_dir(p): os.makedirs(p, exist_ok=True)
def round_bin(series, step): return np.round(series / step) * step

# ═════════════════════ sweep‑parameter parsing ══════════════════════════════
SWEEP_RE = re.compile(r"(\d+Hz_\d+Hz_\d+s)")
def parse_sweep_params(name):
    m = SWEEP_RE.search(name)
    if not m:
        raise ValueError(f"{name}: needs ‘…Hz_…Hz_…s’ substring")
    s, e, t = m.group(1).split('_')
    return float(s[:-2]), float(e[:-2]), float(t[:-1])

# ═════════════════════ CSV loader (Flow OR Pressure) ════════════════════════
MEASUREMENT_NAME = UNIT = None  # set on first CSV

def load_and_process_csv(folder, s_f, e_f, sweep_s):
    global MEASUREMENT_NAME, UNIT

    # pick CSV
    csv = next((c for c in glob.glob(os.path.join(folder, "trimmed_*.csv"))), None)
    if csv:
        df = pd.read_csv(csv, sep=';', engine='python')
        if len(df) < 2:
            csv = None
    if csv is None:
        cand = [c for c in glob.glob(os.path.join(folder, "*.csv"))
                if not os.path.basename(c).startswith("trimmed_")]
        if not cand:
            raise FileNotFoundError(f"No CSV in {folder}")
        csv = cand[0]
        df = pd.read_csv(csv, sep=';', engine='python')

    # drop setpoints & find measurement column
    df.drop(columns=[c for c in df.columns if 'Setpoint' in c],
            inplace=True, errors='ignore')
    time_col = df.columns[0]

    flow_col  = next((c for c in df.columns
                      if 'Flowboard' in c and 'Air' not in c), None)
    press_col = next((c for c in df.columns if re.search(r'IPS\s*\(\d+\)', c)), None) \
             or next((c for c in df.columns if 'Pressureboard' in c), None)

    if flow_col:
        meas_col, name, unit = flow_col, "Flow Rate", "µL/min"
    elif press_col:
        meas_col, name, unit = press_col, "Pressure", "mbar"
    else:
        raise ValueError(f"{csv}: no Flowboard / IPS / Pressureboard column")

    if MEASUREMENT_NAME is None:
        MEASUREMENT_NAME, UNIT = name, unit
    elif name != MEASUREMENT_NAME:
        raise ValueError("Mixed measurement types in one batch.")

    # assemble dataframe
    df = df[[time_col, meas_col]].copy()
    df.columns = ['Time', 'Value']
    df['Time'] = pd.to_datetime(df['Time'])
    df['Elapsed_s'] = (df['Time'] - df['Time'].iloc[0]).dt.total_seconds()

    span = e_f - s_f
    df['Frequency'] = df['Elapsed_s'].apply(
        lambda t: s_f + span * (((t - TIME_OFFSET) % sweep_s) / sweep_s)
    )
    return df

# ═════════════════════ label detection (patch / config) ═════════════════════
PATCH_RE = re.compile(r'patch[_\s\-]?([123])', re.I)

def detect_patch(folder):
    hit = PATCH_RE.search(folder)
    if hit:
        return f"Patch {hit.group(1)}"  # space added
    txt = next((t for t in glob.glob(os.path.join(folder, "*.txt"))), "")
    hit = PATCH_RE.search(os.path.basename(txt))
    return f"Patch {hit.group(1)}" if hit else None

def config_from_txt(folder):
    txt = next((t for t in glob.glob(os.path.join(folder, "*.txt"))), "")
    nm = os.path.basename(txt).lower()
    pp = "push" if "push" in nm else "pull" if "pull" in nm else "unknown"
    dd = "opposite" if "opposite" in nm else "original" if "original" in nm else "unknown"
    ff = "air" if "air" in nm else "water" if "water" in nm else "unknown"
    return f"{pp}_system_{dd}_direction_{ff}"

# ═════════════════════ per‑folder processing ════════════════════════════════
def process_folder(path):
    s_f, e_f, sweep_s = parse_sweep_params(os.path.basename(path))
    df = load_and_process_csv(path, s_f, e_f, sweep_s)

    # binning
    df['FreqBin'] = round_bin(df['Frequency'], FREQ_BIN_HZ)
    g = df.groupby('FreqBin')['Value']
    df_avg = pd.DataFrame({
        'Frequency': g.mean().index,
        'Value':      g.mean().values,
        'StdDev':     g.std(ddof=1).fillna(0).values
    })

    # efficiency (Flow mode only)
    if MEASUREMENT_NAME == "Flow Rate" and consumption_piecewise:
        dens, A = 1000, math.pi * (0.001 / 2) ** 2
        Q = df['Value'] * 1e-9 / 60
        df['Efficiency'] = 0.5 * dens * (Q ** 3) / A ** 2 \
                           / consumption_piecewise(df['Frequency'])
        g2 = df.groupby('FreqBin')['Efficiency']
        df_eff = pd.DataFrame({
            'Frequency': g2.mean().index,
            'Efficiency': g2.mean().values,
            'EffStd': g2.std(ddof=1).fillna(0)
        })
    else:
        df_eff = pd.DataFrame()

    label = detect_patch(path) or config_from_txt(path)
    return label, df_avg, df_eff

# ═════════════════════ recursive collector ═════════════════════════════════
RUN_RE = re.compile(r'\d+Hz_\d+Hz_\d+s')

def looks_like_run(name):
    return bool(RUN_RE.search(name))

def collect_all(root):
    info = {}
    for dirpath, dirnames, _ in os.walk(root):
        if os.path.basename(dirpath) == "plots":
            dirnames[:] = []
            continue
        if not looks_like_run(os.path.basename(dirpath)):
            continue
        try:
            label, avg, eff = process_folder(dirpath)
        except Exception as exc:
            print("[skip]", os.path.relpath(dirpath, root), "→", exc)
            continue
        info.setdefault(label, []).append((os.path.relpath(dirpath, root), avg, eff))
    return info

# ═════════════════════ plotting helpers ════════════════════════════════════
def draw(ax, x, y, yerr, label, errbar):
    if errbar:
        ax.errorbar(x, y, yerr=yerr, fmt='o-', markersize=POINT_SIZE, capsize=4, label=label)
    else:
        ax.plot(x, y, marker='o', linestyle='-', markersize=POINT_SIZE, label=label)

def make_plots(df_map, base_dir, errbar, compare_patches):
    """
    Create all figure families (overview, per-test, etc.).
    Titles now follow the unified structure used in your first script.
    """
    # helper for consistent titles ------------------------------------------
    def fmt_title(context):
        return f"{MEASUREMENT_NAME} vs Frequency  |  {context}  |  2.5 V"

    tag  = 'with_errorbars' if errbar else 'without_errorbars'
    root = os.path.join(base_dir, tag)
    ensure_dir(root)
    ylabel = f"{MEASUREMENT_NAME} [{UNIT}]"

    # ── 1) overview ────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 6))
    for lbl in sorted(df_map):
        frames = [d for _, d, __ in df_map[lbl] if not d.empty]
        if not frames:
            continue
        merged = pd.concat(frames, ignore_index=True)
        g = merged.groupby('Frequency')['Value']
        draw(ax, g.mean().index, g.mean().values, g.std(ddof=1).fillna(0).values, lbl, errbar)


        context = "All patches – averaged values comparison"

    ax.set(title=fmt_title(context), xlabel="Frequency [Hz]", ylabel=ylabel)
    ax.grid(alpha=.3)
    ax.legend(fontsize='small')
    fig.tight_layout()
    ensure_dir(os.path.join(root, 'overview'))
    fig.savefig(os.path.join(root, 'overview', 'overview.png'), dpi=300)
    plt.close(fig)

    # ── 2) all tests per label ────────────────────────────────────────────
    if PLOT_ALL_TESTS_PER_CONFIG:
        for lbl, lst in df_map.items():
            fig, ax = plt.subplots(figsize=(10, 6))
            for sub, avg, _ in lst:
                if avg.empty:
                    continue
                tn = os.path.basename(sub).split('_')[0]
                draw(ax, avg['Frequency'], avg['Value'], avg['StdDev'],
                     f"Test {tn}", errbar)
            ctx = f"{lbl} – individual tests"
            ax.set(title=fmt_title(ctx), xlabel="Frequency [Hz]", ylabel=ylabel)
            ax.grid(alpha=.3)
            ax.legend(fontsize='small')
            fig.tight_layout()
            out = os.path.join(root, 'all_tests', lbl)
            ensure_dir(out)
            fig.savefig(os.path.join(out, 'all_tests.png'), dpi=300)
            plt.close(fig)

    # ── 3) averaged per label ─────────────────────────────────────────────
    if PLOT_AVG_PER_LABEL:
        for lbl, lst in df_map.items():
            frames = [d for _, d, __ in lst if not d.empty]
            if not frames:
                continue
            merged = pd.concat(frames, ignore_index=True).groupby('Frequency')
            avg = merged['Value'].mean()
            std = merged['Value'].std(ddof=1).fillna(0)
            fig, ax = plt.subplots(figsize=(10, 6))
            draw(ax, avg.index, avg.values, std.values, lbl, errbar)
            ctx = f"{lbl} – averaged values"
            ax.set(title=fmt_title(ctx), xlabel="Frequency [Hz]", ylabel=ylabel)
            ax.grid(alpha=.3)
            ax.legend([lbl])
            fig.tight_layout()
            out = os.path.join(root, 'avg_per_label', lbl)
            ensure_dir(out)
            fig.savefig(os.path.join(out, 'avg.png'), dpi=300)
            plt.close(fig)

    if compare_patches:
        return  # stop here for patch-comparison mode

    # ── 4) parameter-effects (configs only) ───────────────────────────────
    if PLOT_PARAMETER_EFFECTS:
        combined = {}
        for cfg, lst in df_map.items():
            pp, dd, ff = cfg.split('_system_')[0], \
                         cfg.split('_direction_')[0].split('_')[-1], \
                         cfg.split('_')[-1]
            frames = [avg for _, avg, __ in lst if not avg.empty]
            if not frames:
                continue
            merged = pd.concat(frames, ignore_index=True).groupby('Frequency')
            combined[(pp, dd, ff)] = pd.DataFrame({
                'Frequency': merged['Value'].mean().index,
                'Value':     merged['Value'].mean().values,
                'StdDev':    merged['Value'].std(ddof=1).fillna(0).values
            })

        def plot_pair(k1, k2, folder, name):
            fig, ax = plt.subplots(figsize=(10, 6))
            for k in (k1, k2):
                if k not in combined:
                    continue
                frame = combined[k]
                draw(ax, frame['Frequency'], frame['Value'],
                     frame['StdDev'], k[0], errbar)
            ax.set(title=fmt_title(name.replace('_', ' ')),
                   xlabel="Frequency [Hz]", ylabel=ylabel)
            ax.grid(alpha=.3)
            ax.legend(fontsize='small')
            fig.tight_layout()
            ensure_dir(folder)
            fig.savefig(os.path.join(folder, f"{name}.png"), dpi=300)
            plt.close(fig)

        outP = os.path.join(root, 'parameter_effects')
        for dd in ['original', 'opposite']:
            for ff in ['water', 'air']:
                plot_pair(('push', dd, ff), ('pull', dd, ff),
                          os.path.join(outP, 'push_vs_pull'),
                          f"push_vs_pull_{dd}_{ff}")
        for pp in ['push', 'pull']:
            for ff in ['water', 'air']:
                plot_pair((pp, 'original', ff), (pp, 'opposite', ff),
                          os.path.join(outP, 'original_vs_opposite'),
                          f"original_vs_opposite_{pp}_{ff}")
        for pp in ['push', 'pull']:
            for dd in ['original', 'opposite']:
                plot_pair((pp, dd, 'air'), (pp, dd, 'water'),
                          os.path.join(outP, 'air_vs_water'),
                          f"air_vs_water_{pp}_{dd}")

    # ── 5) global effects (configs only) ──────────────────────────────────
    if PLOT_GLOBAL_EFFECTS:
        def merge_frames(keyword):
            acc = []
            for cfg, lst in df_map.items():
                if keyword in cfg:
                    acc += [avg for _, avg, __ in lst if not avg.empty]
            if not acc:
                return None
            merged = pd.concat(acc, ignore_index=True).groupby('Frequency')
            return pd.DataFrame({
                'Frequency': merged['Value'].mean().index,
                'Value':     merged['Value'].mean().values,
                'StdDev':    merged['Value'].std(ddof=1).fillna(0).values
            })

        outG = os.path.join(root, 'global_effects')
        for a, b in [('push', 'pull'),
                     ('air', 'water'),
                     ('original', 'opposite')]:
            fa, fb = merge_frames(a), merge_frames(b)
            if fa is None or fb is None:
                continue
            fig, ax = plt.subplots(figsize=(10, 6))
            draw(ax, fa['Frequency'], fa['Value'], fa['StdDev'], a, errbar)
            draw(ax, fb['Frequency'], fb['Value'], fb['StdDev'], b, errbar)
            ctx = f"global {a} vs {b}"
            ax.set(title=fmt_title(ctx), xlabel="Frequency [Hz]", ylabel=ylabel)
            ax.grid(alpha=.3)
            ax.legend(fontsize='small')
            fig.tight_layout()
            ensure_dir(outG)
            fig.savefig(os.path.join(outG, f"global_{a}_vs_{b}.png"), dpi=300)
            plt.close(fig)


# ═════════════════════ MAIN ═════════════════════════════════════════════════
def main():
    if not os.path.isdir(PARENT_FOLDER):
        print("❌  Bad PARENT_FOLDER:", PARENT_FOLDER)
        return

    plots_root = os.path.join(PARENT_FOLDER, "plots")
    ensure_dir(plots_root)

    info = collect_all(PARENT_FOLDER)
    if not info:
        print("❌  No run folders found beneath",PARENT_FOLDER, "\n    (expected directories whose names contain e.g. ‘1Hz_100Hz_1s’)")
        return

    compare_patches = any(re.match(r'Patch[123]', k, re.I) for k in info)

    base = os.path.join(plots_root,
        'pressure_plots' if MEASUREMENT_NAME.startswith('Pressure') else 'flow_rate_plots')
    ensure_dir(base)

    for errbar in (True, False):
        make_plots(info, base, errbar, compare_patches)

    if MEASUREMENT_NAME == "Flow Rate" and consumption_piecewise:
        print("Efficiency plots kept from earlier version.")
    elif MEASUREMENT_NAME == "Flow Rate":
        print("Efficiency plots skipped – power_consumption_code not found.")

    mode = "patch comparison" if compare_patches else "configuration comparison"
    print(f"✓ All {mode} {MEASUREMENT_NAME.lower()} plots generated.")
    print("  Frequency bin width =", FREQ_BIN_HZ, "Hz")
    print("  Output root        =", plots_root)

if __name__ == "__main__":
    main()
