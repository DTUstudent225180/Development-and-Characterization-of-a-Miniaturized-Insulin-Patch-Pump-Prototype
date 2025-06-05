#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 27 00:01:30 2025

@author: oskarjohnbruunsmith
"""

"""
Flow / Pressure sweep processing – quiet version
(All plots are written under: <parent_folder>/plots_2/)
"""


"""
flow_and_pressure_analysis.py

Overview:
This script processes raw flow (or pressure) data recorded during frequency sweeps,
computes summary statistics, and produces a comprehensive suite of analytical plots
used in Chapter 4 of the thesis.

Purpose:
  • Ingest multiple subfolders whose names encode a frequency sweep (e.g., “10Hz_500Hz_60s”),
    each containing a CSV file with timestamped flow (µL/min) or pressure (mbar) measurements.
  • Clean and preprocess each CSV:
      – Filter out rows flagged as “Air detected” (±2 seconds around any air detection event).
      – Convert timestamps into elapsed seconds since the start of each run.
      – Map elapsed time into a “Frequency” column using (start_freq, end_freq, sweep_duration),
        so each measurement is assigned a frequency based on its position within each repeated sweep.
      – Tag each row with a “Sweep” index (integer number of complete sweeps elapsed).
  • Compute “binned” summary statistics at user-defined bin widths:
      – BIN_STEP_INDIVIDUAL (for per-test, high-granularity plots)
      – BIN_STEP_COMBINED   (for combined multi-test plots)
      – BIN_STEP_AVERAGE    (for average/aggregate plots)
  • Generate per-test visualizations:
      1. Flow vs. Frequency (or Pressure vs. Frequency) with mean±std at each bin, one trace per sweep.
      2. Flow vs. Elapsed Time (or Pressure vs. Elapsed Time) with mean±std, one trace per sweep.
      3. Highlight the “best” bin (maximum mean flow or pressure) with a star marker.
      4. Residual diagnostic plots (Normal Q–Q, residuals vs. fitted, histogram, residuals vs. run order).
  • Generate combined multi-test visualizations:
      1. Overlay all tests’ binned points for Flow vs. Frequency (or Time), with and without error bars.
      2. Compute a global average curve across all tests and plot mean±std vs. Frequency (or Time).
      3. Highlight the global maximum flow or pressure bin with a star marker.
      4. Plot “all raw points” scatter (every data point) vs. Frequency (or Time) for each test.
      5. Plot “relative” curves: each test’s flow (or pressure) values normalized to its own maximum (0–100%).
  • Identify high-flow/high-stability “cloud” regions:
      – For each test’s binned (mean, std), find the top 25% by mean-flow (or pressure),
        then within that subset take the bottom 10% by std. Highlight those bins and the single
        most stable bin with a star marker.
      – Repeat the same on the combined data across all tests to produce an “Overall” cloud.
  • Perform residual diagnostics on both per-test and combined data:
      – Generate a 2×2 grid of plots (Normal Q–Q, residuals vs. fitted, histogram, residuals vs. run)
        to check distributional assumptions and identify outliers.
  • Compute hydraulic efficiency (flow mode only):
      1. Take the aggregated binned-average flow (µL/min) vs. Frequency.
      2. Convert to volumetric flow (m³/s), compute velocity (v) through a fixed-area tube.
      3. Compute hydraulic power P_hyd = 0.5·ρ·Q·v² and propagate uncertainty from flow std to P_hyd_std.
      4. Divide by a constant electrical input P_elec (W) to obtain efficiency η = P_hyd / P_elec (unitless),
         along with its uncertainty η_std.
      5. Plot η vs. Frequency both with error bars (η ± η_std) and as point-only markers.
      6. Create dual-axis plots overlaying mean flow vs. Frequency on one axis and η vs. Frequency on the
         other axis. Highlight:
         – The frequency at which η is maximized (star marker, with η±η_std in the legend).
         – The efficiency at 483 Hz (star marker) if that frequency exists in the data.
      7. Align the zero points on the two y-axes so that the baseline of flow and efficiency axes coincide.
  • Analyze variability (standard deviation) vs. mean flow:
      – For each test and for the combined data, bin data by frequency, compute (mean, std) pairs.
      – Create scatter plots of std vs. mean, color-coded by frequency bin.
      – Perform linear regression (std = m·mean + b) and compute Pearson r and Spearman ρ.
      – Annotate each plot’s title with r and ρ values.
  • Conduct Exploratory Data Analysis (EDA):
      – Per-test raw-data histograms (one histogram of flow or pressure values per test).
      – Per-test boxplots of raw flow or pressure values.
      – Across-all-tests binned-average histogram (distribution of mean values across bins).
      – Across-all-tests boxplot (distribution of all raw flow or pressure points in one boxplot),
        with CSV and TXT summaries of quartiles, whiskers, and outlier counts.
      – Correlation heatmaps:
         1. Raw correlation matrix of (Frequency, Elapsed Time, Flow or Pressure).
         2. Aggregated correlation matrix of (Frequency, mean, std) for the binned data.
      – Combined side-by-side boxplots showing the distribution of flow or pressure for each individual test,
        with numeric summary CSV and TXT output.
  • All plots are saved under a new “plots_2” directory in the same parent folder, organized into subfolders:
      – freq_individual_with_err, freq_individual_no_err, time_individual_with_err, time_individual_no_err
      – freq_combined_series, freq_combined_series_err, time_combined_series, time_combined_series_err
      – freq_average_with_err, freq_average_no_err, time_average_with_err, time_average_no_err
      – all_points_with_err/frequency, all_points_with_err/time, all_points_without_err/frequency, all_points_without_err/time
      – relative/frequency
      – high_flow_high_stability
      – residuals
      – efficiency_plots
      – std_vs_flow
      – exploratory_plots/individual_histograms, individual_boxplots, average_histogram, average_boxplot, raw_correlation, aggregate_correlation, combined_boxplots

How It Works (Step-by-Step):
1. Configuration & Constants
   – Set parent_folder to the root directory containing subfolders named like “<start>Hz_<end>Hz_<duration>s”.
   – Choose BIN_STEP_INDIVIDUAL, BIN_STEP_COMBINED, BIN_STEP_AVERAGE to control bin widths (Hz or seconds).
   – Define PRESSURE_TEST = False for flow mode (µL/min) or True for pressure mode (mbar).
   – Set fluid density ρ, tube diameter, supply voltage, and electrical input power P_elec.

2. Helper Functions
   – load_and_process_csv(folder, s_f, e_f, sweep_s):
       • Load the CSV (prefer “trimmed_*.csv” if exists; else any “*.csv”).
       • Rename columns to [Time, Value], detect “Air detected” column, drop rows within ±2 s of any Air flag.
       • Convert “Time” to pandas datetime, compute “Elapsed_s” from start.
       • Compute “Sweep” index = floor(Elapsed_s / sweep_s).
       – Compute “Frequency” = s_f + (e_f – s_f)·((Elapsed_s mod sweep_s) / sweep_s).
       • Return DataFrame with columns [Time, Flow (or Pressure), Elapsed_s, Frequency, Sweep].

   – _bin(value, step):
       • If step <= 0, return value unchanged; else return round(value/step)*step.

   – add_best_marker(ax, x_best, y_best, std_best=None):
       • Plot a gold star at (x_best, y_best), with optional vertical error bar ±std_best.
       • Label the star as “★ <y> ±<std> @ <x> Hz”.

   – residual_plots(resid, out, title):
       • Generate a 2×2 grid: Normal Q–Q, residuals vs. fitted (fitted=0), histogram of resid, residuals vs. run order.
       • Save as PNG to the specified out path.

   – print_data_summary(tests):
       • Concatenate all test DataFrames.
       • Print variable names and dtypes, total rows, total measurement time in minutes/seconds, and any missing values.

3. Per-Test Plotting
   – indiv_plot_freq(tests, out_dir, err=True|False):
       • For each (folder, df) in tests:
           – Extract test_id, s_f, e_f, sweep_duration via extract_test_info.
           – For each group of rows sharing the same “Sweep” index:
               • Compute bin = _bin(Frequency, BIN_STEP_INDIVIDUAL).
               • Compute aggregated DataFrame ag = group.groupby('bin')['Flow'].agg(['mean','std']).
               • Plot ag['bin'] vs. ag['mean'] ± ag['std'] (if err=True) or just ag['mean'] (if err=False).
               • Keep a proxy handle for the legend (“Sweep i”).
           – Compute overall best bin: group entire df by bin, find bin with maximum mean. Plot a gold star via add_best_marker.
           – Set labels: X=“Frequency [Hz]”, Y=YLABEL, Title includes test_id, s_f→e_f Hz, sweep_duration s, supply voltage.
           – Save to out_dir/test_<test_id>.png (dpi=300).

   – indiv_plot_time(tests, out_dir, err=True|False):
       • Same as indiv_plot_freq, but bin on “Elapsed_s” instead of “Frequency”.
       • X‐axis = “Elapsed Time [s]”, title updated accordingly.
       • Save to out_dir/test_<test_id>.png.

4. Combined Series Plotting
   – comb_series(tests, out_dir, vs_freq=True|False, error_bars=True|False):
       • If vs_freq=True: key = “Frequency”; else key = “Elapsed_s”.
       • step = BIN_STEP_COMBINED.
       • Initialize a single figure and axes.
       • For each (folder, df) in tests:
           – Compute df['bin'] = _bin(df[key], step).
           – Compute grp_mean = df.groupby('bin')['Flow'].mean().
           – Compute grp_std  = df.groupby('bin')['Flow'].std().
           – If error_bars=True:  
               • Plot grp_mean vs. grp_std as errorbar.  
               • Create a large proxy error‐bar for the legend (“Test <test_id>”).  
             Else:  
               • Plot grp_mean as scatter.  
               • Create a large colored dot for the legend.
       • Concatenate all tests’ (key,Flow) pairs, find global‐max flow bin, plot a gold star via add_best_marker.
       • Set X label (Frequency or Time), Y label = YLABEL, Title indicating combined series, supply voltage, and whether error bars are used.
       • Legend with three columns of proxy handles.
       • Save to out_dir/combined_series_{err|noerr}.png.

5. Average Across All Tests
   – avg_plot(tests, out_dir, vs_freq=True|False, err=True|False):
       • key = “Frequency” if vs_freq else “Elapsed_s”.
       • Concatenate all tests’ DataFrames, select only [key, Flow].
       • Compute all_df['bin'] = _bin(all_df[key], BIN_STEP_AVERAGE).
       • Compute ag = all_df.groupby('bin')['Flow'].agg(['mean','std']).reset_index().
       • If err=True:  
           • Plot ag['bin'], ag['mean'] ± ag['std'] as error bars.  
         Else:  
           • Plot ag['bin'], ag['mean'] as line+marker.
       • Identify best bin where mean is maximal, plot a star via add_best_marker.
       • Set X, Y labels, Title indicating average of all tests, with/without error bars, supply voltage.
       • Legend shows only the star’s annotation.
       • Save as out_dir/average.png.

6. All Raw Points & Relative Plots
   – all_points(tests, err=True|False):
       • For each (folder, df):
           – Extract test_id, s_f, e_f, sweep_duration.
           – Frequency scatter: for each Sweep group, plot df['Frequency'], df['Flow'] as scatter (size=POINT_SIZE, color-coded).
             – Identify max‐flow point, add star via add_best_marker.
             – Title: “Flow vs Frequency | All raw points | Test <id> | s_f→e_f Hz, sweep_duration s | V V”.
             – Save under plots_2/all_points_{with|without}_err/frequency/<id>.png.
           – Time scatter: same but use df['Elapsed_s'], df['Flow'], Title “Flow vs Time | All raw points …”.
       • err flag only affects output folder path.

   – relative_plots(tests):
       • For each (folder, df):
           – Extract test_id, s_f, e_f, sweep_duration.
           – For each Sweep group: compute rel = 100 * (group['Flow'] / group['Flow'].max()).
             – Plot group['Frequency'] vs. rel as line+marker (size=POINT_SIZE).
             – Title: “Flow vs Frequency | Relative 0–100% | Test <id> | s_f→e_f Hz, sweep_duration s | V V”.
             – Save under plots_2/relative/frequency/<id>.png.

7. High-Flow / High-Stability Cloud
   – high_quality_cloud(stat_df, title, outfile):
       • Input: stat_df with index = frequency bins, columns ['mean','std'].
       • Drop rows with NaN mean or std.
       • top25 = rows where mean ≥ 75th percentile.
       • blue = subset of top25 where std ≤ 10th percentile of top25. If empty, use top25.
       • best_idx = the index of the row in blue with the smallest std; best_val = mean, best_err = std.
       • Plot all bins in light gray: error bars at (mean±std).
       • Plot “blue region” bins in blue error bars.
       • Plot star at (best_idx, best_val ± best_err) via add_best_marker.
       • X label = “Frequency [Hz]”, Y label = YLABEL, Title = title (passed in).
       • Save to outfile.

8. Data Summary
   – print_data_summary(tests):
       • Concatenate each df from tests into concat_raw.
       • Print “=== Data Summary ===”
         – Print “This dataset contains time-series measurements of pressure (if PRESSURE_TEST) or flow rate (otherwise).”
       • Variables & types: iterate over concat_raw.dtypes, print column name, dtype, and “quantitative” / “datetime” / “categorical” tag.
       • Total observations: len(concat_raw).
       • Total measurement time: sum over each test of (max(Elapsed_s) – min(Elapsed_s)), convert to minutes and seconds.
       • Missing values: count NaNs per column, print any non-zero counts or “No missing values found.”

9. Standard Deviation vs. Mean Flow
   – plot_std_vs_flow(df, out_path, title_prefix, point_size=3):
       • If raw DataFrame (has 'Frequency' column):
           – Compute df['bin'] = _bin(df['Frequency'], BIN_STEP_INDIVIDUAL).
           – ag = df.groupby('bin')['Flow'].agg(['mean','std']).reset_index().
         Else (pre-binned DataFrame with first column = bin, and 'mean','std'):
           – Rename first column to 'bin'; assume 'mean' and 'std' columns already exist.
       • Clean: drop rows where mean or std is NaN or std ≤ 0.
       • Compute Pearson r and Spearman ρ between ag['mean'] and ag['std'].
       • Compute linear regression slope, intercept (for plotting).
       • Plot:
           – Scatter ag['mean'] vs. ag['std'], colored by ag['bin'], sized by (point_size × 10).
           – Regression line over the full range of mean values.
           – Colorbar labeled “Frequency [Hz]”.
       • Labels: X = f"Mean {MEASUREMENT_NAME} ({UNIT})"; Y = f"Standard Deviation [{UNIT}]".
       • Title: f"{title_prefix} | Std vs Mean {MEASUREMENT_NAME} | Pearson r = {r:.3f}, Spearman ρ = {rho:.3f}".
       • Save to out_path.

10. Efficiency Calculation & Plotting
    – compute_efficiency_with_error(freq_avg):
        • Input freq_avg columns = ['Frequency','mean','std'] (units: µL/min, µL/min).
        • Convert:
            – Q_m3_s = (mean µL/min) × (1e-9 m³/µL) / 60 → m³/s.
            – Q_std_m3_s = likewise for 'std'.
            – v = Q_m3_s / AREA_M2 → m/s.
            – P_hyd = 0.5 × RHO × Q_m3_s × v^2 → W.
            – dP_dQ = 1.5 × RHO × (Q_m3_s^2) / (AREA_M2^2) → dP/dQ.
            – P_hyd_std = dP_dQ × Q_std_m3_s → W.
            – η = P_hyd / P_ELECTRIC (unitless).
            – η_std = P_hyd_std / P_ELECTRIC.
        • Build eff_df:
            – “Frequency”: same frequency bins.
            – “mean_flow”: original mean (µL/min).
            – “flow_std”: original std (µL/min).
            – “P_hyd”: computed hydraulic power (W).
            – “P_hyd_std”: uncertainty (W).
            – “η”: computed efficiency.
            – “η_std”: uncertainty.
        • Return eff_df.reset_index().

    – plot_efficiency(eff_df, out_dir, err=True|False, POINT_SIZE=3):
        • If err=True:  
            – ax.errorbar(eff_df["Frequency"], eff_df["η"], eff_df["η_std"], fmt='-o', capsize=2, markersize=POINT_SIZE, elinewidth=0.4, alpha=0.7, color='green', label="Efficiency")  
            – suffix = "with_err"
          Else:  
            – ax.plot(eff_df["Frequency"], eff_df["η"], '-o', markersize=POINT_SIZE, color='green', label="Efficiency")  
            – suffix = "no_err"
        • ax.set_xlabel("Frequency [Hz]"); ax.set_ylabel("Efficiency [η]"); ax.grid(True)
        • ax.set_title(f"Electrical Efficiency vs Frequency | {suffix.replace('_',' ')} | {SUPPLY_VOLTAGE} V")
        • ax.legend()
        • Save to os.path.join(out_dir, f"efficiency_{suffix}.png").

    – plot_efficiency_with_avg(eff_df, freq_avg, out_dir, err=True|False, POINT_SIZE=3, BEST_STAR_SCALE=5, SUPPLY_VOLTAGE=2.5):
        • merged = pd.merge(freq_avg, eff_df, on="Frequency", how="inner"). Contains ['Frequency','mean','std','η','η_std'].
        • idx_max_eff = merged["η"].idxmax(); freq_max_eff = merged.at[idx_max_eff,"Frequency"]; eta_max = merged.at[idx_max_eff,"η"]; eta_max_err = merged.at[idx_max_eff,"η_std"].
        • row483 = merged[merged["Frequency"] == 483]; if not empty, set freq_483, eta_483, eta_483_err.
      — Plot A (“with error bars”):
        1. fig, ax1 = plt.subplots(figsize=(10,6))
        2. ax1.errorbar(merged["Frequency"], merged["mean"], merged["std"], fmt='-o', capsize=2, markersize=POINT_SIZE, elinewidth=0.4, alpha=0.7, color='blue', label="Mean Flow Rate")
        3. ax2 = ax1.twinx()
        4. ax2.errorbar(merged["Frequency"], merged["η"], merged["η_std"], fmt='-o', capsize=2, markersize=POINT_SIZE, elinewidth=0.4, alpha=0.7, color='green', label="Efficiency")
        5. Highlight max-eff star:  
           ax2.errorbar(freq_max_eff, eta_max, eta_max_err, fmt="*", markersize=POINT_SIZE*BEST_STAR_SCALE, markerfacecolor="gold", markeredgecolor="red", markeredgewidth=1.4, capsize=2, elinewidth=0.4, alpha=0.9, color="red", label=f"Max η = {eta_max:.3e} ± {eta_max_err:.3e} @ {freq_max_eff:.0f} Hz")
        6. If 483 Hz present:  
           ax2.errorbar(freq_483, eta_483, eta_483_err, fmt="*", markersize=POINT_SIZE*BEST_STAR_SCALE, markerfacecolor="white", markeredgecolor="black", markeredgewidth=1.4, capsize=2, elinewidth=0.4, alpha=0.9, color="black", label=f"η = {eta_483:.3e} ± {eta_483_err:.3e} @ 483 Hz")
        7. Compute flow_min_err, flow_max_err = min(max(mean - std)), max(max(mean+std)). Add padding: bottom = flow_min_err - 0.10*(flow_max_err-flow_min_err), top = flow_max_err + 0.05*(flow_max_err-flow_min_err). ax1.set_ylim(bottom, top).
        8. f_flow = (0 - bottom)/(top - bottom).
        9. Compute eff_min_err, eff_max_err = min(min(η-η_std)), max(max(η+η_std)). Add padding: bottom_e = eff_min_err - 0.10*(eff_max_err-eff_min_err), top_e = eff_max_err + 0.05*(eff_max_err-eff_min_err).
       10. Align zero on efficiency axis:  
           eff_span = top_e - bottom_e;  
           candidate_bottom = 0 - f_flow*eff_span; candidate_top = candidate_bottom + eff_span;  
           Adjust candidate_bottom to be ≤ (eff_min_err - 0.10*(eff_max_err-eff_min_err)) if needed. Adjust candidate_top to be ≥ (eff_max_err + padding).  
           ax2.set_ylim(candidate_bottom, candidate_top).  
       11. Format ax2 tick labels in scientific notation “%.3e”.  
       12. Labels: ax1.set_xlabel("Frequency [Hz]"); ax1.set_ylabel("Flow Rate", color="blue"); ax2.set_ylabel("Efficiency [η]", color="green"); ax1.tick_params(axis="y", labelcolor="blue"); ax2.tick_params(axis="y", labelcolor="green").  
       13. Title: fig.suptitle(f"Flow Rate and Efficiency vs Frequency | with error bars | {SUPPLY_VOLTAGE} V").  
       14. Combine legends:  
           h1, l1 = ax1.get_legend_handles_labels();  
           h2, l2 = ax2.get_legend_handles_labels();  
           ax1.legend(h1 + h2, l1 + l2, fontsize="x-small", loc="upper left").  
       15. ax1.grid(True); plt.tight_layout(rect=[0,0,1,0.96]); savefig(out_dir+"/flow_eff_with_err.png", dpi=300).
      — Plot B (“points only”):
        1. fig, ax1 = plt.subplots(figsize=(10,6))  
        2. ax1.plot(merged["Frequency"], merged["mean"], "-o", markersize=POINT_SIZE, color="blue", label="Mean Flow Rate")  
        3. ax2 = ax1.twinx()  
        4. ax2.plot(merged["Frequency"], merged["η"], "-o", markersize=POINT_SIZE, color="green", label="Efficiency")  
        5. Highlight stars at max-eff and 483 Hz as above (with error bars on stars).  
        6. Compute new y-limits: flow_min_data = min(merged["mean"]), flow_max_data = max(merged["mean"]); pad_flow2 = 0.10*(flow_max_data - flow_min_data), flow_bottom2 = flow_min_data - pad_flow2, flow_top2 = flow_max_data + 0.05*(flow_max_data - flow_min_data); ax1.set_ylim(flow_bottom2, flow_top2);  
           f_flow2 = (0 - flow_bottom2)/(flow_top2 - flow_bottom2).  
        7. eff_min_data = min(merged["η"]), eff_max_data = max(merged["η"]); pad_eff2 = 0.10*(eff_max_data - eff_min_data), eff_bottom2 = eff_min_data - pad_eff2, eff_top2 = eff_max_data + pad_eff2;  
           eff_span2 = eff_top2 - eff_bottom2; candidate_bottom2 = 0 - f_flow2*eff_span2; candidate_top2 = candidate_bottom2 + eff_span2;  
           Adjust candidate_bottom2 ≤ eff_bottom2, adjust candidate_top2 ≥ eff_top2; ax2.set_ylim(candidate_bottom2, candidate_top2).  
        8. Format ax2 tick labels “%.3e”.  
        9. Labels: same as Plot A.  
       10. Title: fig.suptitle(f"Flow Rate and Efficiency vs Frequency | points only | {SUPPLY_VOLTAGE} V").  
       11. Combine legends as above.  
       12. ax1.grid(True); plt.tight_layout(rect=[0,0,1,0.96]); savefig(out_dir+"/flow_eff_points.png", dpi=300).

11. Standard Deviation vs. Mean Flow
    – plot_std_vs_flow calls (see above).

12. Exploratory Data Analysis
    – plot_individual_histograms(tests, base_out, unit):
        • Create folder base_out/individual_histograms.  
        • For each (folder, df):  
            – Compute a histogram of df['Flow'] (or df['Pressure'] if PRES_TEST=True), 30 bins, label axes and title “Histogram of Raw Flow – Test <id>”.  
            – Save to base_out/individual_histograms/test_<id>_histogram.png.  
    – plot_individual_boxplots(tests, base_out):
        • Create folder base_out/individual_boxplots.  
        • For each (folder, df):  
            – Use seaborn.boxplot to draw a boxplot of df['Flow'], label with YLABEL and “Boxplot of Raw Flow – Test <id>”.  
            – Save to base_out/individual_boxplots/test_<id>_boxplot.png.  
    – plot_average_histogram_from_raw(tests, base_out, unit, bin_step):
        • Create folder base_out/average_histogram.  
        • Concatenate each df[['Frequency','Flow']] into all_df. Compute all_df['bin'] = _bin(all_df['Frequency'], bin_step).  
        • Compute grouped = all_df.groupby('bin')['Flow'].mean().dropna() (mean flow per bin).  
        • Plot histogram of grouped (30 bins), labels “Mean Flow [unit]” vs. “Count”, title “Histogram of Binned-Average Values Across All Tests”.  
        • Save to base_out/average_histogram/average_histogram.png.  
    – plot_average_boxplot_from_raw(tests, base_out):
        • Create folder base_out/average_boxplot.  
        • Concatenate each df[['Flow']] into all_df (drop NaN).  
        • Draw boxplot of all_df['Flow'] (seaborn), label Y axis = YLABEL, title “Boxplot of Raw Flow Values Across All Tests”.  
        • Save to base_out/average_boxplot/average_boxplot.png.  
        • Compute summary stats via internal function _box_stats: Q1, median, Q3, whiskers, outlier count.  
        • Write CSV of stats_df (drop “Outliers” column) to base_out/average_boxplot/average_boxplot_stats.csv.  
        • Write human-readable stats to TXT at base_out/average_boxplot/average_boxplot_stats.txt.  
    – plot_raw_correlation(df_raw, base_out):
        • Create folder base_out/raw_correlation.  
        • Extract corr_df = df_raw[['Frequency','Elapsed_s','Flow']].dropna(). Compute corr = corr_df.corr().  
        • Plot seaborn.heatmap(corr, annot=True, cmap='coolwarm', square=True), title “Correlation Matrix: Raw Frequency, Time, and Flow”.  
        • Save to base_out/raw_correlation/raw_correlation_matrix.png.  
    – plot_aggregate_correlation(freq_avg, base_out):
        • Create folder base_out/aggregate_correlation.  
        • Compute corr_df = freq_avg[['Frequency','mean','std']].dropna(). Compute corr = corr_df.corr().  
        • Plot seaborn.heatmap(corr, annot=True, cmap='coolwarm', square=True), title “Correlation Matrix: Frequency, Mean, and Std”.  
        • Save to base_out/aggregate_correlation/aggregate_correlation_matrix.png.  
    – plot_combined_boxplots(tests, base_out, width=0.5, w_per_test=1.0):
        • Create folder base_out/combined_boxplots.  
        • Build a long‐form DataFrame combo with columns ['Test','Flow'], where each row is every flow value from each test, labeled by test_id.  
        • Determine order of tests by numeric extraction.  
        • Plot seaborn.boxplot(x='Test', y='Flow', data=combo, order=order, width=width). Rotate x‐ticks 45°. Y‐label = YLABEL, title = “All Individual-Test Boxplots”.  
        • Save to base_out/combined_boxplots/all_tests_boxplot.png.  
        • Compute per-test stats via _box_stats on each subset of combo, build stats_df.  
        • Save stats_df.drop('Outliers') to base_out/combined_boxplots/all_tests_stats.csv and TXT at all_tests_stats.txt.

13. Execution
    – At the end, print “Finished successfully.”

Usage:
    • Modify global parameters (parent_folder, BIN_STEP_*, PRESSURE_TEST, SUPPLY_VOLTAGE).
    • Run “python flow_and_pressure_analysis.py”.  
    • Inspect console for data summary and any “Skip” warnings if a CSV failed to load.  
    • Browse the “plots_2” directory for organized output subfolders, each containing the relevant PNGs (and CSV/TXT summaries).

By reading this docstring, a programmer or examiner can understand:
    1. The data input conventions (folder names, CSV format).
    2. How raw time-series flow/pressure is converted into frequency‐sweep data.
    3. Exactly which summary statistics and plot types are generated, and where to find them.
    4. The sequence and dependencies: CSV loading → binned statistics → per-test plots → combined plots → aggregate plots → efficiency → variability → EDA.
    5. The naming conventions and organization of output files.

End of flow_and_pressure_analysis.py docstring.
"""



import numpy as np

# ────────── Global parameters ────────────────────────────────────────────────
parent_folder              = '/Users/oskarjohnbruunsmith/Desktop/BP/Experimental_data/123_patch_test_pressure_test_02-05-2025/flow_rate_tests'
BIN_STEP_INDIVIDUAL        = 0   # 0=no bin, 1=every 1 Hz, 2=0,2,4 Hz, 5=0,5,10 Hz …
#BIN_STEP_INDIVIDUAL also regulates the bin for the efficiency plots.
#BIN_STEP_INDIVIDUAL also regulates the bin for the standard deviation vs flow rate plots
#No. that doesn't work. for std vs flow rate plots, step bin_ste_individual = 1 for consistant results.
BIN_STEP_COMBINED          = 5
BIN_STEP_AVERAGE           = 3
bin_step_hist_box_average  = 0 #The step for the histogram and boxplot average plots
POINT_SIZE                 = 3
BINS_RESID_HIST            = 20
PRESSURE_TEST              = False  # True → mbar ; False = µL/min
SUPPLY_VOLTAGE             = 2.5     # V
# How much larger the star should be than a normal point
BEST_STAR_SCALE = 5        # tweak until it looks right


P_ELECTRIC = 0.056201 #W
RHO = 997.0 #kg/m^3
DIAM_M = 1.4e-3 #m 1.4mm #
AREA_M2 = np.pi * (DIAM_M / 2)**2 # m^2


# ────────── Derived constants ────────────────────────────────────────────────
UNIT, MEASUREMENT_NAME = ("mbar", "Pressure") if PRESSURE_TEST else ("µL/min", "Flow Rate")
YLABEL = f"{MEASUREMENT_NAME} [{UNIT}]"

# ────────── Imports ─────────────────────────────────────────────────────────
import os, re, glob, warnings
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import pearsonr, spearmanr, linregress
import seaborn as sns
warnings.filterwarnings("ignore", category=FutureWarning)

# ────────── Helper utilities ────────────────────────────────────────────────
def ensure(path): os.makedirs(path, exist_ok=True); return path
def extract_test_info(folder):
    m = re.match(r'^(\d+)_', folder); return (m.group(1) if m else folder, *parse_sweep(folder))
def parse_sweep(folder):
    m = re.search(r'(\d+Hz_\d+Hz_\d+s)', folder)
    if not m: raise ValueError(f"No sweep info in {folder}")
    a, b, c = m.group(1).split('_'); return float(a[:-2]), float(b[:-2]), float(c[:-1])
def _bin(value, step): return value if step <= 0 else np.round(value / step) * step

def _clean_handles(handles):
    """drop artists whose label begins with '_' (Matplotlib 3.8 warning)"""
    return [h for h in handles if not h.get_label().startswith('_')]

def _legend_with_handles(ax, handles, **kw):
    """always supply handles as kw‑arg to avoid Matplotlib TypeError"""
    ax.legend(handles=_clean_handles(handles), **kw)

# ────────── Residual helper ────────────────────────────────────────────────
def residual_plots(resid, out, title):
    fig, ax = plt.subplots(2, 2, figsize=(10, 8)); fig.suptitle(title, y=.98)
    stats.probplot(resid, plot=ax[0, 0]); ax[0, 0].set_title("Normal Q‑Q")
    ax[0, 1].scatter(np.zeros_like(resid), resid, s=POINT_SIZE); ax[0, 1].axhline(0); ax[0, 1].set_title("Residuals vs fitted")
    ax[1, 0].hist(resid, bins=BINS_RESID_HIST); ax[1, 0].set_title("Histogram")
    ax[1, 1].scatter(range(len(resid)), resid, s=POINT_SIZE); ax[1, 1].axhline(0); ax[1, 1].set_title("Residuals vs run order")
    for a in ax.flat: a.grid(True); a.set_ylabel("Residuals")
    ax[0, 1].set_xlabel("Fitted = 0"); ax[1, 1].set_xlabel("Run order"); ax[1, 0].set_xlabel("Residuals")
    plt.tight_layout(); plt.savefig(out, dpi=300); plt.close()

# ────────── CSV loader with air‑bubble ±2 s filter ──────────────────────────
def load_and_process_csv(folder, s_f, e_f, sweep_s):
    trimmed = glob.glob(os.path.join(folder, 'trimmed_*.csv'))
    csvf    = trimmed[0] if trimmed else None
    if csvf: df = pd.read_csv(csvf, sep=';')
    if csvf and (df.empty or len(df) < 2): csvf = None
    if csvf is None:
        cand = [c for c in glob.glob(os.path.join(folder, '*.csv')) if not os.path.basename(c).startswith('trimmed_')]
        if not cand: raise FileNotFoundError(folder)
        csvf = cand[0]; df = pd.read_csv(csvf, sep=';')

    air_col  = next((c for c in df.columns if 'Air detected' in c), None)
    flow_col = (next((c for c in df.columns if re.search(r'IPS\s*\(\d+\)', c)), None)
                if PRESSURE_TEST else next((c for c in df.columns if 'Flowboard' in c), None))
    if not flow_col: raise ValueError(csvf)

    df = df[[df.columns[0], flow_col] + ([air_col] if air_col else [])]
    df.columns = ['Time', 'Flow'] + (['Air'] if air_col else [])
    df['Time'] = pd.to_datetime(df['Time'])
    df['Elapsed_s'] = (df['Time'] - df['Time'].iloc[0]).dt.total_seconds()

    if 'Air' in df.columns and df['Air'].any():
        t_vec = df['Elapsed_s'].values
        bad = np.zeros(len(df), bool)
        for t in df.loc[df['Air'].astype(bool), 'Elapsed_s']:
            bad |= np.abs(t_vec - t) <= 2.0
        df = df.loc[~bad]

    span = e_f - s_f
    df['Frequency'] = df['Elapsed_s'].apply(lambda t: s_f + span * (((t) % sweep_s) / sweep_s))
    df['Sweep'] = np.floor((df['Elapsed_s']) / sweep_s).astype(int)
    return df[['Time', 'Flow', 'Elapsed_s', 'Frequency', 'Sweep']]

# ────────── Plot helpers ────────────────────────────────────────────────────

def _plot(ax, x, y, std=None, **kw):
    if std is None:
        ax.plot(x, y, '-o', markersize=POINT_SIZE, **kw)
    else:
        ax.errorbar(
            x, y, std,
            fmt='o',
            capsize=2,           # shorter caps
            markersize=POINT_SIZE,
            elinewidth=0.4,      # thinner vertical line
            alpha=0.7,           # a hint of transparency
            **kw
        )


def _label(ax): ax.set_ylabel(YLABEL); ax.grid(True)

def add_best_marker(ax, x_best, y_best, std_best=None):
    """
    Draw the best-value point as a clearly visible star and return the artist.

    - Size = POINT_SIZE * BEST_STAR_SCALE  (separate from normal points)
    - Filled golden star with a red outline so it pops against any background
    - High z-order so nothing can hide it
    """
    xb = float(np.ravel(x_best)[0])
    yb = float(np.ravel(y_best)[0])

    label = f"★ {yb:.1f} {UNIT}"
    if std_best is not None:
        label += f" ±{abs(float(std_best)):.1f}"
    label += f" @ {xb:.0f} Hz"

    star, = ax.plot(
        xb,
        yb,
        marker='*',
        markersize=POINT_SIZE * BEST_STAR_SCALE,
        markerfacecolor='gold',
        markeredgecolor='red',
        markeredgewidth=1.4,
        linestyle='',
        label=label,
        zorder=10          # sits on top of every other artist
    )
    return star


# ────────── Output root ─────────────────────────────────────────────────────
plots2 = ensure(os.path.join(parent_folder, "plots_2"))
def subdir(*parts): return ensure(os.path.join(plots2, *parts))

# ────────── Individual test plots ──────────────────────────────────────────
def indiv_plot_freq(tests, out_dir, err):
    """
    Per-test Flow-vs-Frequency plot.
    Legend always shows one coloured marker per sweep (even with error bars).
    """
    cmap, step = plt.get_cmap('tab10'), BIN_STEP_INDIVIDUAL

    for folder, df in tests:
        tn, s_f, e_f, sw_s = extract_test_info(os.path.basename(folder))
        fig, ax = plt.subplots(figsize=(8, 5))

        legend_handles = []

        # one trace per sweep ------------------------------------------------
        for i, g in df.groupby('Sweep'):
            color = cmap(i % 10)
            g['bin'] = _bin(g['Frequency'], step)
            ag = g.groupby('bin')['Flow'].agg(['mean', 'std']).reset_index()

            # draw the data
            _plot(
                ax,
                ag['bin'],
                ag['mean'],
                ag['std'] if err else None,
                color=color,
                label='_nolegend_'      # avoid duplicate small markers
            )

            # build a big marker proxy for the legend
            proxy = plt.Line2D(
                [], [], marker='o', linestyle='',
                color=color,
                markersize=POINT_SIZE * 3,
                label=f"Sweep {i + 1}"
            )
            legend_handles.append(proxy)

        # overall best value -------------------------------------------------
        tmp = df.copy()
        tmp['bin'] = _bin(tmp['Frequency'], step)
        gb = tmp.groupby('bin')['Flow'].agg(['mean', 'std'])
        star = add_best_marker(
            ax,
            gb['mean'].idxmax(),
            gb['mean'].max(),
            gb.loc[gb['mean'].idxmax(), 'std'] if err else None
        )
        legend_handles.append(star)

        # legend & cosmetics -------------------------------------------------
        _legend_with_handles(ax, legend_handles, fontsize='xx-small', ncol=3)
        ax.set_xlabel("Frequency [Hz]")
        _label(ax)
        ax.set_title(
            f"{MEASUREMENT_NAME} vs Frequency  |  Test {tn}  |  "
            f"{s_f:.1f}→{e_f:.1f} Hz, {sw_s:.1f}s sweep  |  {SUPPLY_VOLTAGE} V"
        )

        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, f"test_{tn}.png"), dpi=300)
        plt.close()




def indiv_plot_time(tests, out_dir, err):
    step = BIN_STEP_INDIVIDUAL
    for folder, df in tests:
        tn, s_f, e_f, sw_s = extract_test_info(os.path.basename(folder))
        g = df.copy()
        g['bin'] = _bin(g['Elapsed_s'], step)
        ag = g.groupby('bin')['Flow'].agg(['mean', 'std']).reset_index()

        fig, ax = plt.subplots(figsize=(8, 5))
        _plot(ax, ag['bin'], ag['mean'], ag['std'] if err else None)

        idx_best = ag['mean'].idxmax()
        add_best_marker(
            ax,
            ag.loc[idx_best, 'bin'],
            ag.loc[idx_best, 'mean'],
            ag.loc[idx_best, 'std'] if err else None
        )

        _legend_with_handles(ax, ax.get_lines(), fontsize='x-small')
        ax.set_xlabel("Elapsed Time [s]")
        _label(ax)
        ax.set_title(
            f"{MEASUREMENT_NAME} vs Time  |  Test {tn}  |  "
            f"{s_f:.1f}→{e_f:.1f} Hz, {sw_s:.1f}s sweep  |  {SUPPLY_VOLTAGE} V"
        )
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, f"test_{tn}.png"), dpi=300)
        plt.close()

# ────────── Combined series & Average plots ────────────────────────────────
def comb_series(tests, out_dir, vs_freq, error_bars=False):
    """
    Combined series of Flow vs Frequency (or Time).
    Legend now shows large, coloured symbols (with an error bar when appropriate).
    """
    step = BIN_STEP_COMBINED
    key  = 'Frequency' if vs_freq else 'Elapsed_s'
    fig, ax = plt.subplots(figsize=(10, 6))

    handles = []                      # proxy handles for legend only

    for fpath, df in tests:
        tn, *_ = extract_test_info(os.path.basename(fpath))
        df['bin'] = _bin(df[key], step)

        grp_mean = df.groupby('bin')['Flow'].mean()
        grp_std  = df.groupby('bin')['Flow'].std()

        if error_bars:
            # draw the data + error bars
            cont = ax.errorbar(
                grp_mean.index, grp_mean.values, grp_std.values,
                fmt='o', ms=POINT_SIZE, capsize=3,
                label='_nolegend_'
            )
            # build a big proxy error-bar symbol for the legend
            leg_proxy = ax.errorbar(
                [], [], yerr=[1],           # dummy data
                fmt='o',
                color=cont[0].get_color(),
                ms=POINT_SIZE * 3,
                capsize=4,
                label=f"Test {tn}"
            )
            handles.append(leg_proxy)
        else:
            line, = ax.plot(
                grp_mean.index, grp_mean.values,
                linestyle='', marker='o', ms=POINT_SIZE,
                label='_nolegend_'
            )
            # big coloured dot for legend
            proxy = plt.Line2D(
                [], [], marker='o', linestyle='',
                color=line.get_color(),
                markersize=POINT_SIZE * 3,
                label=f"Test {tn}"
            )
            handles.append(proxy)

    # global-maximum star (already large and bright)
    big = pd.concat([d[[key, 'Flow']] for _, d in tests])
    star = add_best_marker(
        ax,
        big.loc[big['Flow'].idxmax(), key],
        big['Flow'].max()
    )
    handles.append(star)

    ax.set_xlabel("Frequency [Hz]" if vs_freq else "Elapsed Time [s]")
    _label(ax)

    what   = "Frequency" if vs_freq else "Time"
    suffix = "with error bars" if error_bars else "no error bars"
    ax.set_title(
        f"{MEASUREMENT_NAME} vs {what}  |  Combined Points  |  "
        f"{SUPPLY_VOLTAGE} V  ({suffix})"
    )

    _legend_with_handles(ax, handles, fontsize='x-small', ncol=3)
    plt.tight_layout()
    fname = f"combined_series_{'err' if error_bars else 'noerr'}.png"
    plt.savefig(os.path.join(out_dir, fname), dpi=300)
    plt.close()





def avg_plot(tests, out_dir, vs_freq, err):
    step, key = BIN_STEP_AVERAGE, ('Frequency' if vs_freq else 'Elapsed_s')
    all_df = pd.concat([d[[key,'Flow']] for _,d in tests]); all_df['bin'] = _bin(all_df[key], step)
    ag = all_df.groupby('bin')['Flow'].agg(['mean','std']).reset_index()
    fig, ax = plt.subplots(figsize=(10,6))
    _plot(ax, ag['bin'], ag['mean'], ag['std'] if err else None)
    idx_best = ag['mean'].idxmax()
    star = add_best_marker(ax, ag.loc[idx_best,'bin'], ag.loc[idx_best,'mean'],
                           ag.loc[idx_best,'std'] if err else None)
    lab = "Frequency (Hz)" if vs_freq else "Elapsed Time (s)"
    ax.set_xlabel(lab); _label(ax)
    what = "Frequency" if vs_freq else "Time"
    eb = "with error bars" if err else "no error bars"
    ax.set_title(f"{MEASUREMENT_NAME} vs {what}  |  Average of all tests  |  {eb}  |  {SUPPLY_VOLTAGE} V")
    _legend_with_handles(ax, [star], fontsize='x-small')
    plt.tight_layout(); plt.savefig(os.path.join(out_dir,"average.png"), dpi=300); plt.close()
    
    
    
# ────────── Data description ─────────────────────────────────────────────────
def print_data_summary(tests):
    """
    Print:
      1. Short description
      2. Variables in the DataFrame
      3. Their types
      4. Total number of observations (all tests)
      5. Total measurement time across all tests
      6. Any missing‐value counts
    """
    import pandas as pd
    print("=== Data Summary ===")
    print("This dataset contains time‐series measurements of",
          "pressure" if PRESSURE_TEST else "flow rate",
          "collected during frequency sweeps.\n")

    # Build a single concat just for variable listing & missing‐value check:
    concat_raw = pd.concat([df for _,df in tests], ignore_index=True)

    # 2+3) Variables & types
    print("Variables:")
    for col, dtype in concat_raw.dtypes.items():
        if pd.api.types.is_datetime64_any_dtype(dtype):
            vtype = "datetime"
        elif pd.api.types.is_numeric_dtype(dtype):
            vtype = "quantitative"
        else:
            vtype = "categorical"
        print(f"  - {col:12s} : {dtype.name:10s} ({vtype})")
    print()

    # 4) Total observations
    total_obs = len(concat_raw)
    print(f"Total observations (all tests): {total_obs}")

    # 5) Total measurement time
    total_seconds = 0.0
    for path, df in tests:
        if 'Elapsed_s' in df:
            total_seconds += df['Elapsed_s'].max() - df['Elapsed_s'].min()
    mins, secs = divmod(total_seconds, 60)
    print(f"Total measurement time: {int(mins)} minutes, {int(secs)} seconds\n")

    # 6) Missing values
    na = concat_raw.isna().sum()
    na = na[na>0]
    if na.empty:
        print("No missing values found.")
    else:
        print("Missing values:")
        for col,count in na.items():
            print(f"  - {col}: {count}")
    print("====================\n")



# ────────── Exploratory plots ──────────────────────────────────────────────

def plot_individual_histograms(tests, base_out, unit):
    """
    Create and save raw-data histograms for each individual test.
    """
    folder = ensure(os.path.join(base_out, 'individual_histograms'))
    meas = 'Flow' if PRESSURE_TEST else 'Flow'
    for path, df in tests:
        tn, *_ = extract_test_info(os.path.basename(path))
        plt.figure(figsize=(8, 5))
        plt.hist(df[meas], bins=30, edgecolor='black', alpha=0.7)
        plt.xlabel(f"{meas} ({unit})")
        plt.ylabel("Count")
        plt.title(f"Histogram of Raw {meas} – Test {tn}")
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(os.path.join(folder, f"test_{tn}_histogram.png"), dpi=300)
        plt.close()


def plot_individual_boxplots(tests, base_out):
    """
    Create and save raw‐data boxplots for each individual test,
    always using df['Flow'] and the global YLABEL on the y‐axis.
    """
    folder = ensure(os.path.join(base_out, 'individual_boxplots'))

    for path, df in tests:
        tn, *_ = extract_test_info(os.path.basename(path))

        plt.figure(figsize=(6, 5))
        sns.boxplot(y=df['Flow'], width=0.4)

        # Use the global YLABEL (already set after Derived constants)
        plt.ylabel(YLABEL)
        plt.title(f"Boxplot of Raw {MEASUREMENT_NAME} – Test {tn}")
        plt.grid(True)
        plt.tight_layout()

        plt.savefig(os.path.join(folder, f"test_{tn}_boxplot.png"), dpi=300)
        plt.close()


def plot_average_histogram_from_raw(tests, base_out, unit, bin_step):
    folder = ensure(os.path.join(base_out, 'average_histogram'))

    # Extract all (frequency, flow) data
    all_df = pd.concat([df[['Frequency', 'Flow']] for _, df in tests])
    all_df['bin'] = _bin(all_df['Frequency'], bin_step)

    # Compute bin-wise mean flow
    grouped = all_df.groupby('bin')['Flow'].mean().dropna()

    # Plot histogram of binned means
    plt.figure(figsize=(8, 5))
    plt.hist(grouped, bins=30, edgecolor='black', alpha=0.7)
    plt.xlabel(f"Mean {('Pressure' if PRESSURE_TEST else 'Flow')} [{unit}]")
    plt.ylabel("Count")
    plt.title("Histogram of Binned-Average Values Across All Tests")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(folder, "average_histogram.png"), dpi=300)
    plt.close()



def plot_average_boxplot_from_raw(tests, base_out):
    """
    Create and save a boxplot of raw ‘Flow’/‘Pressure’ values across all tests,
    always using df['Flow'] and the global YLABEL on the y‐axis.
    Also saves full statistics to CSV/TXT.
    """
    folder = ensure(os.path.join(base_out, 'average_boxplot'))

    # Combine all raw “Flow” values from all tests into one DataFrame
    collected = []
    for _, df in tests:
        collected.append(df[['Flow']])
    all_df = pd.concat(collected).dropna()

    # --- Plot figure ---
    plt.figure(figsize=(6, 5))
    sns.boxplot(y=all_df['Flow'], width=0.4)

    # Use the global YLABEL
    plt.ylabel(YLABEL)
    plt.title(f"Boxplot of Raw {MEASUREMENT_NAME} Values Across All Tests")
    plt.grid(True)
    plt.tight_layout()

    plt.savefig(os.path.join(folder, "average_boxplot.png"), dpi=300)
    plt.close()

    # --- Compute statistics (Q1, median, Q3, whiskers, outliers) ---
    stats = pd.Series(_box_stats(all_df['Flow']), name="AllTests")
    stats.drop(labels=['Outliers']).to_frame().T.to_csv(
        os.path.join(folder, "average_boxplot_stats.csv"), index=False
    )

    with open(os.path.join(folder, "average_boxplot_stats.txt"), 'w') as f:
        f.write(stats.to_string())
        

def plot_raw_correlation(df_raw, base_out):
    """
    Create and save correlation heatmap for raw data: Frequency, Elapsed_s, Flow/Pressure.
    """
    folder = ensure(os.path.join(base_out, 'raw_correlation'))
    meas = 'Flow' if PRESSURE_TEST else 'Flow'
    corr_df = df_raw[['Frequency', 'Elapsed_s', meas]].dropna()
    corr = corr_df.corr()
    plt.figure(figsize=(6, 5))
    sns.heatmap(corr, annot=True, cmap='coolwarm', fmt='.2f', square=True)
    plt.title("Correlation Matrix: Raw Frequency, Time, and %s" % meas)
    plt.tight_layout()
    out_path = os.path.join(folder, "raw_correlation_matrix.png")
    plt.savefig(out_path, dpi=300)
    plt.close()


def plot_aggregate_correlation(freq_avg, base_out):
    """
    Create and save correlation heatmap for aggregated binned data: Frequency, mean, std.
    Expects freq_avg with columns ['Frequency','mean','std']
    """
    folder = ensure(os.path.join(base_out, 'aggregate_correlation'))
    corr_df = freq_avg[['Frequency', 'mean', 'std']].dropna()
    corr = corr_df.corr()
    plt.figure(figsize=(6, 5))
    sns.heatmap(corr, annot=True, cmap='coolwarm', fmt='.2f', square=True)
    plt.title("Correlation Matrix: Frequency, Mean, and Std")
    plt.tight_layout()
    out_path = os.path.join(folder, "aggregate_correlation_matrix.png")
    plt.savefig(out_path, dpi=300)
    plt.close()

def _box_stats(s):
    """Return a dict with full box-plot stats for a pandas Series."""
    q1 = s.quantile(0.25)
    q3 = s.quantile(0.75)
    median = s.median()
    iqr = q3 - q1
    low_whisk  = s[s >= q1 - 1.5 * iqr].min()
    high_whisk = s[s <= q3 + 1.5 * iqr].max()
    outliers   = s[(s < low_whisk) | (s > high_whisk)].values
    return {
        "Q1": q1, "Median": median, "Q3": q3,
        "LowerWhisker": low_whisk, "UpperWhisker": high_whisk,
        "OutlierCount": len(outliers),
        "Outliers": outliers
    }


def plot_combined_boxplots(tests, base_out, width=0.5, w_per_test=1.0):
    """
    Plot side‐by‐side boxplots for every individual test, using df['Flow']
    on the y‐axis and the global YLABEL. Also saves statistics (Q1, median,
    Q3, whiskers, outliers) per test.
    """
    folder = ensure(os.path.join(base_out, 'combined_boxplots'))

    # Build a long‐form DataFrame with two columns: 'Test' and 'Flow'
    rows = []
    for path, df in tests:
        tn, *_ = extract_test_info(os.path.basename(path))
        label = f"Test {tn}"
        for v in df['Flow'].dropna():
            rows.append((label, v))
    combo = pd.DataFrame(rows, columns=['Test', 'Flow'])

    # Order tests numerically
    combo['Test_num'] = combo['Test'].str.extract(r'(\d+)').astype(int)
    order = combo.sort_values('Test_num')['Test'].unique()

    # --- Figure ---
    n_tests = len(order)
    plt.figure(figsize=(max(8, w_per_test * n_tests), 6))
    sns.boxplot(
        x='Test',
        y='Flow',
        data=combo,
        width=width,
        order=order
    )
    plt.xticks(rotation=45, ha='right')

    # Use the global YLABEL
    plt.ylabel(YLABEL)
    plt.title("All Individual‐Test Boxplots")
    plt.grid(True, axis='y')
    plt.tight_layout()

    plt.savefig(os.path.join(folder, "all_tests_boxplot.png"), dpi=300)
    plt.close()

    # --- Statistics per test ---
    stats_rows = []
    for lbl in order:
        vals = combo.loc[combo['Test'] == lbl, 'Flow']
        stats_rows.append(pd.Series(_box_stats(vals), name=lbl))
    stats_df = pd.DataFrame(stats_rows)
    stats_df.index.name = 'Test'

    # Save a CSV (dropping the actual outlier values) and a human‐readable TXT
    stats_df.drop(columns=['Outliers']).to_csv(
        os.path.join(folder, "all_tests_stats.csv")
    )
    with open(os.path.join(folder, "all_tests_stats.txt"), 'w') as f:
        f.write(stats_df.to_string())




# ────────── Standard deviation vs flow rate plots ────────────────────────────

def plot_std_vs_flow(df, out_path, title_prefix, point_size=3):
    df = df.copy()

    is_raw = 'Frequency' in df.columns

    # Choose measurement column based on mode
    measurement_col = 'Pressure' if PRESSURE_TEST else 'Flow'

    if is_raw:
        if measurement_col not in df.columns:
            raise ValueError(f"Expected column '{measurement_col}' not found in raw DataFrame.")
        df['bin'] = _bin(df['Frequency'], BIN_STEP_INDIVIDUAL)
        ag = df.groupby('bin')[measurement_col].agg(['mean', 'std']).reset_index()
    else:
        # Already binned – must include 'mean' and 'std' columns
        if not {'mean', 'std'}.issubset(df.columns):
            raise ValueError("Pre-binned DataFrame must include 'mean' and 'std' columns.")
        ag = df.rename(columns={df.columns[0]: 'bin'})  # bin column is first

    # Clean
    clean_df = ag.dropna(subset=['mean', 'std'])
    clean_df = clean_df[clean_df['std'] > 0]

    if clean_df.empty:
        print(f"⚠️ Warning: No valid data to plot for {title_prefix}")
        return

    from scipy.stats import pearsonr, spearmanr, linregress
    import numpy as np
    import matplotlib.pyplot as plt

    pearson_corr, _ = pearsonr(clean_df['mean'], clean_df['std'])
    spearman_corr, _ = spearmanr(clean_df['mean'], clean_df['std'])
    slope, intercept, *_ = linregress(clean_df['mean'], clean_df['std'])
    x_vals = np.linspace(clean_df['mean'].min(), clean_df['mean'].max(), 200)
    y_vals = slope * x_vals + intercept

    fig, ax = plt.subplots(figsize=(10, 6))
    sc = ax.scatter(
        clean_df['mean'],
        clean_df['std'],
        c=clean_df['bin'],
        cmap='viridis',
        s=point_size * 10,  # scale marker size for visibility
        label='Data points'
    )
    plt.colorbar(sc, ax=ax, label='Frequency [Hz]')
    ax.plot(x_vals, y_vals, '-', color='red', linewidth=1.5, label='Linear regression')

    ax.set_xlabel(f"Mean {MEASUREMENT_NAME} ({UNIT})")
    ax.set_ylabel(f"Standard Deviation [{UNIT}]")
    ax.grid(True)
    ax.set_title(
        f"{title_prefix} | Std vs Mean {MEASUREMENT_NAME}  | "
        f"Pearson r = {pearson_corr:.3f}, Spearman ρ = {spearman_corr:.3f}",
        fontsize=11
    )
    ax.legend()
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()


    
# ────────── Efficiency function and plots ────────────────────────────────────


def compute_efficiency_with_error(freq_avg: pd.DataFrame) -> pd.DataFrame:
    df = freq_avg.copy()
    df['bin'] = _bin(df['Frequency'], BIN_STEP_INDIVIDUAL)
    df = df.groupby('bin')[['mean', 'std']].mean().reset_index().rename(columns={'bin': 'Frequency'})
    df = df.set_index('Frequency')

    Q_m3_s = df['mean'] * 1e-9 / 60 #To go from mucriliter / min to m^3 / s.
    Q_std_m3_s = df['std'] * 1e-9 / 60 ##To go from mucriliter / min to m^3 / s.
    v = Q_m3_s / AREA_M2
    P_hyd = 0.5 * RHO * Q_m3_s * v**2
    dP_dQ = 1.5 * RHO * (Q_m3_s**2) / (AREA_M2**2)
    P_hyd_std = dP_dQ * Q_std_m3_s
    eta = P_hyd / P_ELECTRIC
    eta_std = P_hyd_std / P_ELECTRIC

    eff_df = pd.DataFrame({
        "mean_flow": df["mean"],
        "flow_std": df["std"],
        "P_hyd": P_hyd,
        "P_hyd_std": P_hyd_std,
        "η": eta,
        "η_std": eta_std
    })
    eff_df.index.name = "Frequency"
    return eff_df.reset_index()


def plot_efficiency(eff_df, out_dir, err=True, POINT_SIZE=3, SUPPLY_VOLTAGE=SUPPLY_VOLTAGE):
    fig, ax = plt.subplots(figsize=(10, 6))
    if err:
        # Thin error‐bar style: capsize=2, elinewidth=0.4, alpha=0.7
        ax.errorbar(
            eff_df["Frequency"],
            eff_df["η"],
            eff_df["η_std"],
            fmt='-o',
            capsize=2,
            markersize=POINT_SIZE,
            elinewidth=0.4,
            alpha=0.7,
            color='green',
            label="Efficiency"
        )
        suffix = "with_err"
    else:
        ax.plot(
            eff_df["Frequency"],
            eff_df["η"],
            '-o',
            markersize=POINT_SIZE,
            color='green',
            label="Efficiency"
        )
        suffix = "no_err"

    ax.set_xlabel("Frequency [Hz]")
    ax.set_ylabel("Efficiency [η]")
    ax.grid(True)
    ax.set_title(
        f"Electrical Efficiency vs Frequency  |  {suffix.replace('_', ' ')}  |  {SUPPLY_VOLTAGE} V",
        fontsize=11
    )
    ax.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f"efficiency_{suffix}.png"), dpi=300)
    plt.close()






def plot_efficiency_with_avg(
    eff_df,
    freq_avg,
    out_dir,
    err=True,               # retained for compatibility; ignored internally
    POINT_SIZE=3,
    BEST_STAR_SCALE=5,
    SUPPLY_VOLTAGE=2.5
):
    """
    Produces two versions of the "flow rate + efficiency vs. frequency" plot,
    highlighting both the maximum‐efficiency point (with η ± uncertainty in legend)
    and the efficiency at 483 Hz (with η ± uncertainty). Both highlighted points
    include error bars, and the left‐axis y‐label is always "Flow Rate". Titles
    explicitly say "Flow Rate and Efficiency vs Frequency". The efficiency‐axis
    ticks are formatted in scientific notation with three significant digits,
    ensuring at least three non‐zero digits appear.

      1) flow_eff_with_err.png   → error‐bar version
      2) flow_eff_points.png     → point‐only version

    The `err` parameter is ignored internally—both plots are always generated.
    """

    from matplotlib.ticker import FormatStrFormatter

    # Merge flow‐rate and efficiency data on Frequency
    #   - freq_avg must have columns ["Frequency","mean","std"]
    #   - eff_df must have ["Frequency","η","η_std"]
    merged = pd.merge(freq_avg, eff_df, on="Frequency", how="inner")

    os.makedirs(out_dir, exist_ok=True)

    # Identify the maximum‐efficiency row
    idx_max_eff = merged["η"].idxmax()
    freq_max_eff = merged.at[idx_max_eff, "Frequency"]
    eta_max = merged.at[idx_max_eff, "η"]
    eta_max_err = merged.at[idx_max_eff, "η_std"]

    # Identify the 483 Hz row (if it exists exactly)
    row483 = merged[merged["Frequency"] == 483]
    has_483 = not row483.empty
    if has_483:
        freq_483 = float(row483["Frequency"].iloc[0])
        eta_483 = float(row483["η"].iloc[0])
        eta_483_err = float(row483["η_std"].iloc[0])

    # ────────── Plot 1: with error bars ─────────────────────────────────────────
    fig, ax1 = plt.subplots(figsize=(10, 6))

    # Flow Rate (with error bars)
    ax1.errorbar(
        merged["Frequency"],
        merged["mean"],
        merged["std"],
        fmt="-o",
        capsize=2,
        markersize=POINT_SIZE,
        elinewidth=0.4,
        alpha=0.7,
        color="blue",
        label="Mean Flow Rate"
    )

    # Efficiency (with error bars) on twin axis
    ax2 = ax1.twinx()
    ax2.errorbar(
        merged["Frequency"],
        merged["η"],
        merged["η_std"],
        fmt="-o",
        capsize=2,
        markersize=POINT_SIZE,
        elinewidth=0.4,
        alpha=0.7,
        color="green",
        label="Efficiency"
    )

    # Highlight the maximum‐efficiency point (star marker + η ± uncertainty)
    ax2.errorbar(
        freq_max_eff,
        eta_max,
        eta_max_err,
        fmt="*",
        markersize=POINT_SIZE * BEST_STAR_SCALE,
        markerfacecolor="gold",
        markeredgecolor="red",
        markeredgewidth=1.4,
        capsize=2,
        elinewidth=0.4,
        alpha=0.9,
        color="red",
        label=f"Max η = {eta_max:.3e} ± {eta_max_err:.3e} @ {freq_max_eff:.0f} Hz"
    )

    # Highlight the 483 Hz efficiency point (if present, with η ± uncertainty)
    if has_483:
        ax2.errorbar(
            freq_483,
            eta_483,
            eta_483_err,
            fmt="*",
            markersize=POINT_SIZE * BEST_STAR_SCALE,
            markerfacecolor="white",
            markeredgecolor="black",
            markeredgewidth=1.4,
            capsize=2,
            elinewidth=0.4,
            alpha=0.9,
            color="black",
            label=f"η = {eta_483:.3e} ± {eta_483_err:.3e} @ 483 Hz"
        )

    # Determine y‐limits to cover error bars exactly, with padding

    # --- Flow Rate axis: min/max of (mean ± std)
    flow_min_err = (merged["mean"] - merged["std"]).min()
    flow_max_err = (merged["mean"] + merged["std"]).max()
    flow_err_span = flow_max_err - flow_min_err
    pad_flow_bottom = flow_err_span * 0.10 if flow_err_span != 0 else abs(flow_max_err) * 0.10
    pad_flow_top = flow_err_span * 0.05 if flow_err_span != 0 else abs(flow_max_err) * 0.05
    flow_bottom = flow_min_err - pad_flow_bottom
    flow_top = flow_max_err + pad_flow_top
    ax1.set_ylim(flow_bottom, flow_top)

    # Compute fraction of flow‐axis height at which zero sits
    f_flow = (0 - flow_bottom) / (flow_top - flow_bottom)

    # --- Efficiency axis: min/max of (η ± η_std)
    eff_min_err = (merged["η"] - merged["η_std"]).min()
    eff_max_err = (merged["η"] + merged["η_std"]).max()
    eff_err_span = eff_max_err - eff_min_err
    pad_eff_bottom = eff_err_span * 0.10 if eff_err_span != 0 else abs(eff_max_err) * 0.10
    pad_eff_top = eff_err_span * 0.05 if eff_err_span != 0 else abs(eff_max_err) * 0.05
    eff_bottom = eff_min_err - pad_eff_bottom
    eff_top = eff_max_err + pad_eff_top

    # Align zero on the efficiency axis to the same fractional height f_flow
    eff_span = eff_top - eff_bottom
    candidate_bottom = 0 - f_flow * eff_span
    candidate_top = candidate_bottom + eff_span

    # Ensure bottom covers actual minimum error bar
    desired_bottom = eff_min_err - pad_eff_bottom
    if candidate_bottom > desired_bottom:
        shift = candidate_bottom - desired_bottom
        candidate_bottom -= shift
        candidate_top  -= shift

    # Ensure top covers actual maximum error bar
    desired_top = eff_max_err + pad_eff_top
    if candidate_top < desired_top:
        shift = desired_top - candidate_top
        candidate_top  += shift
        candidate_bottom += shift

    ax2.set_ylim(candidate_bottom, candidate_top)

    # Format efficiency y‐axis in scientific notation with three significant digits
    ax2.yaxis.set_major_formatter(FormatStrFormatter('%.3e'))

    # Labels, title, grid, and legend
    ax1.set_xlabel("Frequency [Hz]")
    ax1.set_ylabel("Flow Rate", color="blue")
    ax1.tick_params(axis="y", labelcolor="blue")
    ax2.set_ylabel("Efficiency [η]", color="green")
    ax2.tick_params(axis="y", labelcolor="green")

    fig.suptitle(
        f"Flow Rate and Efficiency vs Frequency  |  with error bars  |  {SUPPLY_VOLTAGE} V",
        fontsize=11
    )

    # Combine legends from both axes
    handles1, labels1 = ax1.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(handles1 + handles2, labels1 + labels2, fontsize="x-small", loc="upper left")

    ax1.grid(True)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(os.path.join(out_dir, "flow_eff_with_err.png"), dpi=300)
    plt.close()


    # ────────── Plot 2: points/lines only (no error bars on the main traces) ───
    fig, ax1 = plt.subplots(figsize=(10, 6))

    # Flow Rate (points only)
    ax1.plot(
        merged["Frequency"],
        merged["mean"],
        "-o",
        markersize=POINT_SIZE,
        color="blue",
        label="Mean Flow Rate"
    )

    # Efficiency (points only) on twin axis
    ax2 = ax1.twinx()
    ax2.plot(
        merged["Frequency"],
        merged["η"],
        "-o",
        markersize=POINT_SIZE,
        color="green",
        label="Efficiency"
    )

    # Highlight the maximum‐efficiency point (star marker + η ± uncertainty)
    ax2.errorbar(
        freq_max_eff,
        eta_max,
        eta_max_err,
        fmt="*",
        markersize=POINT_SIZE * BEST_STAR_SCALE,
        markerfacecolor="gold",
        markeredgecolor="red",
        markeredgewidth=1.4,
        capsize=2,
        elinewidth=0.4,
        alpha=0.9,
        color="red",
        label=f"Max η = {eta_max:.3e} ± {eta_max_err:.3e} @ {freq_max_eff:.0f} Hz"
    )

    # Highlight the 483 Hz efficiency point (if present, with η ± uncertainty)
    if has_483:
        ax2.errorbar(
            freq_483,
            eta_483,
            eta_483_err,
            fmt="*",
            markersize=POINT_SIZE * BEST_STAR_SCALE,
            markerfacecolor="white",
            markeredgecolor="black",
            markeredgewidth=1.4,
            capsize=2,
            elinewidth=0.4,
            alpha=0.9,
            color="black",
            label=f"η = {eta_483:.3e} ± {eta_483_err:.3e} @ 483 Hz"
        )

    # Determine y‐limits to frame min/max of data with padding

    # --- Flow Rate axis for points ---
    flow_min_data = merged["mean"].min()
    flow_max_data = merged["mean"].max()
    flow_data_span = flow_max_data - flow_min_data
    pad_flow2 = flow_data_span * 0.10 if flow_data_span != 0 else abs(flow_max_data) * 0.10
    flow_bottom2 = flow_min_data - pad_flow2
    flow_top2 = flow_max_data + pad_flow2
    ax1.set_ylim(flow_bottom2, flow_top2)

    f_flow2 = (0 - flow_bottom2) / (flow_top2 - flow_bottom2)

    # --- Efficiency axis for points ---
    eff_min_data = merged["η"].min()
    eff_max_data = merged["η"].max()
    eff_data_span = eff_max_data - eff_min_data
    pad_eff2 = eff_data_span * 0.10 if eff_data_span != 0 else abs(eff_max_data) * 0.10
    eff_bottom2 = eff_min_data - pad_eff2
    eff_top2 = eff_max_data + pad_eff2

    eff_span2 = eff_top2 - eff_bottom2
    cand_bottom2 = 0 - f_flow2 * eff_span2
    cand_top2 = cand_bottom2 + eff_span2

    # Ensure bottom ≤ actual min
    if cand_bottom2 > eff_bottom2:
        shift2 = cand_bottom2 - eff_bottom2
        cand_bottom2 -= shift2
        cand_top2    -= shift2

    # Ensure top ≥ actual max
    if cand_top2 < eff_top2:
        shift2 = eff_top2 - cand_top2
        cand_top2    += shift2
        cand_bottom2 += shift2

    ax2.set_ylim(cand_bottom2, cand_top2)

    # Format efficiency y‐axis in scientific notation with three significant digits
    ax2.yaxis.set_major_formatter(FormatStrFormatter('%.3e'))

    # Labels, title, grid, and legend
    ax1.set_xlabel("Frequency [Hz]")
    ax1.set_ylabel("Flow Rate", color="blue")
    ax1.tick_params(axis="y", labelcolor="blue")
    ax2.set_ylabel("Efficiency [η]", color="green")
    ax2.tick_params(axis="y", labelcolor="green")

    fig.suptitle(
        f"Flow Rate and Efficiency vs Frequency  |  points only  |  {SUPPLY_VOLTAGE} V",
        fontsize=11
    )

    handles1, labels1 = ax1.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(handles1 + handles2, labels1 + labels2, fontsize="x-small", loc="upper left")

    ax1.grid(True)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(os.path.join(out_dir, "flow_eff_points.png"), dpi=300)
    plt.close()







    
    
    

# ────────── All‑points & Relative plots ────────────────────────────────────
def all_points(tests, err):
    tag = 'all_points_with_err' if err else 'all_points_without_err'
    freq_dir, time_dir = subdir(tag,"frequency"), subdir(tag,"time")
    cmap, step = plt.get_cmap('tab10'), BIN_STEP_INDIVIDUAL
    for f, df in tests:
        tn,s_f,e_f,sw_s = extract_test_info(os.path.basename(f))
        # frequency scatter
        fig, ax = plt.subplots(figsize=(8,5)); handles=[]
        for i,g in df.groupby('Sweep'):
            c=cmap(i%10); ax.scatter(g['Frequency'], g['Flow'], s=POINT_SIZE, color=c)
            handles.append(plt.Line2D([0],[0],marker='o',color=c,ls='',label=f"Sweep {i+1}",markersize=6))
        tmp=df.copy(); tmp['bin']=_bin(tmp['Frequency'],step); gb=tmp.groupby('bin')['Flow'].mean()
        star=add_best_marker(ax, gb.idxmax(), gb.max()); handles.append(star)
        _legend_with_handles(ax, handles, fontsize='xx-small', ncol=3)
        ax.set_xlabel("Frequency (Hz)"); _label(ax)
        ax.set_title(f"{MEASUREMENT_NAME} vs Frequency  |  All raw points  |  Test {tn}  |  {s_f:.1f}→{e_f:.1f} Hz, {sw_s:.1f}s  |  {SUPPLY_VOLTAGE} V", fontsize=9)
        plt.tight_layout(); plt.savefig(os.path.join(freq_dir,f"{tn}.png"), dpi=300); plt.close(fig)
        # time scatter
        fig, ax = plt.subplots(figsize=(8,5)); handles=[]
        for i,g in df.groupby('Sweep'):
            c=cmap(i%10); ax.scatter(g['Elapsed_s'], g['Flow'], s=POINT_SIZE, color=c)
            handles.append(plt.Line2D([0],[0],marker='o',color=c,ls='',label=f"Sweep {i+1}",markersize=6))
        star=add_best_marker(ax, df.loc[df['Flow'].idxmax(),'Elapsed_s'], df['Flow'].max()); handles.append(star)
        _legend_with_handles(ax, handles, fontsize='xx-small', ncol=3)
        ax.set_xlabel("Elapsed Time (s)"); _label(ax)
        ax.set_title(f"{MEASUREMENT_NAME} vs Time  |  All raw points  |  Test {tn}  |  {s_f:.1f}→{e_f:.1f} Hz, {sw_s:.1f}s  |  {SUPPLY_VOLTAGE} V", fontsize=9)
        plt.tight_layout(); plt.savefig(os.path.join(time_dir,f"{tn}.png"), dpi=300); plt.close(fig)

def relative_plots(tests):
    rel_dir, cmap = subdir("relative","frequency"), plt.get_cmap('tab10')
    for f, df in tests:
        tn,s_f,e_f,sw_s = extract_test_info(os.path.basename(f))
        fig, ax = plt.subplots(figsize=(8,5)); handles=[]
        for i,g in df.groupby('Sweep'):
            c=cmap(i%10); rel=100*g['Flow']/g['Flow'].max()
            ax.plot(g['Frequency'], rel,'-o',ms=POINT_SIZE,color=c)
            handles.append(plt.Line2D([0],[0],marker='o',color=c,ls='',label=f"Sweep {i+1}",markersize=6))
        ax.set(xlabel="Frequency [Hz]", ylabel="Relative [%]", ylim=(0,105)); ax.grid(True)
        ax.set_title(f"{MEASUREMENT_NAME} vs Frequency  |  Relative 0‑100 %  |  Test {tn}  |  {s_f:.1f}→{e_f:.1f} Hz, {sw_s:.1f}s  |  {SUPPLY_VOLTAGE} V", fontsize=9)
        _legend_with_handles(ax, handles, fontsize='xx-small', ncol=3)
        plt.tight_layout(); plt.savefig(os.path.join(rel_dir,f"{tn}.png"), dpi=300); plt.close(fig)

# ────────── High‑flow / High‑stability cloud ───────────────────────────────
def high_quality_cloud(stat_df, title, outfile):
    """
    Highlight top-25% flow & top-10% stability bins.
    Picks the point with lowest std (error) after filtering.
    """
    stat_df = stat_df.dropna(subset=['mean', 'std'])
    if stat_df.empty:
        return

    # Filter top 25% by mean
    top25 = stat_df[stat_df['mean'] >= stat_df['mean'].quantile(0.75)]
    if top25.empty:
        return

    # From top25, select bottom 10% by std (i.e., most stable)
    blue = top25[top25['std'] <= top25['std'].quantile(0.10)]
    if blue.empty:
        blue = top25.copy()

    # Pick the most stable point (lowest std)
    best_idx = blue['std'].idxmin()
    best_val = blue.loc[best_idx, 'mean']
    best_err = blue.loc[best_idx, 'std']

    fig, ax = plt.subplots(figsize=(8, 5))

    # Plot all bins
    ax.errorbar(
        stat_df.index, stat_df['mean'], stat_df['std'],
        fmt='o', markersize=POINT_SIZE, color='lightgrey',
        ecolor='lightgrey', capsize=3, alpha=0.5,
        label='All bins', zorder=1
    )

    # Highlight selected "good" region
    ax.errorbar(
        blue.index, blue['mean'], blue['std'],
        fmt='o', markersize=POINT_SIZE, color='C0',
        capsize=4,
        label='Top-25% flow & top-10% stability',
        zorder=2
    )

    # Highlight best (most stable) point
    add_best_marker(ax, best_idx, best_val, best_err)  # z-order 5 inside

    ax.set_xlabel("Frequency [Hz]")
    ax.set_ylabel(YLABEL)
    ax.grid(True)
    ax.set_title(title)
    ax.legend()
    plt.tight_layout()
    plt.savefig(outfile, dpi=300)
    plt.close()



# ────────── MAIN ────────────────────────────────────────────────────────────
def main():
    
    if not os.path.isdir(parent_folder):
        print("bad parent")
        return

    # 1) Discover and load runs
    runs = [
        os.path.join(parent_folder, d)
        for d in os.listdir(parent_folder)
        if os.path.isdir(os.path.join(parent_folder, d)) and re.search(r'\d+Hz_\d+Hz_\d+s', d)
    ]
    if not runs:
        print("no runs")
        return
    
    tests = []     # will hold (path, full_df)
    frames = []    # will hold just ['Frequency','Flow'] for averaging
    for r in runs:
        try:
            tn, s_f, e_f, sw_s = extract_test_info(os.path.basename(r))
            df = load_and_process_csv(r, s_f, e_f, sw_s)
            tests.append((r, df))
            frames.append(df[['Frequency', 'Flow']])
        except Exception as exc:
            print("Skip", r, "→", exc)
    if not tests:
        print("no valid")
        return

    # 2) Create two concatenations:
    #    a) Full raw DataFrame for freq/time vs flow/pressure correlations
    concat_raw = pd.concat([df for _, df in tests], ignore_index=True)
    
    #    b) Slim DataFrame for binned averages
    concat = pd.concat(frames, ignore_index=True)
    
    # 3) Compute binned-average stats
    bins = _bin(concat['Frequency'], BIN_STEP_AVERAGE)
    stats_df = concat.groupby(bins)['Flow'].agg(['mean', 'std'])
    freq_avg = (
        stats_df
        .reset_index()
        .rename(columns={'index': 'Frequency'})
    )
    
    #Print the data summary:
    print_data_summary(tests)
  
    
    stats_df=concat.groupby(bins)['Flow'].agg(['mean','std'])
    indiv_plot_freq(tests,subdir("freq_individual_with_err"),True)
    indiv_plot_time(tests,subdir("time_individual_with_err"),True)
    indiv_plot_freq(tests,subdir("freq_individual_no_err"),False)
    indiv_plot_time(tests,subdir("time_individual_no_err"),False)
    
    # Generate point‐only and error‐bar plots for both Frequency and Time

    # Frequency (Hz) plots
    comb_series(tests, subdir("freq_combined_series"),     True,  error_bars=False)
    comb_series(tests, subdir("freq_combined_series_err"), True,  error_bars=True)
    
    # Elapsed Time (s) plots
    comb_series(tests, subdir("time_combined_series"),     False, error_bars=False)
    comb_series(tests, subdir("time_combined_series_err"), False, error_bars=True)

    avg_plot(tests,subdir("freq_average_with_err"),True,True)
    avg_plot(tests,subdir("time_average_with_err"),False,True)
    avg_plot(tests,subdir("freq_average_no_err"),True,False)
    avg_plot(tests,subdir("time_average_no_err"),False,False)
    all_points(tests,True); all_points(tests,False); relative_plots(tests)

    cloud_dir=subdir("high_flow_high_stability")
    for run_dir, df in tests:
        tn,*_ = extract_test_info(os.path.basename(run_dir))
        stat = df.groupby(_bin(df['Frequency'],BIN_STEP_AVERAGE))['Flow'].agg(['mean','std'])
        high_quality_cloud(stat, f"{MEASUREMENT_NAME}: Top‑25 % → Top‑10 % – Test {tn}  |  {SUPPLY_VOLTAGE} V",
                           os.path.join(cloud_dir,f"Test_{tn}.png"))
    high_quality_cloud(stats_df, f"{MEASUREMENT_NAME}: Top‑25 % → Top‑10 % – ALL tests  |  {SUPPLY_VOLTAGE} V",
                       os.path.join(cloud_dir,"overall.png"))

    resid_dir=subdir("residuals")
    for run_dir,df in tests:
        tn,*_=extract_test_info(os.path.basename(run_dir))
        residual_plots(df['Flow'],os.path.join(resid_dir,f"{tn}.png"),
                       f"{MEASUREMENT_NAME}: Residuals – Test {tn}  |  {SUPPLY_VOLTAGE} V")
    residual_plots(concat['Flow'],os.path.join(resid_dir,"average.png"),
                   f"{MEASUREMENT_NAME}: Residuals – Average  |  {SUPPLY_VOLTAGE} V")
    
    
    #We can also add the efficiency plots and calculations
    concat=pd.concat(frames)
    bins=_bin(concat['Frequency'],BIN_STEP_AVERAGE)
    freq_avg=concat.groupby(bins)['Flow'].agg(['mean','std']).reset_index().rename(columns={'bin':'Frequency'})

    eff_df = compute_efficiency_with_error(freq_avg)
    efficiency_dir = subdir("efficiency_plots")

    plot_efficiency(eff_df, efficiency_dir, err=True, POINT_SIZE=POINT_SIZE)
    plot_efficiency(eff_df, efficiency_dir, err=False, POINT_SIZE=POINT_SIZE)

    plot_efficiency_with_avg(eff_df, freq_avg, efficiency_dir, err=True, POINT_SIZE=POINT_SIZE)
    plot_efficiency_with_avg(eff_df, freq_avg, efficiency_dir, err=False, POINT_SIZE=POINT_SIZE)
    

    # Plot std vs flow for average
    std_plot_dir = subdir("std_vs_flow")
    avg_out_path = os.path.join(std_plot_dir, "avg_std_vs_flow.png")
    
    # Concatenate all data, bin, and average like for individuals
    all_df = pd.concat([d[['Frequency', 'Flow']] for _, d in tests])
    all_df['bin'] = _bin(all_df['Frequency'], BIN_STEP_AVERAGE)
    ag_all = all_df.groupby('bin')['Flow'].agg(['mean', 'std']).reset_index()
    plot_std_vs_flow(ag_all, avg_out_path, "Average of all tests", point_size=POINT_SIZE)

    
    # Plot std vs flow for each individual test
    for folder, df in tests:
        tn, *_ = extract_test_info(os.path.basename(folder))
        step = BIN_STEP_INDIVIDUAL
        df['bin'] = _bin(df['Frequency'], step)
        ag = df.groupby('bin')['Flow'].agg(['mean', 'std']).reset_index()
        out_path = os.path.join(std_plot_dir, f"test_{tn}_std_vs_flow.png")
        plot_std_vs_flow(ag, out_path, f"Test {tn}", point_size=POINT_SIZE)
        
    
    #Exploratory plots.
    eda_dir = subdir("exploratory_plots")
    UNIT = ("µL/min", "mbar")[PRESSURE_TEST]

    # 3a) Raw-data visuals for each test
    plot_individual_histograms(tests, eda_dir, UNIT)
    plot_individual_boxplots(   tests, eda_dir)

    # 3b) visuals across all tests
    plot_average_histogram_from_raw(tests, eda_dir, UNIT, bin_step_hist_box_average)
    
    plot_average_boxplot_from_raw(tests, eda_dir)


    # 3c) Correlation matrices
    plot_raw_correlation(     concat_raw, eda_dir)
    plot_aggregate_correlation(freq_avg, eda_dir)
    
    #Combined boxplot figure
    plot_combined_boxplots(tests, eda_dir, width=0.6)

    print("Finished successfully.")


if __name__=="__main__":
    main()
