#!/usr/bin/env python3

import os
import glob
import re
import argparse
import numpy as np
import matplotlib.pyplot as plt

try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    HAS_PLOTLY = True
except ImportError:
    HAS_PLOTLY = False

def detect_fe_column(filename: str) -> int:
    """
    Detects the column index for free energy ('TI_F' or 'cumulative_FE').
    Returns the column index or -1 if not found.
    """
    with open(filename, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line.startswith("#"):
                # Stop if we are past the header comments
                break
            # Look for a header line with column names
            if "lambda" in line and ("TI_F" in line or "cumulative_FE" in line):
                parts = line.replace("#", "").split()
                for i, p in enumerate(parts):
                    if p in ("TI_F", "cumulative_FE"):
                        return i
    return -1

def load_sweep_data(path: str):
    """
    Load data for interactive plotting.
    Parses header columns dynamically and extracts distance, lambda, TI_F, TI_err, TI_dE/dlambda,
    and any hb_*_fraction columns.
    """
    col_map = {}
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line.startswith("#"):
                break
            if "lambda" in line:
                parts = line.lstrip("#").split()
                for i, p in enumerate(parts):
                    col_map[p] = i

    d = np.loadtxt(path, comments="#")
    if d.ndim == 1:
        d = d[None, :]

    ncols = d.shape[1]
    
    # Identify columns
    x_col = min(col_map.get("distance", col_map.get("lambda", 0)), ncols - 1)
    f_col = min(col_map.get("TI_F", 3), ncols - 1)
    e_col = min(col_map.get("TI_err", 4), ncols - 1)
    de_col = min(col_map.get("TI_dE/dlambda", 1), ncols - 1)
    lam_col = min(col_map.get("lambda", 0), ncols - 1)

    x = d[:, x_col]
    ti = d[:, f_col]
    err = d[:, e_col]
    de = d[:, de_col]
    lam = d[:, lam_col]

    # Shift TI to start at 0
    valid = ~np.isnan(ti)
    if np.any(valid):
        ti = ti - ti[valid][0]

    # Detect hb_fraction columns
    hb_frac = None
    hb_cols = sorted([col for col in col_map if col.startswith("hb_") and col.endswith("_fraction")])
    if hb_cols:
        # Use the first hydrogen bond fraction column
        hb_frac = d[:, col_map[hb_cols[0]]]

    return x, ti, err, de, lam, hb_frac

def main():
    parser = argparse.ArgumentParser(description="Plot difference of max(F) - min(F) vs Temperature and generate interactive HTML sweep plots")
    parser.add_argument("--results_dir", default=".", help="Directory containing the temperature subdirectories (e.g., 'results_sweep')")
    args = parser.parse_args()

    results_dir = args.results_dir
    
    # Recursively find all .dat files in the directory
    files = glob.glob(os.path.join(results_dir, "**", "*.dat"), recursive=True)
    
    if not files:
        print(f"No .dat files found recursively in '{results_dir}'.")
        return

    # Filter to only keep files matching the most common basename (ignoring spurious files)
    from collections import Counter
    basenames = [os.path.basename(f) for f in files]
    most_common_basename, count = Counter(basenames).most_common(1)[0]
    print(f"Filtering to the most common output file name: '{most_common_basename}' (found {count} times)")
    files = [f for f in files if os.path.basename(f) == most_common_basename]

    runs = []
    print(f"Found {len(files)} matching '.dat' files. Processing...")

    for f in files:
        # Extract temperature from the directory name (e.g., matching "temp_300K")
        m = re.search(r'temp_(\d+(?:\.\d+)?)K', f)
        if not m:
            print(f"Skipping {f}: could not determine temperature from path.")
            continue
        
        T = float(m.group(1))
        runs.append((T, f))

    if not runs:
        print(f"No valid temperature runs identified from the found files in '{results_dir}'.")
        return

    # Sort runs by temperature
    runs.sort(key=lambda r: r[0])

    temperatures = []
    delta_Fs = []
    valid_sweep_data = []

    for T, f in runs:
        try:
            fe_col_idx = detect_fe_column(f)
            if fe_col_idx == -1:
                print(f"Warning: Could not detect a Free Energy column ('TI_F' or 'cumulative_FE') in the header of '{f}'. Skipping.")
                continue

            data = np.loadtxt(f, comments="#")
            
            if data.size == 0:
                continue
            if data.ndim == 1:
                data = data.reshape(1, -1)
            
            if fe_col_idx >= data.shape[1]:
                print(f"Skipping {f}: detected FE column {fe_col_idx} is out of bounds for data shape {data.shape}.")
                continue
                
            # Use the detected column for Free Energy
            F = data[:, fe_col_idx]
            delta_F = np.max(F) - np.min(F)
            
            temperatures.append(T)
            delta_Fs.append(delta_F)
            
            # Load full sweep details for Plotly
            x, ti, err, de, lam, hb_frac = load_sweep_data(f)
            valid_sweep_data.append((T, x, ti, err, de, lam, hb_frac))
            
        except Exception as e:
            print(f"Skipping {f} due to error: {e}")
            
    if not temperatures:
        print(f"No valid data was processed from the found files in '{results_dir}'.")
        return
        
    temperatures = np.array(temperatures)
    delta_Fs = np.array(delta_Fs)
    
    print(f"Plotting data for {len(temperatures)} temperatures.")

    # 1. Static Plot (Matplotlib)
    plt.figure(figsize=(8, 5))
    plt.plot(temperatures, delta_Fs, marker='o', linestyle='-', linewidth=2, color='#2196F3')
    plt.xlabel('Temperature [K]', fontsize=12)
    plt.ylabel(r'$\Delta F = \max(F) - \min(F)$ [eV]', fontsize=12)
    plt.title('Maximum Free Energy Difference vs Temperature', fontsize=14)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    
    out_path = os.path.join(results_dir, "delta_F_vs_T.png")
    plt.savefig(out_path, dpi=150)
    print(f"Static plot saved to {out_path}")
    plt.close()

    # 2. Interactive Plot (Plotly HTML)
    if HAS_PLOTLY and valid_sweep_data:
        # Temperature limits for color mapping
        T_min, T_max = (temperatures.min(), temperatures.max()) if len(temperatures) > 1 else ((temperatures[0], temperatures[0] + 1) if len(temperatures) > 0 else (0, 1))

        def temp_color(T):
            if T_max > T_min:
                t = (T - T_min) / (T_max - T_min)
            else:
                t = 0.5
            r = int(255 * t)
            b = int(255 * (1 - t))
            return f"rgb({r},0,{b})"

        fig = make_subplots(
            rows=3, cols=1,
            subplot_titles=(
                "Cumulative Free Energy ΔF vs distance",
                "dE/dλ vs λ",
                "Hydrogen Bond Fraction vs distance"
            ),
            vertical_spacing=0.08
        )

        for T, x, ti, err, de, lam, hb_frac in valid_sweep_data:
            label = f"{T:.1f} K" if T % 1 != 0 else f"{T:.0f} K"
            color = temp_color(T)
            errarr = err if not np.all(np.isnan(err)) else None
            err_color = None if errarr is None else color.replace("rgb", "rgba").replace(")", ",0.3)")

            # Trace 1: Free Energy
            fig.add_trace(
                go.Scatter(
                    x=x, y=ti,
                    name=label,
                    legendgroup=label,
                    mode="lines",
                    line=dict(color=color, width=2),
                    error_y=dict(type="data", array=errarr, visible=(errarr is not None), color=err_color),
                    hovertemplate=f"<b>{label} Free Energy</b><br>dist = %{{x:.3f}} Å<br>ΔF = %{{y:.4f}} eV<extra></extra>"
                ),
                row=1, col=1
            )

            # Trace 2: Derivative
            fig.add_trace(
                go.Scatter(
                    x=lam, y=de,
                    name=label,
                    legendgroup=label,
                    mode="lines",
                    line=dict(color=color, width=1.5),
                    showlegend=False,
                    hovertemplate=f"<b>{label} dE/dλ</b><br>λ = %{{x:.3f}}<br>dE/dλ = %{{y:.4f}} eV<extra></extra>"
                ),
                row=2, col=1
            )

            # Trace 3: HBond Fraction
            hb_y = hb_frac if hb_frac is not None else np.full_like(x, np.nan)
            fig.add_trace(
                go.Scatter(
                    x=x, y=hb_y,
                    name=label,
                    legendgroup=label,
                    mode="lines+markers" if hb_frac is not None else "lines",
                    marker=dict(size=4),
                    line=dict(color=color, width=1.5),
                    showlegend=False,
                    hovertemplate=f"<b>{label} HB Fraction</b><br>dist = %{{x:.3f}} Å<br>Fraction = %{{y:.4f}}<extra></extra>"
                ),
                row=3, col=1
            )

        fig.update_xaxes(title_text="End-to-end distance [Å]", row=1, col=1)
        fig.update_yaxes(title_text="ΔF [eV]", row=1, col=1)
        fig.update_xaxes(title_text="λ", row=2, col=1)
        fig.update_yaxes(title_text="dE/dλ [eV]", row=2, col=1)
        fig.update_xaxes(title_text="End-to-end distance [Å]", row=3, col=1)
        fig.update_yaxes(title_text="HB Fraction", row=3, col=1)

        fig.update_layout(
            title=f"Temperature Sweep Analysis – {os.path.basename(os.path.abspath(results_dir))}",
            height=1000,
            template="plotly_white",
            hovermode="x unified",
            legend=dict(traceorder="normal")
        )

        out_html = os.path.join(results_dir, "temperature_sweep_interactive.html")
        fig.write_html(out_html)
        print(f"Interactive plot saved to {out_html}")
    elif not HAS_PLOTLY:
        print("plotly is not installed – skipping interactive HTML sweep generation.")

    # Display the static plot if in a headless environment
    if 'DISPLAY' in os.environ or os.name == 'nt':
        try:
            plt.show()
        except Exception as e:
            print(f"Could not display plot interactively: {e}")
    else:
        print("No display found. Static plot was saved to file.")

if __name__ == "__main__":
    main()
