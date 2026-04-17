#!/usr/bin/env python3

import argparse
import json
import os
import re
import subprocess
import shutil
import time

import numpy as np

try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    HAS_PLOTLY = True
except ImportError:
    HAS_PLOTLY = False

def run_cmd(cmd, cwd, log_path):
    print(f"Running: {' '.join(cmd)}")
    with open(log_path, "w") as f:
        p = subprocess.Popen(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        for line in p.stdout:
            f.write(line)
            # print(line, end="") # Optional: print to console
        return p.wait()

def main():
    parser = argparse.ArgumentParser(description="General Free Energy Sweep")
    parser.add_argument("--config", default="default_free_energy_config.json", help="JSON config file")
    args = parser.parse_args()

    with open(args.config, "r") as f:
        cfg = json.load(f)

    output_root = cfg.get("output_root", "results_sweep")
    common = cfg.get("common", {})
    param_sets = cfg.get("parameter_sets", [{}])

    os.makedirs(output_root, exist_ok=True)
    
    # Copy config to output for record
    shutil.copy2(args.config, os.path.join(output_root, "config.json"))

    summary_file = os.path.join(output_root, "summary.txt")
    with open(summary_file, "w") as sf:
        sf.write("Run Name | Status | Wall Time (s)\n")
        sf.write("-" * 40 + "\n")

    for pset in param_sets:
        # Merge common and specific parameters
        params = common.copy()
        params.update(pset)
        
        run_name = params.get("name", "unnamed_run")
        run_dir = os.path.join(output_root, run_name)
        os.makedirs(run_dir, exist_ok=True)
        
        print(f"\n>>> Starting run: {run_name}")
        
        cmd = ["python3", "run_ES.py"]
        # Add arguments from params
        for key, val in params.items():
            if key == "name": continue
            if isinstance(val, bool):
                if val: cmd.append(f"--{key}")
            elif key == "temperature":
                cmd.extend(["-T", str(val)])
            else:
                cmd.extend([f"--{key}", str(val)])

        log_path = os.path.join(run_dir, "run.log")
        t0 = time.time()
        rc = run_cmd(cmd, cwd=os.getcwd(), log_path=log_path)
        dt = time.time() - t0
        
        status = "OK" if rc == 0 else "FAILED"
        print(f"<<< Finished {run_name} with status {status} in {dt:.2f}s")
        
        with open(summary_file, "a") as sf:
            sf.write(f"{run_name} | {status} | {dt:.2f}\n")

        # Find the generated .dat file
        xyz_base = os.path.splitext(os.path.basename(params.get("xyz_name", "")))[0]
        dat_file = None
        for f in os.listdir("."):
            if f.endswith(".dat") and (xyz_base in f or "free_energy" in f):
                dat_file = f
                break

        if dat_file and rc == 0:
            print(f"Generating plot for {dat_file}...")
            plot_prefix = os.path.join(run_dir, f"{run_name}_plot")
            plot_cmd = [
                "python3", "plot_F_interactive.py",
                "--input", dat_file,
                "--output", plot_prefix,
                "--T", str(params.get("temperature", 300.0))
            ]
            # If we have N atoms information, we could pass it here as well
            run_cmd(plot_cmd, cwd=os.getcwd(), log_path=os.path.join(run_dir, "plot.log"))

        # Move all relevant files to run_dir
        for f in os.listdir("."):
            if f.endswith(".dat") or f.endswith(".html"):
                if xyz_base in f or "free_energy" in f or run_name in f:
                    target = os.path.join(run_dir, f)
                    if os.path.exists(target): os.remove(target)
                    shutil.move(f, target)

    # After all runs: generate combined comparison plot
    plot_sweep_comparison(output_root)

def _parse_header_cols(path):
    """Return dict {name: col_index} from the first '#'-comment line containing 'lambda'."""
    col_map = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line.startswith('#'):
                break
            if 'lambda' in line:
                parts = line.lstrip('#').split()
                for i, p in enumerate(parts):
                    col_map[p] = i
    return col_map


def _load_fe_dat(path):
    """Load a free-energy .dat file; returns (distance, TI_F, TI_err) arrays."""
    col_map = _parse_header_cols(path)
    data = np.loadtxt(path, comments='#')
    if data.ndim == 1:
        data = data[None, :]
    ncols = data.shape[1]

    x_col  = min(col_map.get('distance', col_map.get('lambda', 0)),  ncols - 1)
    f_col  = min(col_map.get('TI_F',    3),  ncols - 1)
    e_col  = min(col_map.get('TI_err',  4),  ncols - 1)
    de_col = min(col_map.get('TI_dE/dlambda', 1), ncols - 1)

    x   = data[:, x_col]
    ti  = data[:, f_col]
    err = data[:, e_col]
    de  = data[:, de_col]
    lam = data[:, col_map.get('lambda', 0)]

    # zero-shift TI_F
    valid = ~np.isnan(ti)
    if np.any(valid):
        ti = ti - ti[valid][0]

    return x, ti, err, de, lam


def plot_sweep_comparison(output_root, output_html=None):
    """Scan output_root for per-run subdirectories, load their .dat files
    and write a combined Plotly comparison to output_html."""
    if not HAS_PLOTLY:
        print("plotly not installed – skipping comparison plot.")
        return

    if output_html is None:
        output_html = os.path.join(output_root, "comparison_temperatures.html")

    # Collect (temperature, path) pairs
    runs = []
    for entry in sorted(os.listdir(output_root)):
        run_dir = os.path.join(output_root, entry)
        if not os.path.isdir(run_dir):
            continue
        dat_files = [f for f in os.listdir(run_dir) if f.endswith('.dat')]
        if not dat_files:
            continue
        # Extract temperature from name like temp_300K or nSys_200
        m = re.search(r'(\d+(?:\.\d+)?)[Kk]', entry)
        T_val = float(m.group(1)) if m else None
        runs.append((T_val, entry, os.path.join(run_dir, dat_files[0])))

    if not runs:
        print(f"No .dat files found under {output_root} – nothing to plot.")
        return

    # Sort by temperature (None last)
    runs.sort(key=lambda r: (r[0] is None, r[0] or 0))

    # Build color scale (blue → red across temperature range)
    temps = [r[0] for r in runs if r[0] is not None]
    T_min, T_max = (min(temps), max(temps)) if len(temps) > 1 else (temps[0], temps[0]+1)
    def temp_color(T):
        if T is None:
            return 'grey'
        t = (T - T_min) / (T_max - T_min) if T_max > T_min else 0.5
        r = int(255 * t)
        b = int(255 * (1 - t))
        return f'rgb({r},0,{b})'

    fig = make_subplots(
        rows=2, cols=1,
        subplot_titles=('Cumulative Free Energy ΔF vs distance', 'dE/dλ vs λ'),
        vertical_spacing=0.12,
    )

    for T_val, name, dat_path in runs:
        try:
            x, ti, err, de, lam = _load_fe_dat(dat_path)
        except Exception as exc:
            print(f"  Skipping {name}: {exc}")
            continue

        label  = f"{T_val:.0f} K" if T_val is not None else name
        color  = temp_color(T_val)
        errarr = err if not np.all(np.isnan(err)) else None

        # Row 1: ΔF vs distance
        fig.add_trace(go.Scatter(
            x=x, y=ti, name=label,
            mode='lines',
            line=dict(color=color, width=2),
            error_y=dict(type='data', array=errarr, visible=(errarr is not None),
                         color=color.replace('rgb', 'rgba').replace(')', ',0.3)')),
            legendgroup=label,
            hovertemplate=f"<b>{label}</b><br>d=%{{x:.2f}} Å<br>ΔF=%{{y:.4f}} eV<extra></extra>",
        ), row=1, col=1)

        # Row 2: dE/dlambda vs lambda
        fig.add_trace(go.Scatter(
            x=lam, y=de, name=label,
            mode='lines',
            line=dict(color=color, width=1.5),
            legendgroup=label,
            showlegend=False,
            hovertemplate=f"<b>{label}</b><br>λ=%{{x:.4f}}<br>dE/dλ=%{{y:.4f}} eV<extra></extra>",
        ), row=2, col=1)

    fig.update_xaxes(title_text="End-to-end distance [Å]", row=1, col=1)
    fig.update_yaxes(title_text="ΔF [eV]", row=1, col=1)
    fig.update_xaxes(title_text="λ", row=2, col=1)
    fig.update_yaxes(title_text="dE/dλ [eV]", row=2, col=1)
    fig.update_layout(
        height=900,
        title=f"Free Energy Comparison – {os.path.basename(output_root)}",
        template='plotly_white',
        hovermode='x unified',
        legend=dict(title="Temperature", orientation='v'),
    )

    fig.write_html(output_html)
    print(f"Comparison plot saved to: {output_html}")


if __name__ == "__main__":
    main()
