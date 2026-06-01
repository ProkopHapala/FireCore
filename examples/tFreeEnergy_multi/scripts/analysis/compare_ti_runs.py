#!/usr/bin/env python3
"""
compare_ti_runs.py  –  Run and compare CPU / GPU-debug / multi-GPU TI for entropic spring.

Runs all three implementations, loads their results, saves one comparison data file,
and writes one interactive HTML plot (Plotly): curves on the left, run parameters on the right.

Usage (from examples/tFreeEnergy_multi/):
    python3 compare_ti_runs.py [-T 300] [--skip_cpu] [--skip_gpu] [--skip_multi]
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
except ImportError:
    print("plotly is required:  pip install plotly")
    sys.exit(1)

HERE    = Path(__file__).parent            # examples/tFreeEnergy_multi/
FE_DIR  = HERE.parent / "tFreeEnergy"     # examples/tFreeEnergy/


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def gaussian_chain_reference(lam, lam0, T, n_segments, bond_length=1.198):
    """Analytical free energy of a Gaussian chain:  F(λ) – F(λ0)  [eV]."""
    k_b   = 8.617333262145e-5
    k_eff = 3.0 * k_b * T / (n_segments * bond_length**2)
    return 0.5 * k_eff * (lam**2 - lam0**2)


def infer_n_segments(xyz_path: Path) -> int:
    m = re.search(r"entropic_spring_(\d+)", xyz_path.name)
    if not m:
        raise ValueError(f"Cannot infer chain size from: {xyz_path.name}")
    return int(m.group(1)) - 1


def load_cpu_gpu(path: Path):
    """Load run_ES.py / run_ES_gpu_debug.py output:  lambda  TI_F  TI_sigma  Ref_F  TI-Ref"""
    d = np.loadtxt(path, comments="#")
    return d[:, 0], d[:, 1], d[:, 2]   # lam, TI_F, sigma


def load_multi(path: Path):
    """Load tFreeEnergy_multi/run_ES.py output.

    The file header is parsed to find the column indices for TI_F, TI_err, and
    distance.  The physical end-to-end *distance* (Å) is used as the x-axis so
    that it aligns with the CPU/GPU-debug curves which also use distance [Å].
    """
    # --- parse header for column mapping ---
    col_map = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line.startswith('#'):
                break
            if 'lambda' in line:
                parts = line.replace('#', '').split()
                for i, p in enumerate(parts):
                    if p == 'lambda':             col_map['lambda']   = i
                    elif p == 'TI_F':             col_map['TI_F']     = i
                    elif p == 'TI_err':           col_map['TI_err']   = i
                    elif p == 'distance':         col_map['distance'] = i

    d = np.loadtxt(path, comments="#")
    if d.ndim == 1:
        d = d[None, :]
    ncols = d.shape[1]

    # Prefer physical distance as x-axis; fall back to lambda param
    x_col = col_map.get('distance', col_map.get('lambda', 0))
    f_col = col_map.get('TI_F',  3)
    e_col = col_map.get('TI_err', 4)
    # clamp to actual column count
    x_col = min(x_col, ncols - 1)
    f_col = min(f_col, ncols - 1)
    e_col = min(e_col, ncols - 1)

    x   = d[:, x_col]
    ti  = d[:, f_col]
    err = d[:, e_col]

    # Ensure TI_F starts at 0 (matches CPU/GPU-debug convention)
    valid = ~np.isnan(ti)
    if np.any(valid):
        ti = ti - ti[valid][0]

    return x, ti, err


def run_script(script: Path, extra_args: list, cwd: Path) -> bool:
    """Run script as subprocess. Returns True on success, False on failure."""
    cmd = [sys.executable, str(script)] + extra_args
    print(f"\n>>> {' '.join(cmd)}")
    r = subprocess.run(cmd, cwd=cwd)
    if r.returncode != 0:
        print(f"\n[WARNING] {script.name} failed (exit {r.returncode}) – skipping this method.\n")
        return False
    return True



# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Run and compare CPU / GPU-debug / multi-GPU TI for entropic spring."
    )
    # Common parameters
    ap.add_argument("--N",        type=int,   default=20,    help="Entropic spring size (number of atoms)")
    ap.add_argument("--lamda1",   type=float, default=1.0,   help="Start of lambda range [Å]")
    ap.add_argument("--lamda2",   type=float, default=11.0,  help="End of lambda range [Å]")
    ap.add_argument("-T", "--temperature", type=float, default=300.0, help="Temperature [K]")
    ap.add_argument("--dt",       type=float, default=0.05,  help="MD time step [fs]")
    ap.add_argument("--t_damp",   type=float, default=100.0, help="Langevin damping time [steps]")
    # Per-method step counts
    ap.add_argument("--cpu_nbStep",    type=int, default=50)
    ap.add_argument("--cpu_nMD",       type=int, default=50000)
    ap.add_argument("--cpu_nEQ",       type=int, default=5000)
    ap.add_argument("--gpu_nSys",      type=int, default=1)
    ap.add_argument("--gpu_nbStep",    type=int, default=25)
    ap.add_argument("--gpu_nMD",       type=int, default=10000)
    ap.add_argument("--gpu_nEQ",       type=int, default=5000)
    ap.add_argument("--multi_nSys",    type=int, default=10)
    ap.add_argument("--multi_nLambda", type=int, default=25)
    ap.add_argument("--multi_nMD",     type=int, default=10000)
    ap.add_argument("--multi_nEQ",     type=int, default=5000)
    # Skip flags
    ap.add_argument("--skip_cpu",   action="store_true", help="Skip CPU run, use existing results")
    ap.add_argument("--skip_gpu",   action="store_true", help="Skip GPU-debug run, use existing results")
    ap.add_argument("--skip_multi", action="store_true", help="Skip multi-GPU run, use existing results")
    args = ap.parse_args()

    T   = args.temperature
    N   = args.N
    xyz = f"../tMMFF/data/entropic_spring_{N}.xyz"
    cons = f"../tMMFF/data/entropic_spring_{N}.cons"
    n_segments = infer_n_segments(Path(FE_DIR.parent / "tMMFF" / "data" / f"entropic_spring_{N}.xyz"))

    # Expected output paths (written by each script)
    path_cpu   = FE_DIR / "results" / f"entropic_spring_{N}_TI_cpu_T{T:.0f}K.dat"
    path_gpu   = FE_DIR / "results" / f"entropic_spring_{N}_TI_gpu_debug_T{T:.0f}K.dat"
    path_multi = HERE / f"entropic_spring_{N}_free_energy.dat"

    common = ["--lamda1", str(args.lamda1), "--lamda2", str(args.lamda2),
              "-T", str(T), "--dt", str(args.dt), "--t_damp", str(args.t_damp)]

    # --- Run simulations ---
    ok_cpu = args.skip_cpu or run_script(FE_DIR / "run_ES.py", common + [
        "--xyz_name", xyz, "--constr_name", cons, "--constraints", "constraints_ES.txt",
        "--nbStep", str(args.cpu_nbStep), "--nMDsteps", str(args.cpu_nMD), "--nEQsteps", str(args.cpu_nEQ),
    ], cwd=FE_DIR)

    ok_gpu = args.skip_gpu or run_script(FE_DIR / "run_ES_gpu_debug.py", common + [
        "--nSys", str(args.gpu_nSys),
        "--xyz_name", xyz, "--constr_name", cons, "--constraints", "constraints_ES.txt",
        "--nbStep", str(args.gpu_nbStep), "--nMDsteps", str(args.gpu_nMD), "--nEQsteps", str(args.gpu_nEQ),
    ], cwd=FE_DIR)

    ok_mul = args.skip_multi or run_script(HERE / "run_ES.py", [
        "--nSys", str(args.multi_nSys),
        "--xyz_name", str(FE_DIR.parent / "tMMFF" / "data" / f"entropic_spring_{N}.xyz"),
        "--nLambda", str(args.multi_nLambda), "--nMDsteps", str(args.multi_nMD), "--nEQsteps", str(args.multi_nEQ),
        "--constraints", str(HERE / "constraints_ES.txt"),
        "-T", str(T), "--dt", str(args.dt), "--t_damp", str(args.t_damp),
        "--soft_dist",   # use distance constraint (matches run_ES.sh behaviour)
    ], cwd=HERE)

    if not (ok_cpu or ok_gpu or ok_mul):
        print("ERROR: all three runs failed – nothing to compare.")
        sys.exit(1)

    # --- Load available results ---
    curves = {}   # name → (lam, ti, sigma)
    if ok_cpu   and path_cpu.exists():   curves["CPU TI"]       = load_cpu_gpu(path_cpu)
    if ok_gpu   and path_gpu.exists():   curves["GPU-debug TI"] = load_cpu_gpu(path_gpu)
    if ok_mul   and path_multi.exists(): curves["multi-GPU TI"] = load_multi(path_multi)

    # Use the densest lambda grid available as the reference grid
    ref_lam = max(curves.values(), key=lambda c: len(c[0]))[0]
    ref = gaussian_chain_reference(ref_lam, args.lamda1, T, n_segments)

    # RMS vs reference for each available curve
    rms = {}
    for name, (lam, ti, _) in curves.items():
        ti_i = np.interp(ref_lam, lam, ti)
        rms[name] = float(np.sqrt(np.mean((ti_i - ref) ** 2)))

    # --- Save comparison data ---
    out_dir = HERE / "results_compare"
    out_dir.mkdir(exist_ok=True)
    out_stem = f"entropic_spring_{N}_TI_compare_T{T:.0f}K"
    out_dat  = out_dir / f"{out_stem}.dat"

    col_names = ["lambda[A]"] + [n.replace(" ", "_") + "[eV]" for n in curves] + ["Ref[eV]"]
    cols = [ref_lam]
    for lam, ti, _ in curves.values():
        cols.append(np.interp(ref_lam, lam, ti))
    cols.append(ref)
    header = "# " + "  ".join(col_names)
    np.savetxt(out_dat, np.column_stack(cols), header=header, comments="", fmt="%.8f")

    # --- Interactive Plotly HTML ---
    out_html = out_dir / f"{out_stem}.html"
    fig = make_subplots(
        rows=1, cols=2, column_widths=[0.70, 0.30],
        specs=[[{"type": "xy"}, {"type": "table"}]],
        subplot_titles=["ΔF(λ) — CPU / GPU-debug / multi-GPU", "Run parameters"],
    )

    # colour per method (keyed by display name)
    CLR = {
        "CPU TI":       "#2196F3",
        "GPU-debug TI": "#4CAF50",
        "multi-GPU TI": "#FF9800",
        "ref":          "#9C27B0",
    }

    for name, (lam, ti, sigma) in curves.items():
        color = CLR.get(name, "#888888")
        fig.add_trace(go.Scatter(
            x=lam, y=ti, name=name, mode="lines+markers",
            marker=dict(size=5), line=dict(color=color, width=2),
            error_y=dict(type="data", array=sigma, visible=True, color=color, thickness=1),
            hovertemplate=f"<b>{name}</b><br>λ = %{{x:.3f}} Å<br>ΔF = %{{y:.4f}} eV<extra></extra>",
        ), row=1, col=1)

    fig.add_trace(go.Scatter(
        x=ref_lam, y=ref, name="Gaussian chain (ref)", mode="lines",
        line=dict(color=CLR["ref"], dash="dash", width=2),
        hovertemplate="<b>Reference</b><br>λ = %{x:.3f} Å<br>ΔF = %{y:.4f} eV<extra></extra>",
    ), row=1, col=1)

    fig.update_xaxes(title_text="λ [Å]",  row=1, col=1)
    fig.update_yaxes(title_text="ΔF [eV]", row=1, col=1)

    # Parameters table – common block + one block per method
    labels = ["System", "N atoms", "T [K]", "λ range [Å]", "dt [fs]", "t_damp"]
    values = [f"entropic_spring_{N}", str(N), f"{T:.0f}", f"{args.lamda1}–{args.lamda2}",
              str(args.dt), str(args.t_damp)]

    method_params = {
        "CPU TI":       ("── CPU ──",       [("nbStep", args.cpu_nbStep), ("nMDsteps", args.cpu_nMD),   ("nEQsteps", args.cpu_nEQ)]),
        "GPU-debug TI": ("── GPU-debug ──", [("nSys", args.gpu_nSys), ("nbStep", args.gpu_nbStep), ("nMDsteps", args.gpu_nMD), ("nEQsteps", args.gpu_nEQ)]),
        "multi-GPU TI": ("── multi-GPU ──", [("nSys", args.multi_nSys), ("nLambda", args.multi_nLambda), ("nMDsteps", args.multi_nMD), ("nEQsteps", args.multi_nEQ)]),
    }
    for name, (header_lbl, params) in method_params.items():
        labels.append(header_lbl); values.append("")
        for k, v in params:
            labels.append(k); values.append(str(v))
        if name in curves:
            lam, ti, _ = curves[name]
            labels += ["ΔF final [eV]", "RMS vs ref [eV]"]
            values += [f"{ti[-1]:.6f}", f"{rms[name]:.6f}"]
        else:
            labels += ["ΔF final [eV]", "RMS vs ref [eV]"]
            values += ["FAILED", "—"]

    labels += ["── Reference ──", "ΔF final [eV]"]
    values += ["", f"{ref[-1]:.6f}"]

    row_colors = ["#E3F2FD" if i % 2 == 0 else "white" for i in range(len(labels))]
    fig.add_trace(go.Table(
        header=dict(values=["Parameter", "Value"],
                    fill_color="#37474F", font=dict(color="white", size=12), align="left"),
        cells=dict(values=[labels, values],
                   fill_color=[row_colors, row_colors],
                   font=dict(size=11), align="left"),
    ), row=1, col=2)

    fig.update_layout(
        title=f"Entropic spring TI comparison — entropic_spring_{N},  T = {T:.0f} K",
        height=650, template="plotly_white",
        legend=dict(x=0.01, y=0.99), hovermode="x unified",
    )
    fig.write_html(str(out_html))

    # --- Summary ---
    print(f"\nFinal ΔF:")
    for name, (lam, ti, _) in curves.items():
        print(f"  {name:<16}: {ti[-1]:.6f} eV  (RMS vs ref: {rms[name]:.6f} eV)")
    for name in ["CPU TI", "GPU-debug TI", "multi-GPU TI"]:
        if name not in curves:
            print(f"  {name:<16}: FAILED (skipped in output)")
    print(f"  {'Reference':<16}: {ref[-1]:.6f} eV")
    print(f"\nData : {out_dat}")
    print(f"Plot : {out_html}")


if __name__ == "__main__":
    main()

