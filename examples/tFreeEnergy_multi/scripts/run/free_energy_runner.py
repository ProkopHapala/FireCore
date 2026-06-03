#!/usr/bin/env python3

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    HAS_PLOTLY = True
except ImportError:
    HAS_PLOTLY = False

HAS_MM = None
mm = None


SCRIPT_DIR = Path(__file__).resolve().parent


def run_cmd(cmd, cwd, log_path):
    print(f"Running: {' '.join(cmd)}")
    with open(log_path, "w") as f:
        p = subprocess.Popen(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        for line in p.stdout:
            f.write(line)
        return p.wait()


def read_xyz(path: Path):
    lines = path.read_text().splitlines()
    natoms = int(lines[0].strip())
    comment = lines[1] if len(lines) > 1 else ""
    atoms = []
    for line in lines[2:2 + natoms]:
        parts = line.split()
        atoms.append(
            {
                "elem": parts[0],
                "xyz": np.array([float(parts[1]), float(parts[2]), float(parts[3])], dtype=np.float64),
                "tail": parts[4:],
            }
        )
    return natoms, comment, atoms


def write_xyz(path: Path, comment: str, atoms):
    with open(path, "w") as f:
        f.write(f"{len(atoms)}\n")
        f.write(f"{comment}\n")
        for a in atoms:
            x, y, z = a["xyz"]
            tail = (" " + " ".join(a["tail"])) if a["tail"] else ""
            f.write(f"{a['elem']:<2s} {x:16.8f} {y:16.8f} {z:16.8f}{tail}\n")


def load_constraints(path: Path):
    rows = []
    for line in path.read_text().splitlines():
        s = line.strip()
        if (not s) or s.startswith("#"):
            continue
        vals = [float(x) for x in s.split()]
        if len(vals) == 6:
            xyz0 = vals[:3]
            xyz1 = vals[3:6]
        elif len(vals) == 7:
            xyz0 = vals[1:4]
            xyz1 = vals[4:7]
        else:
            raise ValueError(f"Expected 6 or 7 floats per constraint row, got: {line}")
        rows.append((xyz0, xyz1))
    return rows


def _normalize_surf_name(surf_name):
    if surf_name is None:
        return None
    s = str(surf_name).strip()
    if s.lower() in ("", "none", "null", "off"):
        return None
    return s


def _surface_tag(surf_name):
    surf_name = _normalize_surf_name(surf_name)
    if surf_name is None:
        return "vacuum"
    stem = os.path.basename(surf_name.rstrip("/"))
    return re.sub(r"[^a-zA-Z0-9_-]+", "_", stem).strip("_") or "surface"


def _expand_surface_variants(common, param_sets):
    surfaces = common.get("surfaces")
    common_base = dict(common)
    common_base.pop("surfaces", None)

    if not surfaces:
        return common_base, param_sets

    expanded = []
    for pset in param_sets:
        base_name = pset.get("name", "unnamed_run")
        for surf in surfaces:
            if isinstance(surf, str):
                surf_params = {"name": _surface_tag(surf), "surf_name": surf}
            elif isinstance(surf, dict):
                surf_params = dict(surf)
            else:
                raise ValueError(f"Unsupported surface entry: {surf!r}")
            surf_label = surf_params.pop("name", _surface_tag(surf_params.get("surf_name")))
            merged = dict(pset)
            merged.update(surf_params)
            merged["name"] = f"{base_name}__{surf_label}"
            expanded.append(merged)
    return common_base, expanded


def _load_cfg(config_path):
    with open(config_path, "r") as f:
        cfg = json.load(f)
    output_root = cfg.get("output_root", "results_sweep")
    common = cfg.get("common", {})
    param_sets = cfg.get("parameter_sets", [{}])
    common, param_sets = _expand_surface_variants(common, param_sets)
    return cfg, output_root, common, param_sets


def _parse_header_cols(path):
    col_map = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line.startswith("#"):
                break
            if "lambda" in line:
                parts = line.lstrip("#").split()
                for i, p in enumerate(parts):
                    col_map[p] = i
    return col_map


def _load_fe_dat(path):
    col_map = _parse_header_cols(path)
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data[None, :]
    ncols = data.shape[1]

    x_col = min(col_map.get("distance", col_map.get("lambda", 0)), ncols - 1)
    f_col = min(col_map.get("TI_F", 3), ncols - 1)
    e_col = min(col_map.get("TI_err", 4), ncols - 1)
    de_col = min(col_map.get("TI_dE/dlambda", 1), ncols - 1)

    x = data[:, x_col]
    ti = data[:, f_col]
    err = data[:, e_col]
    de = data[:, de_col]
    lam = data[:, col_map.get("lambda", 0)]

    valid = ~np.isnan(ti)
    if np.any(valid):
        ti = ti - ti[valid][0]

    return x, ti, err, de, lam


def perform_scan(mode_name, b_relaxed, nCVs, initial_positions, final_positions, nLambda, niter, Fconv, K, natoms, si_indices, other_si_idx, lambdas, target_pos, target_dist):
    ppos = np.zeros((nLambda, natoms, 3), dtype=np.float32)
    Es, _ = mm.scan_Milan(
        nCVs,
        initial_positions,
        final_positions,
        nLambda,
        nsteps=niter,
        Fconv=Fconv,
        K=K,
        bRelaxed=b_relaxed,
        Es=np.zeros(nLambda, dtype=np.float64),
        ppos=ppos,
    )

    if len(si_indices) > 0 and other_si_idx is not None:
        pulled_pos = ppos[:, si_indices[0], :]
        other_pos = ppos[:, other_si_idx, :]
        distances = np.linalg.norm(pulled_pos - other_pos, axis=1)
        dist_err = np.abs(distances - target_dist)
    else:
        pulled_pos = np.zeros((nLambda, 3))
        other_pos = np.zeros((nLambda, 3))
        distances = target_dist.copy()
        dist_err = np.zeros_like(target_dist)

    return {
        "name": mode_name,
        "Es": Es,
        "pulled_pos": pulled_pos,
        "other_pos": other_pos,
        "distances": distances,
        "dist_err": dist_err,
        "target_pos": target_pos,
        "lambdas": lambdas,
    }


def run_scan_bundle(output_root, common):
    global HAS_MM, mm
    if HAS_MM is None:
        _root_dir = os.path.abspath(os.path.join(SCRIPT_DIR, "../../.."))
        sys.path.append(_root_dir)
        try:
            from pyBall import MMFF_multi as mm_mod
            mm = mm_mod
            HAS_MM = True
        except ImportError:
            HAS_MM = False

    scan_mode = str(common.get("scan_mode", "none")).lower()
    if scan_mode not in ["rigid", "relaxed", "both"] or (not HAS_MM):
        return [], None, []

    xyz_in = Path(common.get("xyz_name", ""))
    constraints_in = Path(common.get("constraints", "constraints.txt"))
    if not xyz_in.exists() or not constraints_in.exists():
        print("Warning: Missing XYZ or constraints file. Skipping potential-energy scans.")
        return [], None, []

    natoms, comment, atoms = read_xyz(xyz_in)
    constraints = load_constraints(constraints_in)
    if len(constraints) == 0:
        print("Warning: No constraints loaded. Skipping potential-energy scans.")
        return [], None, []

    target_start = np.array(constraints[0][0], dtype=np.float64)
    target_end = np.array(constraints[0][1], dtype=np.float64)
    nLambda = int(common.get("nLambda", 31))
    lambdas = np.linspace(0.0, 1.0, nLambda)
    target_pos = (1.0 - lambdas[:, None]) * target_start[None, :] + lambdas[:, None] * target_end[None, :]

    if len(constraints) > 1:
        target_dist = np.linalg.norm(
            ((1.0 - lambdas[:, None]) * np.array(constraints[0][0], dtype=np.float64)[None, :] + lambdas[:, None] * np.array(constraints[0][1], dtype=np.float64)[None, :]) -
            ((1.0 - lambdas[:, None]) * np.array(constraints[1][0], dtype=np.float64)[None, :] + lambdas[:, None] * np.array(constraints[1][1], dtype=np.float64)[None, :]),
            axis=1,
        )
    else:
        target_dist = np.linalg.norm(target_pos, axis=1)

    tmp_xyz = Path(output_root) / "scan_tmp.xyz"
    write_xyz(tmp_xyz, f"{comment} | unified scan mode", atoms)

    surf_name = _normalize_surf_name(common.get("surf_name"))
    init_kwargs = dict(
        nSys_=int(common.get("nSys", 1)),
        xyz_name=str(tmp_xyz),
        sElementTypes="../../cpp/common_resources/ElementTypes.dat",
        sAtomTypes="../../cpp/common_resources/AtomTypes.dat",
        sBondTypes="../../cpp/common_resources/BondTypes.dat",
        sAngleTypes="../../cpp/common_resources/AngleTypes.dat",
        bMMFF=True,
        bEpairs=False,
        nPBC=(0, 0, 0),
        T=0.0,
        gamma=0.0,
    )
    if surf_name is not None:
        init_kwargs["surf_name"] = surf_name
        init_kwargs["GridFF"] = 6

    mm.init(**init_kwargs)
    mm.getBuffs()
    dt = float(common.get("dt", 0.05))
    mm.set_opt(
        dt_max=dt,
        dt_min=max(dt * 0.1, 1.0e-6),
        damp_max=0.2,
        finc=1.1,
        fdec=0.5,
        falpha=0.8,
        minLastNeg=5,
        cvf_min=-0.1,
        cvf_max=0.1,
    )

    si_indices = [ia for ia in range(mm.natoms) if mm.getTypeName(ia).strip() == "Si"]
    other_si_idx = si_indices[1] if len(si_indices) >= 2 else None

    initial_positions = []
    final_positions = []
    for init_pos, final_pos in constraints:
        initial_positions.extend(init_pos)
        final_positions.extend(final_pos)

    nCVs = len(constraints)
    niter = int(common.get("scan_nsteps", common.get("scan_niter", common.get("nMDsteps", 20000))))
    Fconv = float(common.get("scan_Fconv", common.get("Fconv", 1e-4)))
    K = float(common.get("scan_K", common.get("K", 10.0)))

    scan_results = []
    if scan_mode in ["relaxed", "both"]:
        print("  -> Running relaxed scan...")
        scan_results.append(perform_scan("relaxed", True, nCVs, initial_positions, final_positions, nLambda, niter, Fconv, K, mm.natoms, si_indices, other_si_idx, lambdas, target_pos, target_dist))
    if scan_mode in ["rigid", "both"]:
        print("  -> Running rigid scan...")
        scan_results.append(perform_scan("rigid", False, nCVs, initial_positions, final_positions, nLambda, niter, Fconv, K, mm.natoms, si_indices, other_si_idx, lambdas, target_pos, target_dist))
    return scan_results, target_dist, lambdas


def plot_sweep_comparison(output_root, output_html=None):
    if not HAS_PLOTLY:
        print("plotly not installed – skipping comparison plot.")
        return

    if output_html is None:
        output_html = os.path.join(output_root, "comparison_temperatures.html")

    runs = []
    for entry in sorted(os.listdir(output_root)):
        run_dir = os.path.join(output_root, entry)
        if not os.path.isdir(run_dir):
            continue
        dat_files = [f for f in os.listdir(run_dir) if f.endswith(".dat")]
        if not dat_files:
            continue
        m = re.search(r"(\d+(?:\.\d+)?)[Kk]", entry)
        T_val = float(m.group(1)) if m else None
        runs.append((T_val, entry, os.path.join(run_dir, dat_files[0])))

    if not runs:
        print(f"No .dat files found under {output_root} – nothing to plot.")
        return

    runs.sort(key=lambda r: (r[0] is None, r[0] or 0))
    temps = [r[0] for r in runs if r[0] is not None]
    T_min, T_max = (min(temps), max(temps)) if len(temps) > 1 else ((temps[0], temps[0] + 1) if temps else (0, 1))

    def temp_color(T):
        if T is None:
            return "grey"
        t = (T - T_min) / (T_max - T_min) if T_max > T_min else 0.5
        r = int(255 * t)
        b = int(255 * (1 - t))
        return f"rgb({r},0,{b})"

    fig = make_subplots(rows=2, cols=1, subplot_titles=("Cumulative Free Energy ΔF vs distance", "dE/dλ vs λ"), vertical_spacing=0.12)

    for T_val, name, dat_path in runs:
        try:
            x, ti, err, de, lam = _load_fe_dat(dat_path)
        except Exception as exc:
            print(f"  Skipping {name}: {exc}")
            continue

        label = f"{T_val:.0f} K" if T_val is not None else name
        color = temp_color(T_val)
        errarr = err if not np.all(np.isnan(err)) else None
        err_color = None if errarr is None else color.replace("rgb", "rgba").replace(")", ",0.3)")

        fig.add_trace(
            go.Scatter(
                x=x,
                y=ti,
                name=label,
                mode="lines",
                line=dict(color=color, width=2),
                error_y=dict(type="data", array=errarr, visible=(errarr is not None), color=err_color),
            ),
            row=1,
            col=1,
        )
        fig.add_trace(
            go.Scatter(
                x=lam,
                y=de,
                name=label,
                mode="lines",
                line=dict(color=color, width=1.5),
                showlegend=False,
            ),
            row=2,
            col=1,
        )

    fig.update_xaxes(title_text="End-to-end distance [Å]", row=1, col=1)
    fig.update_yaxes(title_text="ΔF [eV]", row=1, col=1)
    fig.update_xaxes(title_text="λ", row=2, col=1)
    fig.update_yaxes(title_text="dE/dλ [eV]", row=2, col=1)
    fig.update_layout(height=900, title=f"Free Energy Comparison – {os.path.basename(output_root)}", template="plotly_white", hovermode="x unified")
    fig.write_html(output_html)
    print(f"Comparison plot saved to: {output_html}")


def plot_combined(output_root, scan_results, target_dist, lambdas, output_html=None):
    if not HAS_PLOTLY:
        print("Plotly not installed - skipping combined plot.")
        return

    if output_html is None:
        output_html = os.path.join(output_root, "temperature_scan.html")

    runs = []
    if os.path.exists(output_root):
        for entry in sorted(os.listdir(output_root)):
            run_dir = os.path.join(output_root, entry)
            if not os.path.isdir(run_dir):
                continue
            dat_files = [f for f in os.listdir(run_dir) if f.endswith(".dat") and "free_energy" in f]
            if not dat_files:
                continue
            m = re.search(r"(\d+(?:\.\d+)?)[Kk]", entry)
            T_val = float(m.group(1)) if m else None
            runs.append((T_val, entry, os.path.join(run_dir, dat_files[0])))

    runs.sort(key=lambda r: (r[0] is None, r[0] or 0))
    temps = [r[0] for r in runs if r[0] is not None]
    if temps:
        T_min, T_max = (min(temps), max(temps)) if len(temps) > 1 else (temps[0], temps[0] + 1)
    else:
        T_min, T_max = 0, 1

    def temp_color(T):
        if T is None:
            return "grey"
        t = (T - T_min) / (T_max - T_min) if T_max > T_min else 0.5
        r = int(255 * t)
        b = int(255 * (1 - t))
        return f"rgb({r},0,{b})"

    fig = go.Figure()

    for sr in scan_results:
        label = f"{sr['name'].capitalize()} Scan"
        color = "black" if sr["name"] == "relaxed" else "gray"
        dash = "solid" if sr["name"] == "relaxed" else "dash"
        Es = sr["Es"]
        valid = ~np.isnan(Es)
        if np.any(valid):
            Es = Es - Es[valid][0]
        fig.add_trace(go.Scatter(x=sr["lambdas"], y=Es, name=label, mode="lines+markers", line=dict(color=color, width=2, dash=dash)))

    for T_val, name, dat_path in runs:
        try:
            _, ti, err, _, lam = _load_fe_dat(dat_path)
        except Exception as exc:
            print(f"  Skipping {name}: {exc}")
            continue
        label = f"TI {T_val:.0f} K" if T_val is not None else f"TI {name}"
        color = temp_color(T_val)
        errarr = err if not np.all(np.isnan(err)) else None
        fig.add_trace(go.Scatter(x=lam, y=ti, name=label, mode="lines", line=dict(color=color, width=2), error_y=dict(type="data", array=errarr, visible=(errarr is not None))))

    tick_lambdas = np.linspace(0.0, 1.0, 11)
    tick_distances = np.interp(tick_lambdas, lambdas, target_dist) if (target_dist is not None and len(lambdas) > 0) else tick_lambdas
    fig.update_layout(
        title="Temperature Scan: Potential Energy vs Free Energy",
        template="plotly_white",
        hovermode="x unified",
        xaxis=dict(title="λ", domain=[0, 1], tickmode="array", tickvals=tick_lambdas, tickformat=".2f"),
        xaxis2=dict(title="Distance [Å]", overlaying="x", side="top", tickmode="array", tickvals=tick_lambdas, ticktext=[f"{d:.2f}" for d in tick_distances]),
        yaxis=dict(title="ΔE / ΔF [eV]"),
        height=800,
    )
    fig.write_html(output_html)
    print(f"Combined plot saved to: {output_html}")


def run_sweep(config_path, include_scans=True, analyze_pairs=None, hbond_cut=None, log_temperature=None):
    _, output_root, common, param_sets = _load_cfg(config_path)
    if analyze_pairs is not None:
        common["analyze_pairs"] = analyze_pairs
    if hbond_cut is not None:
        common["hbond_cut"] = hbond_cut
    if log_temperature is not None:
        common["log_temperature"] = log_temperature

    os.makedirs(output_root, exist_ok=True)
    shutil.copy2(config_path, os.path.join(output_root, "config.json"))

    scan_results = []
    target_dist = None
    lambdas = []
    if include_scans:
        scan_results, target_dist, lambdas = run_scan_bundle(output_root, common)

    summary_file = os.path.join(output_root, "summary.txt")
    with open(summary_file, "w") as sf:
        sf.write("Run Name | Status | Wall Time (s)\n")
        sf.write("-" * 40 + "\n")

    for pset in param_sets:
        params = common.copy()
        params.update(pset)

        run_name = params.get("name", "unnamed_run")
        run_dir = os.path.join(output_root, run_name)
        os.makedirs(run_dir, exist_ok=True)

        print(f"\n>>> Starting run: {run_name}")
        cmd = ["python3", "scripts/run/run_ES.py"]
        for key, val in params.items():
            if key in ["name", "scan_mode", "surfaces", "scan_nsteps"]:
                continue
            if isinstance(val, bool):
                if val:
                    cmd.append(f"--{key}")
            elif key == "temperature":
                cmd.extend(["-T", str(val)])
            else:
                cmd.extend([f"--{key}", str(val)])

        log_path = os.path.join(run_dir, "run.log")
        t0 = time.time()
        rc = run_cmd(cmd, cwd=os.getcwd(), log_path=log_path)
        dt = time.time() - t0
        
        # Plot interactive HTML for this individual run if data was created
        xyz_base = os.path.splitext(os.path.basename(params.get("xyz_name", "")))[0]
        out_dat = f"{xyz_base}_free_energy.dat"
        if rc == 0 and os.path.exists(out_dat):
            plot_cmd = ["python3", "scripts/analysis/plot_F_interactive.py", "--input", out_dat]
            run_cmd(plot_cmd, cwd=os.getcwd(), log_path=log_path + ".plot")

        status = "OK" if rc == 0 else "FAILED"
        print(f"<<< Finished {run_name} with status {status} in {dt:.2f}s")
        with open(summary_file, "a") as sf:
            sf.write(f"{run_name} | {status} | {dt:.2f}\n")

        for f in os.listdir("."):
            if f.endswith(".dat") or f.endswith(".html"):
                if xyz_base in f or "free_energy" in f or run_name in f:
                    target = os.path.join(run_dir, f)
                    if os.path.exists(target):
                        os.remove(target)
                    shutil.move(f, target)

    plot_sweep_comparison(output_root)
    if include_scans and scan_results:
        plot_combined(output_root, scan_results, target_dist, lambdas)


def main():
    parser = argparse.ArgumentParser(description="Unified free-energy sweep runner")
    parser.add_argument("--config", default="default_free_energy_config.json", help="JSON config file")
    parser.add_argument("--no-scans", action="store_true", help="Disable rigid/relaxed scans even if scan_mode is enabled")
    parser.add_argument("--analyze_pairs", type=str, default=None, help="Atom index pairs to measure, overriding config")
    parser.add_argument("--hbond_cut", type=float, default=None, help="Hydrogen-bond cutoff, overriding config")
    parser.add_argument("--log_temperature", action="store_true", default=None, help="Enable measured temperature logging, overriding config")
    args = parser.parse_args()
    run_sweep(
        args.config,
        include_scans=(not args.no_scans),
        analyze_pairs=args.analyze_pairs,
        hbond_cut=args.hbond_cut,
        log_temperature=args.log_temperature
    )


if __name__ == "__main__":
    main()
