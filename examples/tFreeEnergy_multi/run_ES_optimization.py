#!/usr/bin/env python3

import argparse
import csv
import itertools
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional

from analyze_free_energy_accuracy import compute_metrics_from_file


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
BUILD_DIR = os.path.join(PROJECT_ROOT, "cpp", "Build", "libs_OCL")


def _parse_int_list(txt: str) -> List[int]:
    return [int(x.strip()) for x in txt.split(",") if x.strip()]


def _parse_float_list(txt: str) -> List[float]:
    return [float(x.strip()) for x in txt.split(",") if x.strip()]


def _parse_mode_list(txt: str) -> List[str]:
    out = [x.strip().upper() for x in txt.split(",") if x.strip()]
    for m in out:
        if m not in ("TI", "JE", "BOTH"):
            raise ValueError(f"Unsupported mode '{m}', use TI/JE/BOTH")
    return out


def _float_tag(x: Optional[float]) -> str:
    if x is None:
        return "na"
    s = f"{x:.1f}"
    return s.replace("-", "m").replace(".", "p")


def _safe_mkdir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def _run_cmd(cmd: List[str], cwd: str, log_path: str, dry_run: bool = False) -> int:
    if dry_run:
        with open(log_path, "w", encoding="utf-8") as f:
            f.write("[DRY RUN] " + " ".join(cmd) + "\n")
        return 0
    with open(log_path, "w", encoding="utf-8") as f:
        p = subprocess.Popen(
            cmd,
            cwd=cwd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        assert p.stdout is not None
        for line in p.stdout:
            f.write(line)
            sys.stdout.write(line)
        return p.wait()


def _build_lib(dry_run: bool = False) -> None:
    if dry_run:
        print(f"[DRY RUN] build in {BUILD_DIR}")
        return
    subprocess.run(["make", "MMFFmulti_lib"], cwd=BUILD_DIR, check=True)


def _parse_run_ocl_timing(log_text: str) -> Dict[str, float]:
    rgx = re.compile(
        r"run_ocl_opt\(.*?\)\s+NOT CONVERGED.*?time\s+([0-9eE+\-.]+)\s+\[ms\]\(\s*([0-9eE+\-.]+)\s+\[us/step\]",
    )
    phase = "unknown"
    phase_ms: Dict[str, float] = {"equilibrating": 0.0, "production": 0.0, "pulling": 0.0, "unknown": 0.0}
    total_ms = 0.0
    total_calls = 0
    us_per_step_vals: List[float] = []

    for line in log_text.splitlines():
        if "Equilibrating " in line:
            phase = "equilibrating"
        elif "Production " in line:
            phase = "production"
        elif "Pulling for " in line:
            phase = "pulling"

        m = rgx.search(line)
        if not m:
            continue
        t_ms = float(m.group(1))
        us_per_step = float(m.group(2))
        total_ms += t_ms
        total_calls += 1
        us_per_step_vals.append(us_per_step)
        phase_ms[phase] = phase_ms.get(phase, 0.0) + t_ms

    avg_us_per_step = float("nan")
    if us_per_step_vals:
        avg_us_per_step = float(sum(us_per_step_vals) / len(us_per_step_vals))
    return {
        "md_ocl_time_ms_total": total_ms,
        "md_ocl_calls": float(total_calls),
        "md_ocl_avg_us_per_step": avg_us_per_step,
        "md_ocl_time_ms_equilibrating": phase_ms.get("equilibrating", 0.0),
        "md_ocl_time_ms_production": phase_ms.get("production", 0.0),
        "md_ocl_time_ms_pulling": phase_ms.get("pulling", 0.0),
        "md_ocl_time_ms_unknown": phase_ms.get("unknown", 0.0),
    }


@dataclass
class Combo:
    mode: str
    nSys: int
    nLambda: int
    nEQsteps: int
    nMDsteps: int
    JE_K: Optional[float]
    dt: float
    t_damp: float
    je_target_trajectories: Optional[int] = None
    n_list: Optional[List[int]] = None
    label: str = ""

    def tag(self) -> str:
        tag_parts = []
        if self.label:
            safe_label = re.sub(r"[^a-zA-Z0-9_-]+", "_", self.label).strip("_")
            if safe_label:
                tag_parts.append(safe_label)
        parts = [
            f"nsys{self.nSys}",
            f"nl{self.nLambda}",
            f"neq{self.nEQsteps}",
            f"nmd{self.nMDsteps}",
            f"dt{_float_tag(self.dt)}",
            f"td{_float_tag(self.t_damp)}",
        ]
        if self.mode == "JE":
            parts.append(f"jek{_float_tag(self.JE_K)}")
        if self.mode == "BOTH":
            parts.append(f"jek{_float_tag(self.JE_K)}")
        tag_parts.extend(parts)
        return "_".join(tag_parts)


def _all_combos(
    modes: Iterable[str],
    nsys_list: Iterable[int],
    nlambda_list: Iterable[int],
    neq_list: Iterable[int],
    nmd_list: Iterable[int],
    je_k_list: Iterable[float],
    dt: float,
    t_damp: float,
) -> List[Combo]:
    out: List[Combo] = []
    for mode in modes:
        for nsys, nl, neq, nmd in itertools.product(nsys_list, nlambda_list, neq_list, nmd_list):
            if mode in ("JE", "BOTH"):
                for jk in je_k_list:
                    out.append(Combo(mode, nsys, nl, neq, nmd, jk, dt, t_damp))
            else:
                out.append(Combo(mode, nsys, nl, neq, nmd, None, dt, t_damp))
    return out


def _to_int(v, key: str) -> int:
    try:
        return int(v)
    except Exception as e:
        raise ValueError(f"Invalid int for '{key}': {v}") from e


def _to_float(v, key: str) -> float:
    try:
        return float(v)
    except Exception as e:
        raise ValueError(f"Invalid float for '{key}': {v}") from e


def _combos_from_param_sets(param_sets: List[dict], common: Dict[str, object]) -> List[Combo]:
    out: List[Combo] = []
    for i, ps in enumerate(param_sets):
        mode = str(ps.get("mode", common.get("mode", ""))).upper()
        if mode not in ("TI", "JE", "BOTH"):
            raise ValueError(f"Set[{i}] has invalid mode '{mode}', expected TI/JE/BOTH")
        nsys = _to_int(ps.get("nSys", common.get("nSys")), f"set[{i}].nSys")
        nlambda = _to_int(ps.get("nLambda", common.get("nLambda")), f"set[{i}].nLambda")
        neq = _to_int(ps.get("nEQsteps", common.get("nEQsteps")), f"set[{i}].nEQsteps")
        nmd_val = ps.get("nMDsteps", common.get("nMDsteps"))
        je_traj_val = ps.get("je_target_trajectories", common.get("je_target_trajectories"))
        je_target_trajectories: Optional[int] = None
        if je_traj_val is not None:
            je_target_trajectories = _to_int(je_traj_val, f"set[{i}].je_target_trajectories")

        if nmd_val is None:
            if mode in ("JE", "BOTH") and je_target_trajectories is not None:
                nmd = int(nlambda * je_target_trajectories)
            else:
                raise ValueError(f"Set[{i}] must provide nMDsteps (or je_target_trajectories for JE/BOTH)")
        else:
            nmd = _to_int(nmd_val, f"set[{i}].nMDsteps")
        dt = _to_float(ps.get("dt", common.get("dt", 0.05)), f"set[{i}].dt")
        t_damp = _to_float(ps.get("t_damp", common.get("t_damp", 150.0)), f"set[{i}].t_damp")

        je_k = ps.get("JE_K", common.get("JE_K"))
        if mode in ("JE", "BOTH"):
            if je_k is None:
                raise ValueError(f"Set[{i}] mode={mode} requires JE_K")
            je_k = _to_float(je_k, f"set[{i}].JE_K")
        else:
            je_k = None
        label = str(ps.get("name", "")).strip()
        out.append(
            Combo(
                mode=mode,
                nSys=nsys,
                nLambda=nlambda,
                nEQsteps=neq,
                nMDsteps=nmd,
                JE_K=je_k,
                dt=dt,
                t_damp=t_damp,
                je_target_trajectories=je_target_trajectories,
                n_list=[_to_int(x, f"set[{i}].N_list[]") for x in ps.get("N_list", [])] if "N_list" in ps else None,
                label=label,
            )
        )
    return out


def _has_numeric_rows(path: str) -> bool:
    if not os.path.isfile(path):
        return False
    try:
        with open(path, "r", encoding="utf-8") as f:
            for line in f:
                s = line.strip()
                if not s or s.startswith("#"):
                    continue
                return True
    except Exception:
        return False
    return False


def _is_completed_run(run_dir: str, base_name: str) -> bool:
    data_file = os.path.join(run_dir, f"{base_name}_free_energy.dat")
    plot_file = os.path.join(run_dir, f"{base_name}_free_energy_F_interactive.html")
    metrics_file = os.path.join(run_dir, "metrics.json")
    if not (_has_numeric_rows(data_file) and os.path.isfile(plot_file) and os.path.isfile(metrics_file)):
        return False
    try:
        with open(metrics_file, "r", encoding="utf-8") as f:
            m = json.load(f)
        return str(m.get("status", "")).lower() == "ok"
    except Exception:
        return False


def _load_config(path: str) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        cfg = json.load(f)
    if not isinstance(cfg, dict):
        raise ValueError("Config root must be a JSON object")
    return cfg


def _is_dominated(a: Dict[str, float], b: Dict[str, float]) -> bool:
    # Minimize both wall time and RMSD
    ta, ra = a["wall_s"], a["rmsd"]
    tb, rb = b["wall_s"], b["rmsd"]
    better_or_equal = (tb <= ta) and (rb <= ra)
    strictly_better = (tb < ta) or (rb < ra)
    return better_or_equal and strictly_better


def _write_csv(path: str, rows: List[Dict[str, object]], fieldnames: List[str]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def _write_references(ref_dir: str, n_list: List[int], n_lambda: int, T: float, b: float) -> None:
    _safe_mkdir(ref_dir)
    kB = 8.617333262145e-5
    lam = [i / (n_lambda - 1) for i in range(n_lambda)]
    # Constraints in constraints_ES.txt imply end-to-end distance 1A -> 11A.
    d0, d1 = 1.0, 11.0
    dist = [d0 + (d1 - d0) * x for x in lam]
    for n_atoms in n_list:
        nseg = n_atoms - 1
        if nseg <= 0:
            continue
        k_spring = kB * T / (nseg * (b ** 2))
        fe_abs = [0.5 * k_spring * d * d for d in dist]
        fe0 = fe_abs[0]
        out = os.path.join(ref_dir, f"entropic_spring_{n_atoms}_ref.dat")
        with open(out, "w", encoding="utf-8") as f:
            f.write("# lambda distance FE_ref_shifted force_ref\n")
            for l, d, fa in zip(lam, dist, fe_abs):
                f.write(f"{l:10.6f} {d:12.8f} {fa-fe0:14.10f} {k_spring*d:14.10f}\n")


def _write_summary_plots(summary_dir: str, rows: List[Dict[str, object]]) -> None:
    _safe_mkdir(summary_dir)
    try:
        import plotly.express as px
    except Exception:
        return

    ok = []
    for r in rows:
        if r.get("status") != "ok":
            continue
        wall = float(r.get("wall_s", float("nan")))
        rmsd = float(r.get("rmsd_profile_eV", float("nan")))
        ferr = float(r.get("final_abs_error_eV", float("nan")))
        wig = float(r.get("je_wiggle_rmse_eV", float("nan")))
        rr = dict(r)
        rr["wall_s"] = wall
        rr["rmsd_profile_eV"] = rmsd
        rr["final_abs_error_eV"] = ferr
        rr["je_wiggle_rmse_eV"] = wig
        ok.append(rr)

    if not ok:
        return

    fig1 = px.scatter(
        ok,
        x="wall_s",
        y="rmsd_profile_eV",
        color="mode",
        hover_data=["N", "nSys", "nLambda", "nEQsteps", "nMDsteps", "JE_K", "run_dir"],
        title="Pareto: Wall Time vs RMSD",
    )
    fig1.write_html(os.path.join(summary_dir, "pareto_time_vs_rmsd.html"))

    fig2 = px.scatter(
        ok,
        x="wall_s",
        y="final_abs_error_eV",
        color="mode",
        hover_data=["N", "nSys", "nLambda", "nEQsteps", "nMDsteps", "JE_K", "run_dir"],
        title="Pareto: Wall Time vs Final |DeltaF-Ref|",
    )
    fig2.write_html(os.path.join(summary_dir, "pareto_time_vs_final_error.html"))

    je = [r for r in ok if r.get("mode") in ("JE", "BOTH") and math.isfinite(r.get("je_wiggle_rmse_eV", float("nan")))]
    if je:
        fig3 = px.scatter(
            je,
            x="wall_s",
            y="je_wiggle_rmse_eV",
            color="N",
            symbol="mode",
            hover_data=["nSys", "nLambda", "nEQsteps", "nMDsteps", "JE_K", "run_dir"],
            title="JE Wiggle Fit RMSE on lambda in [0,0.05]",
        )
        fig3.write_html(os.path.join(summary_dir, "je_wiggle_fit_lambda_0_0p05.html"))


def main() -> None:
    ap = argparse.ArgumentParser(description="Parameter sweep for TI/JE entropic spring free-energy runs.")
    ap.add_argument("--N-list", default="10,20,30", help="Comma-separated atom counts (e.g. 10,20,30)")
    ap.add_argument("--modes", default="TI,JE", help="Comma-separated: TI,JE")
    ap.add_argument("--nSys-list", default="100", help="Comma-separated nSys values")
    ap.add_argument("--nLambda-list", default="1000", help="Comma-separated nLambda values")
    ap.add_argument("--nEQsteps-list", default="5000", help="Comma-separated nEQsteps values")
    ap.add_argument("--nMDsteps-list", default="1000000", help="Comma-separated nMDsteps values")
    ap.add_argument("--je-k-list", default="2.0", help="Comma-separated JE_K values (JE mode only)")
    ap.add_argument("--output-root", default="bench_ES", help="Output directory under examples/tFreeEnergy_multi")
    ap.add_argument("--constraints", default="constraints_ES.txt")
    ap.add_argument("--Fconv", type=float, default=1e-6)
    ap.add_argument("--dt", type=float, default=0.05)
    ap.add_argument("--temperature", type=float, default=300.0)
    ap.add_argument("--t_damp", type=float, default=150.0)
    ap.add_argument("--segment-b", type=float, default=1.198, help="Entropic spring segment length for reference")
    ap.add_argument("--wiggle-lambda-max", type=float, default=0.05)
    ap.add_argument("--ref-nLambda", type=int, default=1000, help="Number of lambda points in reference files")
    ap.add_argument("--config", default=None, help="JSON config file with explicit parameter sets")
    ap.add_argument("--skip-build", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--stop-on-error", action="store_true")
    ap.add_argument("--skip-existing", action="store_true", help="Skip parameter sets already computed in output folder")
    args = ap.parse_args()

    if args.config:
        cfg = _load_config(args.config)
        common_cfg = cfg.get("common", {})
        if not isinstance(common_cfg, dict):
            raise ValueError("'common' in config must be an object")

        if "output_root" in cfg:
            args.output_root = str(cfg["output_root"])
        if "constraints" in common_cfg:
            args.constraints = str(common_cfg["constraints"])
        if "Fconv" in common_cfg:
            args.Fconv = float(common_cfg["Fconv"])
        if "dt" in common_cfg:
            args.dt = float(common_cfg["dt"])
        if "temperature" in common_cfg:
            args.temperature = float(common_cfg["temperature"])
        if "t_damp" in common_cfg:
            args.t_damp = float(common_cfg["t_damp"])
        if "segment_b" in common_cfg:
            args.segment_b = float(common_cfg["segment_b"])
        if "wiggle_lambda_max" in common_cfg:
            args.wiggle_lambda_max = float(common_cfg["wiggle_lambda_max"])
        if "ref_nLambda" in common_cfg:
            args.ref_nLambda = int(common_cfg["ref_nLambda"])
        if "skip_build" in common_cfg:
            args.skip_build = bool(common_cfg["skip_build"])
        if "stop_on_error" in common_cfg:
            args.stop_on_error = bool(common_cfg["stop_on_error"])
        if "skip_existing" in common_cfg:
            args.skip_existing = bool(common_cfg["skip_existing"])

        n_list = [int(x) for x in cfg.get("N_list", _parse_int_list(args.N_list))]

        param_sets = cfg.get("parameter_sets")
        if param_sets is None:
            # Backward-compatible config style
            param_sets = []
            for s in cfg.get("ti_sets", []):
                x = dict(s)
                x["mode"] = "TI"
                param_sets.append(x)
            for s in cfg.get("je_sets", []):
                x = dict(s)
                x["mode"] = "JE"
                param_sets.append(x)
            for s in cfg.get("both_sets", []):
                x = dict(s)
                x["mode"] = "BOTH"
                param_sets.append(x)
        if not isinstance(param_sets, list) or len(param_sets) == 0:
            raise ValueError("Config must provide non-empty 'parameter_sets' (or ti_sets/je_sets/both_sets).")
        combos = _combos_from_param_sets(param_sets, common_cfg)
        modes = sorted(list({c.mode for c in combos}))
    else:
        n_list = _parse_int_list(args.N_list)
        modes = _parse_mode_list(args.modes)
        nsys_list = _parse_int_list(args.nSys_list)
        nlambda_list = _parse_int_list(args.nLambda_list)
        neq_list = _parse_int_list(args.nEQsteps_list)
        nmd_list = _parse_int_list(args.nMDsteps_list)
        je_k_list = _parse_float_list(args.je_k_list)
        combos = _all_combos(modes, nsys_list, nlambda_list, neq_list, nmd_list, je_k_list, args.dt, args.t_damp)

    out_root = os.path.join(SCRIPT_DIR, args.output_root)
    runs_root = os.path.join(out_root, "runs")
    summary_root = os.path.join(out_root, "summary")
    refs_root = os.path.join(out_root, "references")
    _safe_mkdir(runs_root)
    _safe_mkdir(summary_root)
    _safe_mkdir(refs_root)
    _write_references(refs_root, n_list, args.ref_nLambda, args.temperature, args.segment_b)

    if not args.skip_build:
        print("Building MMFFmulti_lib once before sweep...")
        _build_lib(dry_run=args.dry_run)

    all_rows: List[Dict[str, object]] = []
    total_jobs = 0
    for n_atoms in n_list:
        for combo in combos:
            if combo.n_list is not None and len(combo.n_list) > 0 and n_atoms not in combo.n_list:
                continue
            total_jobs += 1
    job_i = 0

    for n_atoms in n_list:
        for combo in combos:
            if combo.n_list is not None and len(combo.n_list) > 0 and n_atoms not in combo.n_list:
                continue
            job_i += 1
            run_tag = combo.tag()
            run_dir = os.path.join(runs_root, f"N{n_atoms}", combo.mode, run_tag)
            _safe_mkdir(run_dir)
            print(f"\n=== [{job_i}/{total_jobs}] N={n_atoms} mode={combo.mode} {run_tag} ===")

            base_name = f"entropic_spring_{n_atoms}"
            cwd_data = os.path.join(SCRIPT_DIR, f"{base_name}_free_energy.dat")
            cwd_work = os.path.join(SCRIPT_DIR, "jarzynski_work.dat")

            if args.skip_existing and _is_completed_run(run_dir, base_name):
                print("  Skipping already completed run.")
                continue

            run_log_path = os.path.join(run_dir, "run.log")
            cmd = [
                "python3",
                "run_ES.py",
                "--mode",
                combo.mode,
                "--nSys",
                str(combo.nSys),
                "--xyz_name",
                f"../tMMFF/data/{base_name}.xyz",
                "--nLambda",
                str(combo.nLambda),
                "--nMDsteps",
                str(combo.nMDsteps),
                "--nEQsteps",
                str(combo.nEQsteps),
                "--Fconv",
                str(args.Fconv),
                "--constraints",
                args.constraints,
                "--JEforceconst",
                str(combo.JE_K if combo.JE_K is not None else 2.0),
                "--dt",
                str(combo.dt),
                "-T",
                str(args.temperature),
                "--t_damp",
                str(combo.t_damp),
            ]

            t0 = time.time()
            rc = _run_cmd(cmd, cwd=SCRIPT_DIR, log_path=run_log_path, dry_run=args.dry_run)
            wall_s = time.time() - t0

            row: Dict[str, object] = {
                "status": "ok" if rc == 0 else "failed",
                "N": n_atoms,
                "mode": combo.mode,
                "param_set": combo.label,
                "nSys": combo.nSys,
                "nLambda": combo.nLambda,
                "nEQsteps": combo.nEQsteps,
                "nMDsteps": combo.nMDsteps,
                "JE_K": combo.JE_K if combo.JE_K is not None else "",
                "dt": combo.dt,
                "t_damp": combo.t_damp,
                "je_target_trajectories": combo.je_target_trajectories if combo.je_target_trajectories is not None else "",
                "je_effective_trajectories": (combo.nMDsteps / combo.nLambda) if combo.mode in ("JE", "BOTH") and combo.nLambda > 0 else "",
                "wall_s": wall_s,
                "run_dir": run_dir,
                "run_log": run_log_path,
            }

            log_txt = ""
            if os.path.isfile(run_log_path):
                with open(run_log_path, "r", encoding="utf-8") as f:
                    log_txt = f.read()
            row.update(_parse_run_ocl_timing(log_txt))

            if rc != 0:
                all_rows.append(row)
                if args.stop_on_error:
                    break
                continue

            if not args.dry_run and not os.path.isfile(cwd_data):
                row["status"] = "failed_missing_data"
                all_rows.append(row)
                if args.stop_on_error:
                    break
                continue

            if not args.dry_run:
                # Keep canonical data filename in each run dir.
                run_data_path = os.path.join(run_dir, f"{base_name}_free_energy.dat")
                shutil.copy2(cwd_data, run_data_path)
                row["data_file"] = run_data_path

                if os.path.isfile(cwd_work):
                    shutil.copy2(cwd_work, os.path.join(run_dir, f"{base_name}_{combo.mode}_{run_tag}_jarzynski_work.dat"))

                if os.path.getsize(run_data_path) == 0:
                    row["status"] = "failed_empty_data"
                    with open(os.path.join(run_dir, "metrics.json"), "w", encoding="utf-8") as f:
                        json.dump(row, f, indent=2)
                    all_rows.append(row)
                    if args.stop_on_error:
                        break
                    continue

                plot_prefix = os.path.join(run_dir, f"{base_name}_free_energy")
                plot_cmd = [
                    "python3",
                    "plot_F_interactive.py",
                    "--input",
                    run_data_path,
                    "--output",
                    plot_prefix,
                    "--N",
                    str(n_atoms),
                    "--T",
                    str(args.temperature),
                    "--b",
                    str(args.segment_b),
                ]
                plot_log_path = os.path.join(run_dir, "plot.log")
                prc = _run_cmd(plot_cmd, cwd=SCRIPT_DIR, log_path=plot_log_path, dry_run=False)
                if prc == 0:
                    row["plot_file"] = f"{plot_prefix}_F_interactive.html"
                else:
                    row["status"] = "failed_plot"

                metric_mode = combo.mode if combo.mode in ("TI", "JE") else "TI"
                try:
                    metrics = compute_metrics_from_file(
                        input_file=run_data_path,
                        mode=metric_mode,
                        n_atoms=n_atoms,
                        temperature=args.temperature,
                        b_segment=args.segment_b,
                        lambda_fit_max=args.wiggle_lambda_max,
                    )
                    row.update(
                        {
                            "rmsd_profile_eV": metrics["rmsd_profile_eV"],
                            "final_abs_error_eV": metrics["final_abs_error_eV"],
                            "final_value_eV": metrics["final_value_eV"],
                            "final_ref_eV": metrics["final_ref_eV"],
                            "je_wiggle_slope_eV_per_lambda": metrics["je_wiggle_slope_eV_per_lambda"],
                            "je_wiggle_r2": metrics["je_wiggle_r2"],
                            "je_wiggle_rmse_eV": metrics["je_wiggle_rmse_eV"],
                        }
                    )
                except Exception:
                    row["status"] = "failed_metrics"

                with open(os.path.join(run_dir, "metrics.json"), "w", encoding="utf-8") as f:
                    json.dump(row, f, indent=2)

            all_rows.append(row)

        if args.stop_on_error and any(r.get("status", "").startswith("failed") for r in all_rows):
            break

    fieldnames = [
        "status",
        "N",
        "mode",
        "param_set",
        "nSys",
        "nLambda",
        "nEQsteps",
        "nMDsteps",
        "JE_K",
        "dt",
        "t_damp",
        "je_target_trajectories",
        "je_effective_trajectories",
        "wall_s",
        "md_ocl_time_ms_total",
        "md_ocl_calls",
        "md_ocl_avg_us_per_step",
        "md_ocl_time_ms_equilibrating",
        "md_ocl_time_ms_production",
        "md_ocl_time_ms_pulling",
        "md_ocl_time_ms_unknown",
        "rmsd_profile_eV",
        "final_abs_error_eV",
        "final_value_eV",
        "final_ref_eV",
        "je_wiggle_slope_eV_per_lambda",
        "je_wiggle_r2",
        "je_wiggle_rmse_eV",
        "run_dir",
        "run_log",
        "data_file",
        "plot_file",
    ]
    all_csv = os.path.join(summary_root, "all_runs.csv")
    _write_csv(all_csv, all_rows, fieldnames)

    ok_rows = [r for r in all_rows if r.get("status") == "ok" and isinstance(r.get("rmsd_profile_eV"), (float, int))]
    best_rows: List[Dict[str, object]] = []
    for n_atoms in n_list:
        for mode in modes:
            subset = [r for r in ok_rows if r.get("N") == n_atoms and r.get("mode") == mode]
            if not subset:
                continue
            subset = [r for r in subset if math.isfinite(float(r.get("rmsd_profile_eV", float("nan"))))]
            if not subset:
                continue
            subset.sort(key=lambda r: (float(r["rmsd_profile_eV"]), float(r["wall_s"])))
            best_rows.append(subset[0])
    _write_csv(os.path.join(summary_root, "best_by_method.csv"), best_rows, fieldnames)

    pareto_rows: List[Dict[str, object]] = []
    for n_atoms in n_list:
        for mode in modes:
            subset = [r for r in ok_rows if r.get("N") == n_atoms and r.get("mode") == mode]
            points = []
            for r in subset:
                wall = float(r.get("wall_s", float("nan")))
                rmsd = float(r.get("rmsd_profile_eV", float("nan")))
                if math.isfinite(wall) and math.isfinite(rmsd):
                    points.append({"wall_s": wall, "rmsd": rmsd, "row": r})
            for a in points:
                dominated = any(_is_dominated(a, b) for b in points if b is not a)
                if not dominated:
                    pareto_rows.append(a["row"])
    _write_csv(os.path.join(summary_root, "best_pareto.csv"), pareto_rows, fieldnames)
    _write_summary_plots(os.path.join(summary_root, "plots"), all_rows)

    try:
        import pandas as pd
        pd.DataFrame(all_rows).to_parquet(os.path.join(summary_root, "all_runs.parquet"), index=False)
    except Exception:
        pass

    print("\nSweep finished.")
    print(f"All runs:      {all_csv}")
    print(f"Best per mode: {os.path.join(summary_root, 'best_by_method.csv')}")
    print(f"Pareto:        {os.path.join(summary_root, 'best_pareto.csv')}")


if __name__ == "__main__":
    main()
