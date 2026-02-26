#!/usr/bin/env python3

import argparse
import csv
import glob
import json
import math
import os
import re
from typing import Dict, List

kB = 8.617333262145e-5  # eV/K


def _write_csv(path: str, rows: List[Dict[str, object]], fieldnames: List[str]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def _is_dominated(a: Dict[str, float], b: Dict[str, float]) -> bool:
    ta, ra = a["wall_s"], a["rmsd"]
    tb, rb = b["wall_s"], b["rmsd"]
    return (tb <= ta and rb <= ra) and (tb < ta or rb < ra)


def _write_summary_plots(summary_dir: str, rows: List[Dict[str, object]]) -> None:
    os.makedirs(summary_dir, exist_ok=True)
    try:
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
    except Exception:
        return

    ok: List[Dict[str, object]] = []
    for r in rows:
        if r.get("status") != "ok":
            continue
        rr = dict(r)
        rr["wall_s"] = float(rr.get("wall_s", float("nan")))
        rr["rmsd_profile_eV"] = float(rr.get("rmsd_profile_eV", float("nan")))
        rr["final_abs_error_eV"] = float(rr.get("final_abs_error_eV", float("nan")))
        rr["je_wiggle_rmse_eV"] = float(rr.get("je_wiggle_rmse_eV", float("nan")))
        ok.append(rr)
    if not ok:
        return

    # Marker style by N, as requested.
    shape_map = {10: "circle", 20: "x", 30: "star"}
    color_map = {"TI": "#1f77b4", "JE": "#d62728", "BOTH": "#2ca02c"}

    # Build per-N 5% reference thresholds analytically (robust to corrupted rows).
    # Reference uses constraints_ES.txt endpoints distance 1A -> 11A.
    ref_thr: Dict[int, float] = {}
    d0, d1 = 1.0, 11.0
    b_segment = 1.198
    temperature = 300.0
    for n in sorted({int(r["N"]) for r in ok}):
        nseg = n - 1
        if nseg <= 0:
            continue
        k_spring = kB * temperature / (nseg * (b_segment ** 2))
        final_ref = 0.5 * k_spring * (d1 * d1 - d0 * d0)
        if final_ref > 0.0 and math.isfinite(final_ref):
            ref_thr[n] = 0.05 * final_ref

    fig = make_subplots(
        rows=1,
        cols=3,
        subplot_titles=(
            "Wall Time vs RMSD",
            "Wall Time vs Final |DeltaF-Ref|",
            "Wall Time vs JE Wiggle RMSE (lambda 0-0.05)",
        ),
        horizontal_spacing=0.08,
    )

    # Keep unique point ids for linked hover across all traces.
    for idx, r in enumerate(ok):
        r["_pid"] = idx

    def mk_hover(r: Dict[str, object], yname: str, yval: float) -> str:
        return (
            f"set={r.get('param_set','')}<br>"
            f"mode={r.get('mode')} N={r.get('N')}<br>"
            f"nSys={r.get('nSys')} nLambda={r.get('nLambda')}<br>"
            f"nEQ={r.get('nEQsteps')} nMD={r.get('nMDsteps')} JE_K={r.get('JE_K')}<br>"
            f"dt={r.get('dt')} t_damp={r.get('t_damp')}<br>"
            f"wall_s={float(r.get('wall_s', float('nan'))):.6g}<br>"
            f"{yname}={yval:.6g}<br>"
            f"run_dir={r.get('run_dir')}<br>"
            f"(click point to open run plot)"
        )

    # Three scatter panels.
    panels = [
        ("rmsd_profile_eV", "RMSD [eV]", 1),
        ("final_abs_error_eV", "Final |DeltaF-Ref| [eV]", 2),
        ("je_wiggle_rmse_eV", "JE Wiggle RMSE [eV]", 3),
    ]

    for key, ylabel, col in panels:
        for mode in ("TI", "JE", "BOTH"):
            for n in sorted({int(x["N"]) for x in ok}):
                xs, ys, cds, htxt = [], [], [], []
                for r in ok:
                    if r.get("mode") != mode or int(r["N"]) != n:
                        continue
                    yv = float(r.get(key, float("nan")))
                    xv = float(r.get("wall_s", float("nan")))
                    if not (math.isfinite(xv) and math.isfinite(yv) and yv > 0.0):
                        continue
                    xs.append(xv)
                    ys.append(yv)
                    pfile = str(r.get("plot_file", ""))
                    purl = f"file://{pfile}" if pfile else ""
                    cds.append([int(r["_pid"]), purl])
                    htxt.append(mk_hover(r, key, yv))
                if not xs:
                    continue
                fig.add_trace(
                    go.Scatter(
                        x=xs,
                        y=ys,
                        mode="markers",
                        name=f"{mode} N{n}",
                        legendgroup=f"{mode}_N{n}",
                        marker=dict(
                            symbol=shape_map.get(n, "circle"),
                            size=10,
                            color=color_map.get(mode, "#444"),
                            opacity=0.9,
                            line=dict(width=1, color="#333"),
                        ),
                        customdata=cds,
                        hovertext=htxt,
                        hoverinfo="text",
                        showlegend=(col == 1),
                    ),
                    row=1,
                    col=col,
                )

        # 5% reference lines (one per N) with dash style linked to marker shape.
        xvals = [float(r.get("wall_s", float("nan"))) for r in ok if math.isfinite(float(r.get("wall_s", float("nan"))))]
        if xvals:
            xmin, xmax = min(xvals), max(xvals)
            dash_by_n = {10: "dot", 20: "dash", 30: "dashdot"}
            for n, thr in sorted(ref_thr.items()):
                if not (math.isfinite(thr) and thr > 0):
                    continue
                fig.add_trace(
                    go.Scatter(
                        x=[xmin, xmax],
                        y=[thr, thr],
                        mode="lines",
                        name=f"N{n} 5% ref",
                        legendgroup=f"refN{n}",
                        line=dict(color="#666", width=1.5, dash=dash_by_n.get(n, "dot")),
                        hoverinfo="skip",
                        showlegend=(col == 1),
                    ),
                    row=1,
                    col=col,
                )

        fig.update_xaxes(title_text="Wall time [s]", row=1, col=col)
        fig.update_yaxes(title_text=ylabel, type="log", row=1, col=col)

    fig.update_layout(
        title="Free-Energy Summary Dashboard",
        height=640,
        width=1900,
        hovermode="closest",
    )

    # Linked-hover behavior: highlight same point id in all subplots, dim others.
    post_script = """
const gd = document.getElementById('{plot_id}');
function setOpacities(activePid){
  const update = {};
  for(let i=0;i<gd.data.length;i++){
    const tr = gd.data[i];
    if(!tr.customdata || !tr.marker || tr.mode!=='markers'){ continue; }
    const ops = [];
    for(let j=0;j<tr.customdata.length;j++){
      const pid = tr.customdata[j][0];
      ops.push((activePid===null || pid===activePid) ? 0.95 : 0.12);
    }
    update[`marker.opacity[${i}]`] = [ops];
  }
  if(Object.keys(update).length>0){ Plotly.restyle(gd, update); }
}
gd.on('plotly_hover', function(evt){
  if(!evt.points || evt.points.length===0){ return; }
  const pt = evt.points[0];
  if(!pt.customdata){ return; }
  setOpacities(pt.customdata[0]);
});
gd.on('plotly_unhover', function(){ setOpacities(null); });
gd.on('plotly_click', function(evt){
  if(!evt.points || evt.points.length===0){ return; }
  const pt = evt.points[0];
  if(!pt.customdata || !pt.customdata[1]){ return; }
  window.open(pt.customdata[1], '_blank');
});

// Build references grouped by legend label and place them under the plot.
const refs = {};
for(let i=0;i<gd.data.length;i++){
  const tr = gd.data[i];
  if(!tr.customdata || !tr.name || tr.mode!=='markers'){ continue; }
  if(!refs[tr.name]) refs[tr.name] = [];
  const seen = new Set(refs[tr.name].map(x => x.url));
  for(let j=0;j<tr.customdata.length;j++){
    const url = tr.customdata[j][1];
    if(!url || seen.has(url)){ continue; }
    seen.add(url);
    const parts = url.split('/');
    const label = parts.length >= 2 ? parts[parts.length-2] : url;
    refs[tr.name].push({url, label});
  }
}
const wrap = document.createElement('div');
wrap.style.marginTop = '10px';
wrap.style.fontFamily = 'sans-serif';
const h = document.createElement('h4');
h.textContent = 'References';
h.style.margin = '4px 0 8px 0';
wrap.appendChild(h);
Object.keys(refs).sort().forEach(name => {
  const row = document.createElement('div');
  row.style.margin = '4px 0';
  const b = document.createElement('span');
  b.textContent = name + ': ';
  b.style.fontWeight = '600';
  row.appendChild(b);
  refs[name].forEach((it, k) => {
    const a = document.createElement('a');
    a.href = it.url;
    a.target = '_blank';
    a.rel = 'noopener noreferrer';
    a.textContent = it.label;
    a.style.marginRight = '10px';
    row.appendChild(a);
  });
  wrap.appendChild(row);
});
gd.parentNode.appendChild(wrap);
"""

    combo_html = os.path.join(summary_dir, "summary_dashboard.html")
    fig.write_html(combo_html, post_script=post_script)

    # Keep backwards-compatible files as lightweight copies of the combined dashboard.
    for name in (
        "pareto_time_vs_rmsd.html",
        "pareto_time_vs_final_error.html",
        "je_wiggle_fit_lambda_0_0p05.html",
    ):
        with open(combo_html, "r", encoding="utf-8") as fin, open(os.path.join(summary_dir, name), "w", encoding="utf-8") as fout:
            fout.write(fin.read())


def _parse_combo(tag: str) -> Dict[str, object]:
    out: Dict[str, object] = {
        "param_set": "",
        "nSys": "",
        "nLambda": "",
        "nEQsteps": "",
        "nMDsteps": "",
        "JE_K": "",
        "dt": "",
        "t_damp": "",
        "je_target_trajectories": "",
        "je_effective_trajectories": "",
    }
    if "_nsys" in tag:
        out["param_set"] = tag.split("_nsys", 1)[0]
    m = re.search(r"nsys(\d+)", tag)
    if m:
        out["nSys"] = int(m.group(1))
    m = re.search(r"_nl(\d+)", tag)
    if m:
        out["nLambda"] = int(m.group(1))
    m = re.search(r"_neq(\d+)", tag)
    if m:
        out["nEQsteps"] = int(m.group(1))
    m = re.search(r"_nmd(\d+)", tag)
    if m:
        out["nMDsteps"] = int(m.group(1))
    m = re.search(r"_jek([0-9mp]+)", tag)
    if m:
        out["JE_K"] = m.group(1).replace("m", "-").replace("p", ".")
    m = re.search(r"_dt([0-9mp]+)", tag)
    if m:
        out["dt"] = m.group(1).replace("m", "-").replace("p", ".")
    m = re.search(r"_td([0-9mp]+)", tag)
    if m:
        out["t_damp"] = m.group(1).replace("m", "-").replace("p", ".")
    return out


def _has_numeric_rows(path: str) -> bool:
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


def main() -> None:
    ap = argparse.ArgumentParser(description="Collect summary from existing bench_ES runs.")
    ap.add_argument("--output-root", default="bench_ES", help="Folder under examples/tFreeEnergy_multi")
    args = ap.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    root = os.path.join(script_dir, args.output_root)
    runs_root = os.path.join(root, "runs")
    summary_root = os.path.join(root, "summary")
    os.makedirs(summary_root, exist_ok=True)

    rows: List[Dict[str, object]] = []
    for run_dir in sorted(glob.glob(os.path.join(runs_root, "N*", "*", "*"))):
        if not os.path.isdir(run_dir):
            continue
        parts = run_dir.split(os.sep)
        if len(parts) < 3:
            continue
        n_part = parts[-3]
        mode = parts[-2]
        tag = parts[-1]
        if not n_part.startswith("N"):
            continue
        try:
            n_atoms = int(n_part[1:])
        except Exception:
            continue

        row: Dict[str, object] = {
            "status": "failed",
            "N": n_atoms,
            "mode": mode,
            "run_dir": run_dir,
            "run_log": os.path.join(run_dir, "run.log"),
            "data_file": os.path.join(run_dir, f"entropic_spring_{n_atoms}_free_energy.dat"),
            "plot_file": os.path.join(run_dir, f"entropic_spring_{n_atoms}_free_energy_F_interactive.html"),
        }
        row.update(_parse_combo(tag))
        try:
            nl = float(row.get("nLambda", "nan"))
            nmd = float(row.get("nMDsteps", "nan"))
            if math.isfinite(nl) and math.isfinite(nmd) and nl > 0.0 and mode in ("JE", "BOTH"):
                row["je_effective_trajectories"] = nmd / nl
        except Exception:
            pass

        mpath = os.path.join(run_dir, "metrics.json")
        if os.path.isfile(mpath):
            try:
                with open(mpath, "r", encoding="utf-8") as f:
                    m = json.load(f)
                row.update(m)
                row["status"] = m.get("status", row["status"])
            except Exception:
                pass
        else:
            if os.path.isfile(row["data_file"]) and _has_numeric_rows(row["data_file"]):
                row["status"] = "ok"
            elif os.path.isfile(row["data_file"]):
                row["status"] = "failed_empty_data"

        rows.append(row)

    fieldnames = [
        "status", "N", "mode", "param_set", "nSys", "nLambda", "nEQsteps", "nMDsteps", "JE_K",
        "dt", "t_damp",
        "je_target_trajectories", "je_effective_trajectories",
        "wall_s", "md_ocl_time_ms_total", "md_ocl_calls", "md_ocl_avg_us_per_step",
        "md_ocl_time_ms_equilibrating", "md_ocl_time_ms_production", "md_ocl_time_ms_pulling", "md_ocl_time_ms_unknown",
        "rmsd_profile_eV", "final_abs_error_eV", "final_value_eV", "final_ref_eV",
        "je_wiggle_slope_eV_per_lambda", "je_wiggle_r2", "je_wiggle_rmse_eV",
        "run_dir", "run_log", "data_file", "plot_file",
    ]
    _write_csv(os.path.join(summary_root, "all_runs.csv"), rows, fieldnames)

    ok_rows = [r for r in rows if r.get("status") == "ok" and math.isfinite(float(r.get("rmsd_profile_eV", float("nan"))))]
    best_rows: List[Dict[str, object]] = []
    for n in sorted({int(r["N"]) for r in rows if str(r.get("N", "")).isdigit()}):
        for mode in ("TI", "JE", "BOTH"):
            subset = [r for r in ok_rows if int(r.get("N", -1)) == n and r.get("mode") == mode]
            if not subset:
                continue
            subset.sort(key=lambda r: (float(r["rmsd_profile_eV"]), float(r.get("wall_s", float("inf")))))
            best_rows.append(subset[0])
    _write_csv(os.path.join(summary_root, "best_by_method.csv"), best_rows, fieldnames)

    pareto_rows: List[Dict[str, object]] = []
    for n in sorted({int(r["N"]) for r in rows if str(r.get("N", "")).isdigit()}):
        for mode in ("TI", "JE", "BOTH"):
            subset = [r for r in ok_rows if int(r.get("N", -1)) == n and r.get("mode") == mode]
            pts = []
            for r in subset:
                w = float(r.get("wall_s", float("nan")))
                rr = float(r.get("rmsd_profile_eV", float("nan")))
                if math.isfinite(w) and math.isfinite(rr):
                    pts.append({"wall_s": w, "rmsd": rr, "row": r})
            for a in pts:
                if not any(_is_dominated(a, b) for b in pts if b is not a):
                    pareto_rows.append(a["row"])
    _write_csv(os.path.join(summary_root, "best_pareto.csv"), pareto_rows, fieldnames)

    try:
        import pandas as pd
        pd.DataFrame(rows).to_parquet(os.path.join(summary_root, "all_runs.parquet"), index=False)
    except Exception:
        pass

    _write_summary_plots(os.path.join(summary_root, "plots"), rows)
    print(f"Collected {len(rows)} runs into {summary_root}")


if __name__ == "__main__":
    main()
