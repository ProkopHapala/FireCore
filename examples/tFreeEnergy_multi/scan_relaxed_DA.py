import argparse
import os
import sys
from pathlib import Path

import numpy as np

sys.path.append("../../")
from pyBall import MMFF_multi as mm


"""
DA scan comparing relaxed and rigid branches of scan_Milan().
"""


ROOT = Path(__file__).resolve().parent
XYZ_IN = (ROOT / "../../cpp/common_resources/xyz/DA.xyz").resolve()
CONSTRAINTS_IN = (ROOT / "constraints_DA.txt").resolve()
TMP_XYZ = ROOT / "DA_scan_relaxed_atom0_Si.xyz"
OUT_RELAXED_XYZ = ROOT / "DA_scan_relaxed_via_scan.xyz"
OUT_RELAXED_DAT = ROOT / "DA_scan_relaxed_via_scan.dat"
OUT_RIGID_XYZ = ROOT / "DA_scan_rigid_via_scan.xyz"
OUT_RIGID_DAT = ROOT / "DA_scan_rigid_via_scan.dat"
OUT_COMPARE_DAT = ROOT / "DA_scan_compare_via_scan.dat"
OUT_HTML = ROOT / "DA_scan_compare_via_scan.html"


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
                "raw": line,
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


def save_scan_xyz(path: Path, elems, positions, energies, lambdas, target_pos, free_si_idx):
    with open(path, "w") as f:
        for il in range(len(lambdas)):
            f.write(f"{len(elems)}\n")
            free_si = positions[il, free_si_idx]
            comment = (
                f"lambda={lambdas[il]:.6f} "
                f"E={energies[il]:.8f} "
                f"targetSi=({target_pos[il,0]:.6f},{target_pos[il,1]:.6f},{target_pos[il,2]:.6f}) "
                f"freeSi=({free_si[0]:.6f},{free_si[1]:.6f},{free_si[2]:.6f})"
            )
            f.write(comment + "\n")
            for elem, p in zip(elems, positions[il]):
                f.write(f"{elem:<2s} {p[0]:16.8f} {p[1]:16.8f} {p[2]:16.8f}\n")


def run_scan(mode_name, b_relaxed, nCVs, initial_positions, final_positions, nLambda, niter, Fconv, K, natoms, si_indices, other_si_idx, lambdas, target_pos, target_dist):
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
    pulled_pos = ppos[:, si_indices[0], :]
    other_pos = ppos[:, other_si_idx, :]
    distances = np.linalg.norm(pulled_pos - other_pos, axis=1)
    dist_err = np.abs(distances - target_dist)
    return {
        "name": mode_name,
        "ppos": ppos,
        "Es": Es,
        "pulled_pos": pulled_pos,
        "other_pos": other_pos,
        "distances": distances,
        "dist_err": dist_err,
        "target_pos": target_pos,
        "lambdas": lambdas,
    }


def main():
    parser = argparse.ArgumentParser(description="DA scan comparing relaxed and rigid branches of scan_Milan()")
    parser.add_argument("--nSys", type=int, default=8, help="Number of parallel replicas in MMFFmulti")
    parser.add_argument("--nLambda", type=int, default=31, help="Number of scan points")
    parser.add_argument("--niter", type=int, default=20000, help="Max FIRE iterations per point")
    parser.add_argument("--dt", type=float, default=0.05, help="FIRE dt_max")
    parser.add_argument("--Fconv", type=float, default=1e-4, help="Force convergence")
    parser.add_argument("--K", type=float, default=10.0, help="Constraint spring constant")
    args = parser.parse_args()

    if not XYZ_IN.exists():
        raise FileNotFoundError(f"Missing input XYZ: {XYZ_IN}")
    if not CONSTRAINTS_IN.exists():
        raise FileNotFoundError(f"Missing constraints file: {CONSTRAINTS_IN}")

    natoms, comment, atoms = read_xyz(XYZ_IN)
    constraints = load_constraints(CONSTRAINTS_IN)
    if len(constraints) == 0:
        raise RuntimeError(f"No constraints loaded from {CONSTRAINTS_IN}")
    target_start = np.array(constraints[0][0], dtype=np.float64)
    target_end   = np.array(constraints[0][1], dtype=np.float64)
    write_xyz(TMP_XYZ, f"{comment} | scan_relaxed endpoint mode", atoms)

    print(f"Input XYZ       : {XYZ_IN}")
    print(f"Temporary XYZ   : {TMP_XYZ}")
    print(f"Constraints     : {CONSTRAINTS_IN}")
    print(f"Target start    : {target_start}")
    print(f"Target end      : {target_end}")
    print("")
    print("Note: this uses the same endpoint-style constraints as computeFreeEnergy().")

    mm.init(
        nSys_=args.nSys,
        xyz_name=str(TMP_XYZ),
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
    mm.getBuffs()
    mm.set_opt(
        dt_max=args.dt,
        dt_min=max(args.dt * 0.1, 1.0e-6),
        damp_max=0.2,
        finc=1.1,
        fdec=0.5,
        falpha=0.8,
        minLastNeg=5,
        cvf_min=-0.1,
        cvf_max=0.1,
    )

    si_indices = [ia for ia in range(mm.natoms) if mm.getTypeName(ia).strip() == "Si"]
    if len(si_indices) < 2:
        raise RuntimeError(f"Need at least 2 Si atoms after reorder, found {len(si_indices)}")
    other_si_idx = si_indices[1]
    ref_positions = mm.apos.copy()
    ref_pulled = ref_positions[si_indices[0]].copy()
    ref_other = ref_positions[other_si_idx].copy()
    ref_dist = np.linalg.norm(ref_pulled - ref_other)

    print(f"First Si index               : {si_indices[0]}")
    print(f"Second Si index              : {other_si_idx}")
    print(f"Reference pulled Si position : {ref_pulled}")
    print(f"Reference second Si position : {ref_other}")
    print(f"Reference Si-Si distance     : {ref_dist:.6f} A")

    lambdas = np.linspace(0.0, 1.0, args.nLambda)
    target_pos = (1.0 - lambdas[:, None]) * target_start[None, :] + lambdas[:, None] * target_end[None, :]
    target_dist = np.linalg.norm(
        ((1.0 - lambdas[:, None]) * target_start[None, :] + lambdas[:, None] * target_end[None, :]) -
        ((1.0 - lambdas[:, None]) * np.array(constraints[1][0], dtype=np.float64)[None, :] +
         lambdas[:, None] * np.array(constraints[1][1], dtype=np.float64)[None, :]),
        axis=1,
    )
    if (len(constraints) % 2) != 0:
        raise RuntimeError(f"Expected an even number of constraint rows, got {len(constraints)}")
    initial_positions = []
    final_positions = []
    for init_pos, final_pos in constraints:
        initial_positions.extend(init_pos)
        final_positions.extend(final_pos)
    nCVs = len(constraints) // 2

    relaxed = run_scan(
        "relaxed",
        True,
        nCVs,
        initial_positions,
        final_positions,
        args.nLambda,
        args.niter,
        args.Fconv,
        args.K,
        mm.natoms,
        si_indices,
        other_si_idx,
        lambdas,
        target_pos,
        target_dist,
    )
    rigid = run_scan(
        "rigid",
        False,
        nCVs,
        initial_positions,
        final_positions,
        args.nLambda,
        args.niter,
        args.Fconv,
        args.K,
        mm.natoms,
        si_indices,
        other_si_idx,
        lambdas,
        target_pos,
        target_dist,
    )

    elems = [mm.getTypeName(i).strip() for i in range(mm.natoms)]

    relaxed_table = np.column_stack(
        [
            np.arange(args.nLambda, dtype=np.float64),
            lambdas,
            target_dist,
            relaxed["distances"],
            relaxed["Es"],
            relaxed["dist_err"],
            relaxed["pulled_pos"],
            relaxed["other_pos"],
        ]
    )
    rigid_table = np.column_stack(
        [
            np.arange(args.nLambda, dtype=np.float64),
            lambdas,
            target_dist,
            rigid["distances"],
            rigid["Es"],
            rigid["dist_err"],
            rigid["pulled_pos"],
            rigid["other_pos"],
        ]
    )
    compare_table = np.column_stack(
        [
            np.arange(args.nLambda, dtype=np.float64),
            lambdas,
            target_dist,
            relaxed["distances"],
            relaxed["Es"],
            relaxed["dist_err"],
            rigid["distances"],
            rigid["Es"],
            rigid["dist_err"],
        ]
    )

    np.savetxt(
        OUT_RELAXED_DAT,
        relaxed_table,
        header="il lambda target_dist dist_sisi E dist_err pulled_x pulled_y pulled_z other_x other_y other_z",
    )
    np.savetxt(
        OUT_RIGID_DAT,
        rigid_table,
        header="il lambda target_dist dist_sisi E dist_err pulled_x pulled_y pulled_z other_x other_y other_z",
    )
    np.savetxt(
        OUT_COMPARE_DAT,
        compare_table,
        header="il lambda target_dist relaxed_dist relaxed_E relaxed_dist_err rigid_dist rigid_E rigid_dist_err",
    )
    save_scan_xyz(OUT_RELAXED_XYZ, elems, relaxed["ppos"], relaxed["Es"], lambdas, target_pos, other_si_idx)
    save_scan_xyz(OUT_RIGID_XYZ, elems, rigid["ppos"], rigid["Es"], lambdas, target_pos, other_si_idx)

    print("")
    print(f"Relaxed table : {OUT_RELAXED_DAT}")
    print(f"Relaxed XYZ   : {OUT_RELAXED_XYZ}")
    print(f"Rigid table   : {OUT_RIGID_DAT}")
    print(f"Rigid XYZ     : {OUT_RIGID_XYZ}")
    print(f"Compare table : {OUT_COMPARE_DAT}")
    print("")
    print(f"{'il':>3s} {'lambda':>8s} {'target[A]':>12s} {'dist_rel[A]':>12s} {'E_rel[eV]':>14s} {'dist_rig[A]':>12s} {'E_rig[eV]':>14s}")
    for il in range(args.nLambda):
        print(
            f"{il:3d} {lambdas[il]:8.4f} {target_dist[il]:12.6f} "
            f"{relaxed['distances'][il]:12.6f} {relaxed['Es'][il]:14.6f} "
            f"{rigid['distances'][il]:12.6f} {rigid['Es'][il]:14.6f}"
        )

    try:
        from plotly.subplots import make_subplots
        import plotly.graph_objects as go

        relaxed_dE = relaxed["Es"] - np.nanmin(relaxed["Es"])
        rigid_dE = rigid["Es"] - np.nanmin(rigid["Es"])

        custom_rel = np.column_stack(
            [
                lambdas,
                target_dist,
                relaxed["dist_err"],
                relaxed["pulled_pos"],
                relaxed["other_pos"],
            ]
        )
        custom_rig = np.column_stack(
            [
                lambdas,
                target_dist,
                rigid["dist_err"],
                rigid["pulled_pos"],
                rigid["other_pos"],
            ]
        )

        fig = make_subplots(
            rows=1,
            cols=3,
            subplot_titles=(
                "Energy profile comparison",
                "Distance error",
                "Energy vs lambda",
            ),
            horizontal_spacing=0.08,
        )

        hover_tmpl = (
            "lambda=%{customdata[0]:.4f}<br>"
            "target d(Si-Si)=%{customdata[1]:.6f} A<br>"
            "dist err=%{customdata[2]:.6f} A<br>"
            "pulled Si=(%{customdata[3]:.6f}, %{customdata[4]:.6f}, %{customdata[5]:.6f})<br>"
            "other Si=(%{customdata[6]:.6f}, %{customdata[7]:.6f}, %{customdata[8]:.6f})<extra>%{fullData.name}</extra>"
        )

        fig.add_trace(
            go.Scatter(
                x=relaxed["distances"],
                y=relaxed_dE,
                mode="lines+markers",
                name="relaxed",
                line=dict(width=2),
                marker=dict(symbol="circle", size=7),
                customdata=custom_rel,
                hovertemplate=hover_tmpl,
            ),
            row=1, col=1
        )
        fig.add_trace(
            go.Scatter(
                x=rigid["distances"],
                y=rigid_dE,
                mode="lines+markers",
                name="rigid",
                line=dict(width=2, dash="dash"),
                marker=dict(symbol="square", size=7),
                customdata=custom_rig,
                hovertemplate=hover_tmpl,
            ),
            row=1, col=1
        )

        fig.add_trace(
            go.Scatter(
                x=lambdas,
                y=relaxed["dist_err"],
                mode="lines+markers",
                name="relaxed",
                line=dict(width=2),
                marker=dict(symbol="circle", size=7),
                customdata=custom_rel,
                hovertemplate=hover_tmpl,
                showlegend=False,
            ),
            row=1, col=2
        )
        fig.add_trace(
            go.Scatter(
                x=lambdas,
                y=rigid["dist_err"],
                mode="lines+markers",
                name="rigid",
                line=dict(width=2, dash="dash"),
                marker=dict(symbol="square", size=7),
                customdata=custom_rig,
                hovertemplate=hover_tmpl,
                showlegend=False,
            ),
            row=1, col=2
        )

        fig.add_trace(
            go.Scatter(
                x=lambdas,
                y=relaxed_dE,
                mode="lines+markers",
                name="relaxed",
                line=dict(width=2),
                marker=dict(symbol="circle", size=7),
                customdata=custom_rel,
                hovertemplate=hover_tmpl,
                showlegend=False,
            ),
            row=1, col=3
        )
        fig.add_trace(
            go.Scatter(
                x=lambdas,
                y=rigid_dE,
                mode="lines+markers",
                name="rigid",
                line=dict(width=2, dash="dash"),
                marker=dict(symbol="square", size=7),
                customdata=custom_rig,
                hovertemplate=hover_tmpl,
                showlegend=False,
            ),
            row=1, col=3
        )

        fig.update_xaxes(title_text="Si-Si distance [A]", row=1, col=1)
        fig.update_yaxes(title_text="Delta E [eV]", row=1, col=1)

        fig.update_xaxes(title_text="lambda", row=1, col=2)
        fig.update_yaxes(title_text="|Si-Si distance - target| [A]", row=1, col=2)

        fig.update_xaxes(title_text="lambda", row=1, col=3)
        fig.update_yaxes(title_text="Delta E [eV]", row=1, col=3)

        fig.update_layout(
            title="DA scan: relaxed vs rigid",
            template="plotly_white",
            width=1600,
            height=520,
            legend=dict(x=0.5, y=1.12, xanchor="center", orientation="h"),
            margin=dict(l=60, r=30, t=80, b=60),
        )

        fig.write_html(str(OUT_HTML), include_plotlyjs="cdn", full_html=True)
        print(f"Plot          : {OUT_HTML}")

    except Exception as exc:
        print(f"Plot skipped: {exc}")


if __name__ == "__main__":
    main()