#!/usr/bin/env python3
"""Multi-replica UFF driver that distributes replicas along a 1D scan."""

import argparse
import glob
import os
import sys

import numpy as np

sys.path.append("../../")
from pyBall import MMFF_multi as uff


def _parse_vec3(text: str) -> np.ndarray:
    raw = text.strip()
    if raw.startswith("(") and raw.endswith(")"):
        raw = raw[1:-1]
    parts = [p.strip() for p in raw.split(",") if p.strip()]
    if len(parts) != 3:
        raise argparse.ArgumentTypeError(f"Expected three comma-separated values, got '{text}'")
    try:
        vec = np.array([float(p) for p in parts], dtype=np.float64)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"Cannot parse vector '{text}'") from exc
    return vec


def _build_offsets(n_sys: int, start: np.ndarray, end: np.ndarray) -> np.ndarray:
    if n_sys < 1:
        raise ValueError("nSys must be positive")
    if n_sys == 1:
        return start[None, :]
    lambdas = np.linspace(0.0, 1.0, n_sys, dtype=np.float64)
    direction = end - start
    return start[None, :] + lambdas[:, None] * direction[None, :]


def _scan_coordinates(offsets: np.ndarray, start: np.ndarray, end: np.ndarray) -> tuple[np.ndarray, float]:
    direction = end - start
    length = float(np.linalg.norm(direction))
    if length > 0.0:
        unit = direction / length
        coords = (offsets - start[None, :]) @ unit
    else:
        coords = np.zeros(offsets.shape[0], dtype=np.float64)
    return coords, length


def _save_energy_table(path: str, coords: np.ndarray, energies: np.ndarray) -> None:
    data = np.column_stack([np.arange(len(coords)), coords, energies])
    header = "index coordinate energy"
    np.savetxt(path, data, header=header)
    print(f"[scan] Energy table written to {path}")


def _plot_energy_profile(path: str, coords: np.ndarray, energies: np.ndarray) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    order = np.argsort(coords)
    xs = coords[order]
    ys = energies[order]

    fig, ax = plt.subplots(figsize=(7.0, 4.5))
    ax.plot(xs, ys, marker="o", linewidth=1.5, markersize=4, label="relaxed energy")
    ax.set_xlabel("scan coordinate [Å]")
    ax.set_ylabel("Energy [eV]")
    ax.set_title("Relaxed UFF energy along scan")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    print(f"[scan] Energy plot saved to {path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run multi-replica UFF scan with per-replica offsets")
    parser.add_argument("--xyz_name",     type=str,   default="data/xyz/xylitol.xyz", help="Path to xyz file (without .xyz suffix)")
    parser.add_argument("--surf_name",    type=str,   default="data/xyz/NaCl_1x1_L3", help="Path to surface file (without .xyz suffix)")
    parser.add_argument("--nSys",         type=int,   default=100,                      help="Number of replicas")
    parser.add_argument("--offset_start", type=str,   default="(0,0,0)",              help="Displacement applied to replica 0 [Ang]")
    parser.add_argument("--offset_end",   type=str,   default="(0,0,2)",              help="Displacement applied to replica nSys-1 [Ang]")
    parser.add_argument("--bGridFF",      type=int,   default=1,                     help="GridFF flag")
    parser.add_argument("--gridnPBC",     type=str,   default="(1,1,0)",              help="GridFF periodicity")
    parser.add_argument("--bUFF",         type=int,   default=1,                      help="Enable UFF kernels")
    parser.add_argument("--bNonBonded",   type=int,   default=1,                      help="Enable non-bonded interactions")
    parser.add_argument("--doSurfAtoms",  type=int,   default=1,                      help="Enable surface atoms task")
    parser.add_argument("--dt",           type=float, default=0.01,                   help="FIRE time step [ps]")
    parser.add_argument("--Fconv",        type=float, default=1e-4,                   help="Force convergence threshold [eV/A]")
    parser.add_argument("--perframe",     type=int,   default=1,                      help="Iterations per MDloop call")
    parser.add_argument("--perVF",        type=int,   default=1,                      help="FIRE evaluations inside kernels")
    parser.add_argument("--loops",        type=int,   default=1,                      help="Number of MDloop calls")
    parser.add_argument("--T",            type=float, default=1500.0,                 help="Exploration temperature [K]")
    parser.add_argument("--gamma",        type=float, default=0.1,                    help="Langevin damping [1/ps]")
    parser.add_argument("--nExplore",     type=int,   default=0,                      help="Exploration length in MD steps")
    parser.add_argument("--nRelax",       type=int,   default=1,                      help="Relaxation length in MD steps")
    parser.add_argument("--traj_base",    type=str,   default="traj_scan",            help="Base name for trajectory output")
    parser.add_argument("--traj_stride",  type=int,   default=500,                    help="Steps between trajectory frames (negative disables)")
    parser.add_argument("--energy_table", type=str,   default="scan_energy.csv",      help="Write relaxed energies to CSV (empty to skip)")
    parser.add_argument("--energy_plot",  type=str,   default="scan_energy.png",      help="Write relaxed energy plot (empty to skip)")
    args = parser.parse_args()

    offset_start = _parse_vec3(args.offset_start)
    offset_end   = _parse_vec3(args.offset_end)
    offsets      = _build_offsets(args.nSys, offset_start, offset_end)

    BASE_DIR       = os.path.dirname(os.path.abspath(__file__))
    UFF_DIR        = os.path.join(BASE_DIR, "common_resources")
    sElementTypes  = os.path.join(UFF_DIR, "ElementTypes.dat")
    sAtomTypes     = os.path.join(UFF_DIR, "AtomTypes.dat")
    sBondTypes     = os.path.join(UFF_DIR, "BondTypes.dat")
    sAngleTypes    = os.path.join(UFF_DIR, "AngleTypes.dat")
    sDihedralTypes = os.path.join(UFF_DIR, "DihedralTypes.dat")

    gridnPBC = tuple(map(int, args.gridnPBC.strip("()").split(",")))

    uff.setVerbosity(0)
    uff.init(
        nSys_=args.nSys,
        xyz_name=args.xyz_name,
        surf_name=args.surf_name,
        sElementTypes=sElementTypes,
        sAtomTypes=sAtomTypes,
        sBondTypes=sBondTypes,
        sAngleTypes=sAngleTypes,
        sDihedralTypes=sDihedralTypes,
        bMMFF=True,
        bUFF=bool(args.bUFF),
        GridFF=args.bGridFF,
        gridnPBC=gridnPBC,
        T=args.T,
        gamma=args.gamma,
        nExplore=args.nExplore,
        nRelax=args.nRelax,
    )

    uff.set_dt_default(args.dt)
    uff.setSwitches2(
        NonBonded=args.bNonBonded,
        SurfAtoms=args.doSurfAtoms,
        GridFF=args.bGridFF,
    )
    if args.bUFF:
        uff.setSwitchesUFF(
            DoBond=1,
            DoAngle=1,
            DoDihedral=1,
            DoInversion=1,
            DoAssemble=1,
            SubtractBondNonBond=-1 if args.bNonBonded else -1,
            ClampNonBonded=-1 if args.bNonBonded else -1,
        )

    uff.download(bForces=False, bVel=False)
    uff.getBuffs_UFF()
    natoms = uff.natoms
    gpu_atoms = getattr(uff, "gpu_atoms", None)
    if gpu_atoms is None:
        gpu_atoms = uff.getfBuff("gpu_atoms", (args.nSys, natoms, 4))
    atoms = gpu_atoms
    if atoms.shape[0] != args.nSys:
        raise RuntimeError(f"gpu_atoms first dimension {atoms.shape[0]} != nSys {args.nSys}")
    base_positions = atoms[0, :natoms, :3].copy()
    for isys in range(args.nSys):
        atoms[isys, :natoms, :3] = base_positions + offsets[isys]
    uff.upload(bParams=False, bForces=False, bVel=False, blvec=False)

    for path in glob.glob(f"{args.traj_base}_*.xyz"):
        try:
            os.remove(path)
        except FileNotFoundError:
            pass
    uff.setTrjName(trj_fname_=args.traj_base, savePerNsteps=args.traj_stride, bDel=False)

    for i in range(args.loops):
        uff.MDloop(perframe=args.perframe, Ftol=args.Fconv, perVF=args.perVF)
        if (i + 1) % 50 == 0:
            print(f"[scan] MDloop {i + 1}/{args.loops}")

    energies = np.zeros(args.nSys, dtype=np.double)
    total = uff.get_uff_energy(energies, isys_choice=-1, download=True)
    print(f"[scan] Total energy (summed over replicas) = {total: .6e} eV")

    scan_coords, scan_length = _scan_coordinates(offsets, offset_start, offset_end)
    print(f"[scan] Scan length = {scan_length: .6f} Å")
    for isys, (coord, energy) in enumerate(zip(scan_coords, energies)):
        print(f"  replica {isys:3d} : coord={coord: .6f} Å  energy={energy: .6e} eV")

    if args.energy_table:
        _save_energy_table(args.energy_table, scan_coords, energies)
    if args.energy_plot:
        _plot_energy_profile(args.energy_plot, scan_coords, energies)