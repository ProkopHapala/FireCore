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


def main() -> None:
    parser = argparse.ArgumentParser(description="Run multi-replica UFF scan with per-replica offsets")
    parser.add_argument("--xyz_name", type=str, default="data/xyz/xylitol.xyz", help="Path to xyz file (without .xyz suffix)")
    parser.add_argument("--surf_name", type=str, default="data/xyz/NaCl_1x1_L3", help="Path to surface file (without .xyz suffix)")
    parser.add_argument("--nSys", type=int, default=4, help="Number of replicas")
    parser.add_argument("--offset_start", type=str, default="(0,0,0)", help="Displacement applied to replica 0 [Ang]")
    parser.add_argument("--offset_end", type=str, default="(0,0,2)", help="Displacement applied to replica nSys-1 [Ang]")
    parser.add_argument("--bGridFF", type=int, default=-1, help="GridFF flag")
    parser.add_argument("--gridnPBC", type=str, default="(1,1,0)", help="GridFF periodicity")
    parser.add_argument("--bUFF", type=int, default=1, help="Enable UFF kernels")
    parser.add_argument("--bNonBonded", type=int, default=1, help="Enable non-bonded interactions")
    parser.add_argument("--doSurfAtoms", type=int, default=0, help="Enable surface atoms task")
    parser.add_argument("--dt", type=float, default=0.05, help="FIRE time step [ps]")
    parser.add_argument("--Fconv", type=float, default=1e-4, help="Force convergence threshold [eV/A]")
    parser.add_argument("--perframe", type=int, default=50, help="Iterations per MDloop call")
    parser.add_argument("--perVF", type=int, default=50, help="FIRE evaluations inside kernels")
    parser.add_argument("--loops", type=int, default=40, help="Number of MDloop calls")
    parser.add_argument("--T", type=float, default=1500.0, help="Exploration temperature [K]")
    parser.add_argument("--gamma", type=float, default=0.1, help="Langevin damping [1/ps]")
    parser.add_argument("--nExplore", type=int, default=0, help="Exploration length in MD steps")
    parser.add_argument("--nRelax", type=int, default=100000, help="Relaxation length in MD steps")
    parser.add_argument("--traj_base", type=str, default="traj_scan", help="Base name for trajectory output")
    parser.add_argument("--traj_stride", type=int, default=500, help="Steps between trajectory frames (negative disables)")
    args = parser.parse_args()

    offset_start = _parse_vec3(args.offset_start)
    offset_end = _parse_vec3(args.offset_end)
    offsets = _build_offsets(args.nSys, offset_start, offset_end)

    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    UFF_DIR = os.path.join(BASE_DIR, "common_resources")
    sElementTypes = os.path.join(UFF_DIR, "ElementTypes.dat")
    sAtomTypes = os.path.join(UFF_DIR, "AtomTypes.dat")
    sBondTypes = os.path.join(UFF_DIR, "BondTypes.dat")
    sAngleTypes = os.path.join(UFF_DIR, "AngleTypes.dat")
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
    atoms = uff.gpu_atoms
    natoms = uff.natoms
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


if __name__ == "__main__":
    main()
