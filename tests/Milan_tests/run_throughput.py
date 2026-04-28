import sys
import numpy as np
import argparse
import os
import glob

sys.path.append("../../")
from pyBall import MMFF_multi as ff

# ========================================================================
# ====================== Argument Parsing ================================
# ========================================================================

parser = argparse.ArgumentParser(description="Run UFF/MMFF throughput MD with configurable parameters.")

# --- Force Field Selection ---
parser.add_argument("--ff", type=str, default="mmff", help="Force field: 'mmff' or 'uff'")

# --- Common arguments ---
parser.add_argument("--xyz_name", type=str, help="Path to xyz file (without .xyz suffix)")
parser.add_argument("--surf_name", type=str, default=None, help="Path to surface file (without .xyz suffix)")
parser.add_argument("--nSys", type=int, default=1, help="Number of systems")
parser.add_argument("--Fconv", type=float, default=1e-4, help="Force convergence")
parser.add_argument("--perframe", type=int, default=10, help="Steps per frame (MDloop nIter)")
parser.add_argument("--perVF", type=int, default=10, help="Vector-field evals inside kernels")
parser.add_argument("--gridnPBC", type=str, default="(1,1,0)", help="gridnPBC")
parser.add_argument("--doSurfAtoms", type=int, default=0, help="doSurfAtoms flag")
parser.add_argument("--GridFF", type=int, default=-1, help="GridFF flag (bGridFF for UFF)")
parser.add_argument("--loops", type=int, default=50000, help="How many times to call MDloop in a row")
parser.add_argument("--elapse_time", type=float, default=0.0, help="Elapse time")
parser.add_argument("--T", type=float, default=300.0, help="Thermostat target temperature [K]")
parser.add_argument("--gamma", type=float, default=0.1, help="Langevin damping [1/ps, in internal units]")
parser.add_argument("--nExplore", type=int, default=1000, help="Exploring duration in MD steps")
parser.add_argument("--nRelax", type=int, default=100000, help="Relaxation duration in MD steps")
parser.add_argument("--bSaveToDatabase", type=int, default=1, help="bSaveToDatabase flag for MMFF")
parser.add_argument("--bUFF", type=int, default=1, help="bUFF flag for UFF")
parser.add_argument("--dt", type=float, default=0.05, help="Time step for FIRE optimizer (UFF)")
parser.add_argument("--bNonBonded", type=int, default=1, help="bNonBonded flag for UFF")

args = parser.parse_args()

# ========================================================================
# ====================== Initialization ==================================
# ========================================================================

ff.setVerbosity(0)

gridnPBC = tuple(map(int, args.gridnPBC.strip('()').split(',')))

sElementTypes, sAtomTypes, sBondTypes, sAngleTypes, sDihedralTypes = None, None, None, None, None
bMMFF = True
bUFF = False

if args.ff == 'uff':
    bUFF = True
    print("Initializing UFF")
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    UFF_DIR = os.path.join(BASE_DIR, "common_resources")
    sElementTypes = os.path.join(UFF_DIR, "ElementTypes.dat")
    sAtomTypes = os.path.join(UFF_DIR, "AtomTypes.dat")
    sBondTypes = os.path.join(UFF_DIR, "BondTypes.dat")
    sAngleTypes = os.path.join(UFF_DIR, "AngleTypes.dat")
    sDihedralTypes = os.path.join(UFF_DIR, "DihedralTypes.dat")
elif args.ff == 'mmff':
    print("Initializing MMFF")
else:
    print(f"Error: Unknown force field '{args.ff}'")
    sys.exit(1)

# these should be set for both
ff.setSwitches(bSaveToDatabase=args.bSaveToDatabase)
ff.setSwitches2(
    NonBonded=args.bNonBonded,
    SurfAtoms=args.doSurfAtoms,
    GridFF=args.GridFF,
)


ff.init(
    nSys_=args.nSys,
    xyz_name=args.xyz_name,
    surf_name=args.surf_name,
    # sElementTypes=sElementTypes,
    # sAtomTypes=sAtomTypes,
    # sBondTypes=sBondTypes,
    # sAngleTypes=sAngleTypes,
    # sDihedralTypes=sDihedralTypes,
    bMMFF=bMMFF,
    bUFF=bUFF,
    GridFF=args.GridFF,
    gridnPBC=gridnPBC,
    T=args.T,
    gamma=args.gamma,
    nExplore=args.nExplore,
    nRelax=args.nRelax
)

if args.ff == 'uff':
    ff.set_dt_default(args.dt)

    if args.bUFF:
        ff.setSwitchesUFF(
            DoBond=1, DoAngle=1, DoDihedral=1, DoInversion=1, DoAssemble=1,
            SubtractBondNonBond=-1 if args.bNonBonded else -1,
            ClampNonBonded=-1 if args.bNonBonded else -1
        )

# Remove old trajectories
for _f in glob.glob("traj_UFF_*.xyz"):
    try:
        os.remove(_f)
    except FileNotFoundError:
        pass

# Set trajectory output (per-system files: traj_UFF_000.xyz, ...)
ff.setTrjName(trj_fname_="traj_UFF", savePerNsteps=-10, bDel=True)

# ========================================================================
# ====================== MD Loop =========================================
# ========================================================================

print(f"Starting MD loop for {args.ff}")
iParalel = 3
for i in range(args.loops):
    ff.MDloop(perframe=args.perframe, Ftol=args.Fconv, iParalel=iParalel, perVF=args.perVF, elapse_time=args.elapse_time)
    
    if (i + 1) % 50 == 0:
        print(f"[{args.ff.upper()}] MDloop {i+1}/{args.loops} done")

print("MD loops finished.")