"""

Key files: 
- `run_throughput_UFF.py` : driver
- `pyBall/MMFF_multi.py` : ctypes interface to C/C++ library 
- `cpp/libs_OCL/MMFFmulti_lib.cpp` : C/C++ library
- `cpp/common/molecular/MolWorld_sp3_multi.h` : molecular world used in C/C++ library 
- `cpp/common/OpenCL/OCL_UFF.h` : OpenCL driver used in MolWorld_sp3_multi.h
- `cpp/common_resources/cl/UFF.cl` : kernels used in OCL_UFF.h

GridFF single-system run: 
python run_throughput_UFF.py --nSys 1 --xyz_name data_UFF/xyz/HHO-h.xyz --bGridFF 1 --gridnPBC "(1,1,0)" --loops 10 --perframe 500 --perVF 100 --Fconv 1e-4

"""

import sys
import numpy as np
import argparse

sys.path.append("../../")
from pyBall import MMFF_multi as uff

# Exploring (minima hopping) controls
parser = argparse.ArgumentParser(description="Run UFF throughput MD with configurable parameters.")
parser.add_argument("--xyz_name",        type=str,   default="data/xyz/xylitol.xyz", help="Path to xyz file (without .xyz suffix)")
parser.add_argument("--surf_name",       type=str,   default="data/xyz/NaCl_1x1_L3",            help="Path to surface file (without .xyz suffix)")
parser.add_argument("--nSys",            type=int,   default=1,         help="Number of systems")
parser.add_argument("--dovdW",           type=int,   default=1,         help="dovdW flag")
parser.add_argument("--doSurfAtoms",     type=int,   default=0,         help="doSurfAtoms flag")
parser.add_argument("--bSaveToDatabase", type=int,   default=-1,        help="bSaveToDatabase flag") # DO NOT CHANGE (but has no effect)
parser.add_argument("--bGridFF",         type=int,   default=-1,        help="bGridFF flag")
parser.add_argument("--bUFF",            type=int,   default=1,         help="bUFF flag")
parser.add_argument("--Fconv",           type=float, default=1e-4,      help="Force convergence")
parser.add_argument("--dt",              type=float, default=0.05,     help="Time step for FIRE optimizer (reduce for light molecules like H2O)")
parser.add_argument("--perframe",        type=int,   default=10,      help="Steps per frame (MDloop nIter)")
parser.add_argument("--perVF",           type=int,   default=10,       help="Vector-field evals inside kernels")
parser.add_argument("--loops",           type=int,   default=5,       help="How many times to call MDloop in a row") #should be set to large number, the duration is set internaly in Molworld_sp3_multi::MDLoop function
parser.add_argument("--gridnPBC",        type=str,   default="(1,1,0)", help="gridnPBC")
parser.add_argument("--bNonBonded",      type=int,   default=1,         help="bNonBonded flag")
parser.add_argument("--T",               type=float, default=1500.0,    help="Thermostat target temperature [K] during exploring")
parser.add_argument("--gamma",           type=float, default=0.1,       help="Langevin damping [1/ps, in internal units]")
parser.add_argument("--nExplore",        type=int,   default=700,       help="Exploring duration in MD steps (perVF steps accumulate)")
parser.add_argument("--nRelax",          type=int,   default=100000,    help="Relaxation duration in MD steps before forcing exploring if not converged")
args = parser.parse_args()

nSys = args.nSys
xyz_name = args.xyz_name
surf_name = args.surf_name
Fconv = args.Fconv
perframe = args.perframe
perVF = args.perVF
loops = args.loops
bNonBonded = args.bNonBonded
bGridFF = args.bGridFF
gridnPBC = tuple(map(int, args.gridnPBC.strip('()').split(',')))
T = args.T
gamma = args.gamma
nExplore = args.nExplore
nRelax = args.nRelax
doSurfAtoms = args.doSurfAtoms

# Reduce verbosity for speed
uff.setVerbosity(0)

# Use UFF types from data_UFF (relative to this script)
import os
import glob

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
UFF_DIR = os.path.join(BASE_DIR, "common_resources")
sElementTypes = os.path.join(UFF_DIR, "ElementTypes.dat")
sAtomTypes    = os.path.join(UFF_DIR, "AtomTypes.dat")
sBondTypes    = os.path.join(UFF_DIR, "BondTypes.dat")
sAngleTypes   = os.path.join(UFF_DIR, "AngleTypes.dat")
sDihedralTypes= os.path.join(UFF_DIR, "DihedralTypes.dat")

# Initialize for UFF (bUFF=True); GridFF off for UFF
uff.init(
    nSys_=nSys,
    xyz_name=xyz_name,
    surf_name=surf_name,
    sElementTypes=sElementTypes,
    sAtomTypes=sAtomTypes,
    sBondTypes=sBondTypes,
    sAngleTypes=sAngleTypes,
    sDihedralTypes=sDihedralTypes,
    bMMFF=True, # this must be true for both UFF and MMFF
    bUFF=args.bUFF,  # if this is false, we use MMFF
    GridFF=bGridFF, gridnPBC=gridnPBC,
    T=T, gamma=gamma, nExplore=nExplore, nRelax=nRelax
)
# Set time step (reduce for light molecules like H2O)
uff.set_dt_default(args.dt)

# First set core switches (enables NonBonded), then UFF component/clamp flags
uff.setSwitches2(
        NonBonded=bNonBonded,
        SurfAtoms=doSurfAtoms,
        GridFF=bGridFF,
    )
if args.bUFF:
    uff.setSwitchesUFF(
        DoBond=1, DoAngle=1, DoDihedral=1, DoInversion=1, DoAssemble=1,
        SubtractBondNonBond = -1 if bNonBonded else -1,
        ClampNonBonded      = -1 if bNonBonded else -1
    )


# Remove old trajectories
for _f in glob.glob("traj_UFF_*.xyz"):
    try:
        os.remove(_f)
    except FileNotFoundError:
        pass

# Set trajectory output (per-system files: traj_UFF_000.xyz, ...)
# set savePerNsteps to change the frequency of saving (negative to disable)
uff.setTrjName(trj_fname_="traj_UFF", savePerNsteps=1000, bDel=True)

# Throughput loop: call MDloop repeatedly
for i in range(loops):
    uff.MDloop(perframe=perframe, Ftol=Fconv, perVF=perVF)
    # Optionally print simple progress every so often
    if (i+1) % 50 == 0:
        print(f"[UFF] MDloop {i+1}/{loops} done")

