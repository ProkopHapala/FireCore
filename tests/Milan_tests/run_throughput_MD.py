import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import MMFF_multi as mmff

import argparse

parser = argparse.ArgumentParser(description="Run throughput MD simulation with configurable parameters.")
parser.add_argument("--nSys", type=int, default=6000, help="Number of systems")
parser.add_argument("--xyz_name", type=str, default="data/xyz/xylitol_for_gridFF", help="Path to xyz file")
parser.add_argument("--surf_name", type=str, default=None, help="Path to surface file")
parser.add_argument("--bNonBonded", type=int, default=-1, help="bNonBonded flag")
parser.add_argument("--doSurfAtoms", type=int, default=0, help="doSurfAtoms flag")
parser.add_argument("--bSaveToDatabase", type=int, default=1, help="bSaveToDatabase flag")
parser.add_argument("--GridFF", type=int, default=1, help="GridFF flag")
parser.add_argument("--Fconv", type=float, default=1e-5, help="Force convergence")
parser.add_argument("--perframe", type=int, default=100, help="Per frame")
parser.add_argument("--perVF", type=int, default=100, help="Per VF")
parser.add_argument("--gridnPBC", type=str, default="(1,1,0)", help="gridnPBC")
parser.add_argument("--elapse_time", type=float, default=0.0, help="Elapse time")
args = parser.parse_args()

nSys = args.nSys
xyz_name = args.xyz_name
surf_name = args.surf_name
GridFF = args.GridFF
bNonBonded = args.bNonBonded
bSaveToDatabase = args.bSaveToDatabase
Fconv = args.Fconv
perframe = args.perframe
perVF = args.perVF
gridnPBC = tuple(map(int, args.gridnPBC.strip('()').split(',')))
elapse_time = args.elapse_time




mmff.setVerbosity(0)
print("Surface name:")
if(surf_name is None):
    print("No surface file specified")
else:
    print("Surface file: ", surf_name)

mmff.setSwitches( bSaveToDatabase=bSaveToDatabase )
mmff.setSwitches2( NonBonded=bNonBonded, SurfAtoms=0, GridFF=GridFF )
mmff.init(nSys_=nSys, xyz_name=xyz_name, surf_name=surf_name, T=300, gamma=0.1, nExplore=1000, nRelax=100000, pos_kick=0.25, vel_kick=1.0, GridFF=GridFF, gridnPBC=gridnPBC)
for i in range(10000):
    mmff.MDloop( perframe=perframe, Ftol=Fconv, iParalel=3, perVF=perVF, elapse_time=elapse_time )