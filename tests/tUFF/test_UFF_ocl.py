#!/usr/bin/env python3
"""Simple UFF validation test comparing CPU and GPU implementations"""

import sys
import os
import argparse

# Add FireCore to path
base_path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(base_path)

data_dir = os.path.join(base_path, "cpp/common_resources")

import numpy as np
from pyBall.AtomicSystem import AtomicSystem 
from pyBall import     MMFF as uff_cpp
from pyBall.OCL import UFF  as uff_ocl

def run_uff_cpp( args ):        
    uff_cpp.init(
        xyz_name=args.molecule,
        sElementTypes=os.path.join(data_dir, "ElementTypes.dat"),
        sAtomTypes=os.path.join(data_dir, "AtomTypes.dat"),
        sBondTypes=os.path.join(data_dir, "BondTypes.dat"),
        sAngleTypes=os.path.join(data_dir, "AngleTypes.dat"),
        sDihedralTypes=os.path.join(data_dir, "DihedralTypes.dat"),
        bMMFF=True,
        bUFF=True,
        b141=True,
    )
    uff_cpp.setSwitches(NonBonded=-1, SurfAtoms=-1, GridFF=-1)
    uff_cpp.setSwitchesUFF( DoAssemble=+1, DoBond=1, DoAngle=-1, DoDihedral=-1, DoInversion=-1,  SubtractBondNonBond=-1, ClampNonBonded=-1 )
    print("-----------\n uff_cpp.getBuffs_UFF()  ")
    uff_cpp.getBuffs_UFF()
    print("-----------\n uff_cpp.eval()  ")
    energy = uff_cpp.eval()
    forces = uff_cpp.fapos.copy()
    return energy, forces

def run_uff_ocl(args):
    mol      = AtomicSystem(fname=args.molecule)
    uff_cl   = uff_ocl.UFF_CL(nloc=32)
    uff_data = uff_cl.toUFF(mol)
    uff_cl.set_component_flags(
        bBonds=args.components['bonds'],
        bAngles=args.components['angles'],
        bDihedrals=args.components['dihedrals'],
        bInversions=args.components['inversions']
    )
    uff_cl.upload_topology_params(uff_data)
    uff_cl.upload_positions(mol.apos)
    uff_cl.run_eval_step(bClearForce=True)
    energy = uff_cl.get_total_energy()
    forces = uff_cl.get_forces()
    return energy, forces

def compare_results(cpu_energy, cpu_forces, gpu_energy, gpu_forces):
    print(f"CPU Energy: {cpu_energy}")
    print(f"GPU Energy: {gpu_energy}")
    print(f"Energy diff: {abs(cpu_energy - gpu_energy)}")
    
    force_diff = np.linalg.norm(cpu_forces - gpu_forces, axis=1)
    print(f"Max force diff: {np.max(force_diff)}")
    print(f"Mean force diff: {np.mean(force_diff)}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='UFF Validation Test Script')
    parser.add_argument('-m','--molecule',      type=str,     default='../../cpp/common_resources/mol/H2O.mol2',             help='Molecule file (.mol, .xyz)'    )
    parser.add_argument('-g','--gpu',           type=int,     default=0,                     help='Use GPU/OpenCL implementation' )
    parser.add_argument('-c','--components',    default='bonds,angles,dihedrals,inversions', help='Comma-separated list of UFF components to test')
    parser.add_argument('-v','--verbose',       action='count', default=1,                   help='Increase verbosity level'    )
    parser.add_argument('-t','--tolerance',     type=float,   default=1e-6,                  help='Energy comparison tolerance' )
    args = parser.parse_args()

    cpu_energy, cpu_forces = run_uff_cpp(args)
    
    if args.gpu:
        gpu_energy, gpu_forces = run_uff_ocl(args)
        compare_results(cpu_energy, cpu_forces, gpu_energy, gpu_forces)
    else:
        print(f"CPU Energy: {cpu_energy}")
        print(f"CPU Forces:\n{cpu_forces}")

    print("=========\nALL DONE")
