#!/usr/bin/python3
import sys
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import argparse

# Remove custom pyopencl from path to use system version
sys.path = [p for p in sys.path if 'SW/pyopencl' not in p]

# Import pyopencl first to ensure system version is loaded
import pyopencl as cl

# Add FireCore to path (don't clear sys.path to preserve system packages)
sys.path.append("../../")

from pyBall.OCL.FittingDriver import FittingDriver
from pyBall.OCL.NonBondFitting import run_energy_imshow

def parse_args():
    parser = argparse.ArgumentParser(description='Evaluate energy using pyOpenCL and plot')
    parser.add_argument('--xyz', type=str, default='wb97m_input/H2O-A1_H2O-D1.xyz', help='XYZ file')
    parser.add_argument('--atom-types', type=str, default='../../cpp/common_resources/AtomTypes.dat', help='Atom types file')
    parser.add_argument('--model', type=str, default='ENERGY_MorseQ_PAIR', help='Energy model name')
    parser.add_argument('--platform', type=int, default=1, help='OpenCL platform index (1=NVIDIA)')
    parser.add_argument('--device', type=int, default=0, help='OpenCL device index')
    parser.add_argument('--kcal', action='store_true', help='Convert energies to kcal/mol')
    parser.add_argument('--show', action='store_true', help='Show plot interactively')
    parser.add_argument('--save', type=str, default=None, help='Save plot to file')
    parser.add_argument('--verbose', type=int, default=1, help='Verbosity level')
    return parser.parse_args()

def main():
    args = parse_args()
    
    xyz_file = os.path.abspath(args.xyz)
    atom_types_file = os.path.abspath(args.atom_types)
    
    print(f"XYZ file: {xyz_file}")
    print(f"Atom types: {atom_types_file}")
    print(f"Model: {args.model}")
    print(f"OpenCL platform: {args.platform}, device: {args.device}")
    
    if not os.path.exists(xyz_file):
        print(f"Error: XYZ file not found: {xyz_file}")
        return
    
    # Use run_energy_imshow with NVIDIA GPU
    print("\nEvaluating energy using pyOpenCL (NVIDIA GPU)...")
    Vref, Vmod, Vdif, rv, Arow = run_energy_imshow(
        model_name=args.model,
        xyz_file=xyz_file,
        atom_types_file=atom_types_file,
        kcal=args.kcal,
        sym=True,
        bColorbar=True,
        verbose=args.verbose,
        show=args.show,
        save=args.save,
        cmap='bwr',
        lines=True
    )
    
    print(f"\nEnergy evaluation complete:")
    print(f"  Reference energy range: [{np.nanmin(Vref):.6f}, {np.nanmax(Vref):.6f}]")
    print(f"  Model energy range: [{np.nanmin(Vmod):.6f}, {np.nanmax(Vmod):.6f}]")
    print(f"  Difference range: [{np.nanmin(Vdif):.6f}, {np.nanmax(Vdif):.6f}]")
    
    if not args.show and args.save is None:
        print("\nPlot saved to default location or closed. Use --show to display interactively or --save to specify output file.")

if __name__ == "__main__":
    main()
