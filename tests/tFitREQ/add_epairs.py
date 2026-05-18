#!/usr/bin/env python3
"""
Add electron pairs and/or sigma holes to multi-frame XYZ dimer files.

Processes each frame separately, splitting at n0 (first n0 atoms = molecule A,
remaining = molecule B). Updates n0 in the header after adding atoms.

AtomTypes.dat is read to determine which atoms get epairs (nepair column)
and which have sigma-hole epair types. The epair distance (Lepair) and
sigma-hole distance are configurable.

Usage examples:
  # Process one file:
  python tests/tFitREQ/add_epairs.py -i scan.xyz -o scan_ep.xyz

  # Process all .xyz files in a directory:
  python tests/tFitREQ/add_epairs.py -i input_dir/ -o output_dir/

  # Custom distances and sigma holes:
  python tests/tFitREQ/add_epairs.py -i input_dir/ -o out_dir/ 
      --mode both --lepair 1.2 --sigma-dist 0.6

  # Use a specific AtomTypes.dat:
  python tests/tFitREQ/add_epairs.py -i scan.xyz -o out.xyz 
      --atypes tests/tFitREQ_PN/data/AtomTypes.dat

  # Show help:
  python tests/tFitREQ/add_epairs.py -h
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from pyBall import atomicUtils as au
from pyBall import elements
from pyBall.AtomicSystem import AtomicSystem
from pyBall.FitREQutils import (read_atom_types, process_one_file, find_data_path )

# ── CLI ──

def main():
    parser = argparse.ArgumentParser( description="Add electron pairs and sigma holes to multi-frame XYZ dimer files", formatter_class=argparse.RawTextHelpFormatter, )
    parser.add_argument("-i", "--input", required=True, help="Input .xyz file or directory of .xyz files")
    parser.add_argument("-o", "--output", default="output",  help="Output file or directory (default: output)")
    parser.add_argument("--mode", choices=["epairs", "sigma", "both"], default="epairs",  help="What to add (default: epairs)")
    parser.add_argument("--lepair", type=float, default=1.0,  help="Epair distance from host in Angstrom (default: 1.0)")
    parser.add_argument("--sigma-dist", type=float, default=0.5,  help="Sigma-hole distance from host in Angstrom (default: 0.5)")
    parser.add_argument("--atypes", default=None,  help="Path to AtomTypes.dat (auto-searched if not given)")
    parser.add_argument("--simple-names", action="store_true",   help="Output simple element names (N, O, H) instead of full type names (N_3, O_2, H_O)")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()

    do_epairs = args.mode in ("epairs", "both")
    do_sigma = args.mode in ("sigma", "both")

    # ── Find AtomTypes.dat ──
    script_dir = os.path.dirname(os.path.abspath(__file__))
    if args.atypes:
        atypes_path = args.atypes
    else:
        data_dir = find_data_path(script_dir)
        if data_dir is None:
            print("Error: cannot find AtomTypes.dat. Use --atypes to specify path.")
            sys.exit(1)
        atypes_path = os.path.join(data_dir, "AtomTypes.dat")

    if not os.path.isfile(atypes_path):
        print(f"Error: {atypes_path} not found")
        sys.exit(1)

    atom_types = read_atom_types(atypes_path)
    print(f"Loaded {len(atom_types)} atom types from {atypes_path}")

    # ── Collect input files ──
    inpath = args.input
    outpath = args.output

    if os.path.isfile(inpath):
        files = [(inpath, outpath)]
    elif os.path.isdir(inpath):
        os.makedirs(outpath, exist_ok=True)
        files = []
        for fn in sorted(os.listdir(inpath)):
            if not fn.endswith(".xyz"):
                continue
            if "_bak" in fn or "_Epairs" in fn:
                continue
            files.append((os.path.join(inpath, fn), os.path.join(outpath, fn)))
        if not files:
            print(f"No .xyz files found in {inpath}")
            sys.exit(1)
        print(f"Found {len(files)} .xyz files in {inpath}")
    else:
        print(f"Error: {inpath} is not a file or directory")
        sys.exit(1)

    # ── Process ──
    total_frames = 0
    total_A = 0
    total_B = 0

    for src, dst in files:
        print(f"Processing: {os.path.basename(src)}")
        nf, nA, nB = process_one_file(  src, dst, atom_types, args.lepair, args.sigma_dist, do_epairs, do_sigma, args.verbose, args.simple_names,)
        if nf > 0: print(f"  → {nf} frames, +{nA}/+{nB} (molA/molB) → {dst}")
        total_frames += nf
        total_A += nA
        total_B += nB

    print(f"\nTotal: {total_frames} frames across all files")
    print(f"  Added atoms: molA={total_A}, molB={total_B}")
    print(f"  Mode: {args.mode}, Lepair={args.lepair}, sigma-dist={args.sigma_dist}")


if __name__ == "__main__":
    main()
