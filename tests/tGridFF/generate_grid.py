#!/usr/bin/env python3
"""
Generate a grid for force field calculations from a substrate XYZ file.

This script processes a substrate structure and generates a grid with the specified voxel size
for force field calculations using OpenCL acceleration.
"""
import sys
import os
import argparse
from pathlib import Path

# Add parent directory to path to allow importing pyBall
sys.path.append("../../")
from pyBall.tests import ocl_GridFF_new as gff
import pyopencl as cl
import os.path

# Create symlinks if they don't exist
if not os.path.exists("common_resources"):
    os.symlink("../../cpp/common_resources", "common_resources")
if not os.path.exists("data"):
    os.symlink("../../cpp/common_resources", "data")

# Define paths using the symlinks
ELEMENT_TYPES = "common_resources/ElementTypes.dat"
XYZ_DIR = "common_resources/xyz"

def parse_arguments():
    parser = argparse.ArgumentParser(description='Generate grid for force field calculations')
    parser.add_argument('--fname', required=True, help='Input XYZ filename (with or without .xyz)')
    parser.add_argument('--desired_voxel', type=float, default=0.1, help='Voxel size in Angstroms')
    # parser.add_argument('--element_types', default="./data/ElementTypes.dat", help='Element types file')
    # parser.add_argument('--output_prefix', default="double3", help='Output file prefix')
    return parser.parse_args()

def main():
    args = parse_arguments()
    
    # Process input filename
    input_path = Path(args.fname)
    if input_path.suffix != '.xyz':
        input_path = input_path.with_suffix('.xyz')
    
    if not input_path.exists():
        # Check in xyz directory if file not found in current directory
        data_path = Path(XYZ_DIR) / input_path.name
        if data_path.exists():
            input_path = data_path
        else:
            print(f"Error: Input file {args.fname} not found in current directory or {XYZ_DIR}/")
            sys.exit(1)
    
    print(f"Processing: {input_path}")
    print(f"Using element types from: {ELEMENT_TYPES}")
    gff.test_gridFF_ocl(
        fname=str(input_path),
        Element_Types_name=ELEMENT_TYPES,
        save_name="double3",
        job="PLQ",
        desired_voxel=args.desired_voxel
    )

if __name__ == "__main__":
    main()





