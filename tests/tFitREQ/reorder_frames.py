#!/usr/bin/env python3
"""
Reorder XYZ frames by angle (y or z field) then distance (x0 field).

Usage:
  python reorder_frames.py -i file.xyz
  python reorder_frames.py -i file.xyz -o reordered.xyz
  python reorder_frames.py --dir /path/to/files/
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from pyBall.FitREQutils import reorder_frames_by_angle

def main():
    parser = argparse.ArgumentParser(description="Reorder XYZ frames by angle then distance")
    parser.add_argument("-i", "--input", required=True, help="Input .xyz file or directory")
    parser.add_argument("-o", "--output", default=None, help="Output .xyz file (default: overwrite input)")
    parser.add_argument("--dir", action="store_true", help="Process all .xyz files in input directory")
    args = parser.parse_args()
    
    if args.dir:
        if not os.path.isdir(args.input):
            print(f"Error: {args.input} is not a directory")
            return
        for fname in sorted(os.listdir(args.input)):
            if fname.endswith('.xyz'):
                inpath = os.path.join(args.input, fname)
                reorder_frames_by_angle(inpath)
    else:
        if not os.path.isfile(args.input):
            print(f"Error: {args.input} not found")
            return
        reorder_frames_by_angle(args.input, args.output)


if __name__ == "__main__":
    main()
