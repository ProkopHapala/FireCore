#!/usr/bin/env python3
import argparse
import sys
import os

sys.path.append("../../")
import pyBall.barrier_utils as bu

def main():
    parser = argparse.ArgumentParser(description="Plot analytical directional barriers")
    parser.add_argument( "-L", type=float, default=1.0  )
    parser.add_argument( "-R", type=float, default=3.0  )
    parser.add_argument( "-E", type=float, default=0.01 )
    parser.add_argument( "-k", type=float, default=1.6  )
    parser.add_argument( "-M", "--bMorse", action='store_true' )
    parser.add_argument( "--out", type=str, default=None )
    args = parser.parse_args()

    out_file = args.out if args.out else f"Directional_Barrier_bMorse{int(args.bMorse)}.png"
    bu.run_directional_barrier(args.L, args.R, args.E, args.k, args.bMorse, out_file)

if __name__ == "__main__":
    main()
