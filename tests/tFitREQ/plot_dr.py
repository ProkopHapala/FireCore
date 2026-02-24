#!/usr/bin/env python3
import argparse
import sys
import os

sys.path.append("../../")
import pyBall.barrier_utils as bu

def main():
    parser = argparse.ArgumentParser(description="Dimensionality Reduction on Trajectories")
    parser.add_argument("--traj_file", type=str, default="path_trajectory.npz", help="Path to trajectory/ensemble NPZ")
    parser.add_argument("--trajectory_npz", type=str, default=None, help="Alias for --traj_file")
    parser.add_argument("--out_prefix", type=str, default="dr")
    parser.add_argument("--pn", action="store_true", help="Use FitREQ_PN instead of FitREQ")
    args = parser.parse_args()

    if args.trajectory_npz:
        args.traj_file = args.trajectory_npz

    bu.run_dr(args.traj_file, args.out_prefix)

if __name__ == "__main__":
    main()
