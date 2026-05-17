#!/usr/bin/env python3
import os
import sys
import argparse

sys.path.append("../../")
import pyBall.barrier_utils as bu

def main():
    parser = argparse.ArgumentParser(description="Grid search over hyperparameters")
    parser.add_argument("--xyz", type=str, default="../tFitREQ_PN/wb97m-split/H2O-A1_H2O-D1-y.xyz")
    parser.add_argument("--dof_file", type=str, default="dofSelection_MorseSR_H2O_epairOnly.dat")
    parser.add_argument("--out_dir", type=str, default="grid_search_out")
    parser.add_argument("--pn", action="store_true", help="Use FitREQ_PN instead of FitREQ")
    args = parser.parse_args()

    bu.init(use_pn=args.pn)
    
    kMorse_vals       = [1.6, 1.7, 1.8]
    Lepairs_vals      = [0.6, 0.8, 1.0, 1.2]
    weight_alpha_vals = [4.0]
    hScale_vals       = [0.0, 0.5, 1.0]

    bu.run_grid_search(args.xyz, args.dof_file, args.out_dir, kMorse_vals, Lepairs_vals, weight_alpha_vals, hScale_vals)

if __name__ == "__main__":
    main()
