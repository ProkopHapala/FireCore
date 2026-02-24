#!/usr/bin/env python3
import argparse
import sys
import os

sys.path.append("../../")
import pyBall.barrier_utils as bu

def main():
    parser = argparse.ArgumentParser(description="Compare two DOF endpoints on a 2D map")
    parser.add_argument("--xyz",        type=str,   default="../tFitREQ_PN/wb97m-split/H2O-A1_H2O-D1-z.xyz")
    parser.add_argument("--dofA",       type=str,   default="dofSelection_saved_Fit1.dat")
    parser.add_argument("--dofB",       type=str,   default="dofSelection_MorseSR_H2O.dat")
    parser.add_argument("--kmorse",     type=float, default=1.7)
    parser.add_argument("--lepairs",    type=float, default=0.9)
    parser.add_argument("--nstep",      type=int,   default=100)
    parser.add_argument("--out_prefix", type=str,   default="compare")
    parser.add_argument("--relax", action="store_true", help="Also run constrained relaxation")
    parser.add_argument("--pn", action="store_true", help="Use FitREQ_PN instead of FitREQ")
    args = parser.parse_args()
    
    bu.init(use_pn=args.pn)
    fit = bu.fit
    
    # Parse reference mapping once
    if hasattr(fit, 'parse_xyz_mapping'):
        Gref, seq, axis, distances, angles = fit.parse_xyz_mapping(args.xyz)
        
        bu.run_compare_endpoints("A", args.dofA, args.xyz, Gref, seq, distances, angles, args.kmorse, args.lepairs, args.nstep, args.out_prefix, do_relax=False)
        bu.run_compare_endpoints("B", args.dofB, args.xyz, Gref, seq, distances, angles, args.kmorse, args.lepairs, args.nstep, args.out_prefix, do_relax=False)
        
        if args.relax:
            bu.run_compare_endpoints("A", args.dofA, args.xyz, Gref, seq, distances, angles, args.kmorse, args.lepairs, args.nstep, args.out_prefix, do_relax=True)
            bu.run_compare_endpoints("B", args.dofB, args.xyz, Gref, seq, distances, angles, args.kmorse, args.lepairs, args.nstep, args.out_prefix, do_relax=True)
    else:
        print("Warning: parse_xyz_mapping is not available in PN backend. Skipping grid plot comparison.")

if __name__ == "__main__":
    main()
