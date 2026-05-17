#!/usr/bin/env python3
import argparse
import sys
import os

sys.path.append("../../")
import pyBall.barrier_utils as bu

def main():
    parser = argparse.ArgumentParser(description="1D String interpolation and Hessian analysis")
    parser.add_argument("--xyz", type=str, default="../tFitREQ_PN/wb97m-split/H2O-A1_H2O-D1-z.xyz")
    parser.add_argument("--dofA", type=str, required=True, help="Path to first DOF file")
    parser.add_argument("--dofB", type=str, required=True, help="Path to second DOF file")
    parser.add_argument("--n_points", type=int, default=11)
    parser.add_argument("--kmorse", type=float, default=1.7)
    parser.add_argument("--lepairs", type=float, default=0.9)
    parser.add_argument("--springK", type=float, default=100.0)
    parser.add_argument("--nstep", type=int, default=50)
    parser.add_argument("--out_prefix", type=str, default="path")
    parser.add_argument("--relax", action="store_true", help="Run constrained relaxation")
    parser.add_argument("--hessian", action="store_true", help="Compute Hessian at alpha=0, 0.5, 1")
    parser.add_argument("--pn", action="store_true", help="Use FitREQ_PN instead of FitREQ")
    args = parser.parse_args()

    bu.init(use_pn=args.pn)

    if not os.path.exists(args.dofA):
        print(f"ERROR: {args.dofA} does not exist.")
        sys.exit(1)
    if not os.path.exists(args.dofB):
        print(f"ERROR: {args.dofB} does not exist.")
        sys.exit(1)

    bu.run_string_scan(
        xyz_path=args.xyz,
        dofA=args.dofA,
        dofB=args.dofB,
        n_points=args.n_points,
        kmorse=args.kmorse,
        lepairs=args.lepairs,
        springK=args.springK,
        nstep=args.nstep,
        out_prefix=args.out_prefix,
        do_relax=args.relax,
        do_hessian=args.hessian
    )

if __name__ == "__main__":
    main()
