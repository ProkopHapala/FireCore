#!/usr/bin/env python3
import argparse
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import FitREQ as fit
import barrier_utils as bu

def run_endpoint(name, dof_file, xyz_path, Gref, seq, distances, angles, args, do_relax=False):
    print(f"\n--- Processing Endpoint {name} (Relax={do_relax}) ---")
    bu.setup_fitreq(dof_file, xyz_path, kMorse=args.kmorse, Lepairs=args.lepairs,  softClamp=True, sc_start=4.0, sc_max=6.0, imodel=7, Regularize=1)
    dofs, hyper = bu.parse_dof_file(dof_file)
    if do_relax: bu.relax_tethered(dofs, kMorse=args.kmorse, Lepairs=args.lepairs, nstep=args.nstep)
    run_params = dict(nstep=args.nstep, Fmax=1e-8, dt=0.05, max_step=0.1, damping=0.0)
    Gmodel = fit.compute_model_grid(xyz_path, seq, Gref.shape, do_fit=False, bAddEpairs=True, run_params=run_params, bOutXYZ=False)
    prefix = f"{args.out_prefix}_{name}_{'relax' if do_relax else 'rigid'}"
    title  = f"Endpoint {name} ({'Relaxed' if do_relax else 'Rigid'})"
    fit.plot_compare(Gref, Gmodel, angles, distances, title, save_prefix=prefix + "_2d",line=True,kcal=True, save_fmt="png")

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
    args = parser.parse_args()
    
    # Parse reference mapping once
    Gref, seq, axis, distances, angles = fit.parse_xyz_mapping(args.xyz)
    
    run_endpoint("A", args.dofA, args.xyz, Gref, seq, distances, angles, args, do_relax=False)
    run_endpoint("B", args.dofB, args.xyz, Gref, seq, distances, angles, args, do_relax=False)
    
    if args.relax:
        run_endpoint("A", args.dofA, args.xyz, Gref, seq, distances, angles, args, do_relax=True)
        run_endpoint("B", args.dofB, args.xyz, Gref, seq, distances, angles, args, do_relax=True)

if __name__ == "__main__":
    main()
