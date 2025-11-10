#!/usr/bin/python3 -u

import sys
import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt
import tkinter
import time
import argparse


sys.path.append("../../")
from pyBall import FitREQ as fit

fit.plt = plt

np.set_printoptions(linewidth=300)

defalt_inputs=[
    #"H2O-A1_H2O-D1-y.xyz",   # sample out-of   plane epairs of O_3 atom in H2O with H2O probe
    #"H2O-A1_HCN-D1-y.xyz",   # sample in-plane plane epairs of O_3 atom in H2O with HCN probe
    #"H2O-A1_HF-D1-y.xyz",    # sample in-plane plane epairs of O_3 atom in H2O with HF probe
    
    "CH2O-A1_H2O-D1-z.xyz",      # sample in-plane plane epairs of O_2 atom in CH2O with H2O probe
    #"CH2O-A1_HCN-D1-z.xyz",      # sample in-plane plane epairs of O_2 atom in CH2O with HCN  probe
    #"CH2O-A1_HF-D1-z.xyz",       # sample in-plane plane epairs of O_2 atom in CH2O with HF   probe
    
    #"HCN-A1_HCN-D1-z.xyz", 
    
    #"NH3-A1_HCN-D1-z.xyz",

    #"HF-A1_H2O-D1-z.xyz",        # sample in-plane plane epairs of F atom in HF with H2O   probe
    #"HF-A1_HCN-D1-z.xyz",        # sample in-plane plane epairs of F atom in HF with HCN   probe
    #"HF-A1_HF-D1-z.xyz",         # sample in-plane plane epairs of F atom in HF with HF    probe
]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot and compare 2D energy surfaces from multiple .xyz files")
    parser.add_argument("-i", "--inputs", nargs='*', default=None, help="List of input .xyz files (single movie). If empty, uses default set")
    #parser.add_argument("--dir", type=str, default='/home/prokophapala/Desktop/CARBSIS/wb97m-split/', help="Directory where input files are located")
    #parser.add_argument("--dir", type=str, default="/home/prokop/Desktop/CARBSIS/PEOPLE/Paolo/HbondFit_small_mols_2025_08_15/confs/wb97m-split/", help="Directory where input files are located")
    parser.add_argument("--dir", type=str, default="/home/niko/work/HBOND/REFERENCE/2-pairs_small_small/4-to_firecore/confs_wb97m/", help="Directory where input files are located")
    parser.add_argument("--dof-selection",         default="dofSelection_MorseSR_H2O_CH2O.dat", help="DOF selection file")
    
    #parser.add_argument("--dof-selection",         default="dofSelection_MorseSR_H2O_CH2O_HF_HCN.dat", help="DOF selection file")

    #parser.add_argument("--dof-selection",         default="dofSelection_MorseSR_H2O_CH2O_HF_HCN-fix2.dat", help="DOF selection file")

    #parser.add_argument("--dof-selection",         default="dofSelection_MorseSR_H2O_CH2O_HF_HCN-fix.dat", help="DOF selection file")
    #parser.add_argument("--dof-selection",         default="dofSelection_MorseSR_nofit_H2O_CH2O_HF_HCN.dat", help="DOF selection file")

    #parser.add_argument("--dof-selection",         default="dofSelection_MorseSR_HF_HCN.dat", help="DOF selection file")

    parser.add_argument("--verbosity", type=int,   default=2,    help="Verbosity for FitREQ")
    parser.add_argument("--nstep",     type=int,   default=1000,   help="Fitting steps")
    parser.add_argument("--fmax",      type=float, default=1e-8, help="Target force max for fitting")
    parser.add_argument("--dt",        type=float, default=0.01, help="Integrator dt")
    parser.add_argument("--max-step",  type=float, default=0.05, help="Max step")
    parser.add_argument("--damping",   type=float, default=0.01,  help="Damping")
    # Global model params
    parser.add_argument("--kMorse",    type=float, default=1.8,  help="Global kMorse parameter")
    parser.add_argument("--Lepairs",   type=float, default=1.0,  help="Global Lepairs parameter")
    parser.add_argument("--lj",        type=int,   default=0,    help="Use LJ instead of Morse presets")
    # Weighting controls
    parser.add_argument("--user_weights",  type=int,   default=0,       help="Enable user weights during scan")
    parser.add_argument("--weight-a",      type=float, default=1.0,     help="Weight amplitude 'a' for exp weight func")
    parser.add_argument("--weight-alpha",  type=float, default=4.0,     help="Weight sharpness 'alpha' for exp weight func")
    parser.add_argument("--emin-min",      type=float, default=-0.02,   help="Emin threshold for weighting segments")
    parser.add_argument("--plot-dir",      type=str,   default='plots', help="Directory to save all output plots")
    parser.add_argument("--data-dir",      type=str,   default='data',  help="Directory to save all output data (2D maps, 1D lines)")
    parser.add_argument("--save-fmt",      choices=["both","npz","gnuplot"], default="both", help="Format for saved data")
    parser.add_argument("--save",          type=str,   default=None,    help="Base name for saving plots (without extension)")
    parser.add_argument("--epairs",        type=int,   default=1,       help="Disable epair terms when loading XYZ")
    parser.add_argument("--show",          type=int,   default=1,       help="show the figure")
    parser.add_argument("--line",          type=int,   default=1,       help="plot r_min(angle) and E_min(angle) lines")
    parser.add_argument("--out-xyz",       type=int,   default=0,       help="Output XYZ with fitted DOFs")
    parser.add_argument("--soft_clamp",    type=int,   default=1,       help="Enable soft clamp during scan")
    
    parser.add_argument("--regularization",    type=int,             default=0,                help="Enable regularization during scan")
    parser.add_argument("--clamp_thresholds",  nargs=2, type=float,  default=[0.2,0.5],       help="Soft clamp thresholds: start max")
    parser.add_argument("--kcal", type=int, default=1, help="Use kcal instead of eV")


    args = parser.parse_args()
    
    # Use default inputs if none provided
    if args.inputs is None or len(args.inputs) == 0: args.inputs = defalt_inputs
    print("Using default input files:", args.inputs )

    imodel = 1 if args.lj else 5

    fit.setVerbosity(args.verbosity, PrintDOFs=1, PrintfDOFs=1, PrintBeforReg=-1, PrintAfterReg=1)
    
    regularize_val = args.regularization
    fit.setup( imodel=imodel, EvalJ=1, WriteJ=1, Regularize=regularize_val, SaveJustElementXYZ=-1, SoftClamp=args.soft_clamp )

    if args.soft_clamp:
        fit.setGlobalParams( kMorse=args.kMorse, Lepairs=args.Lepairs, softClamp_start=args.clamp_thresholds[0], softClamp_max=args.clamp_thresholds[1] )
    else:
        fit.setGlobalParams( kMorse=args.kMorse, Lepairs=args.Lepairs )

    fit.loadTypes()
    fit.loadDOFSelection(fname=args.dof_selection)

    # --- Batch Loading
    sample_counts = []
    all_Erefs = []
    all_x0s = []
    for i, f in enumerate(args.inputs):
        fname = os.path.join(args.dir, f)
        print(f"Loading {fname}")
        bAppend = (i > 0)
        n = fit.loadXYZ(fname, bAddEpairs=(args.epairs == 1), bAppend=bAppend)
        sample_counts.append(n)
        if args.user_weights:
            Erefs, x0s = fit.read_xyz_data(fname)
            all_Erefs.append(Erefs)
            all_x0s.append(x0s)

    fit.getBuffs()
    
    if args.user_weights:
        Erefs_all = np.concatenate(all_Erefs)
        x0s_all   = np.concatenate(all_x0s)
        n_before = 5 if args.lj else 100
        weight_func = lambda E: fit.exp_weight_func(E, a=args.weight_a, alpha=args.weight_alpha)
        weights0, lens = fit.split_and_weight_curves(  Erefs_all, x0s_all,  n_before_min=n_before, weight_func=weight_func, EminMin=args.emin_min )
        fit.setWeights( weights0 )

    # --- Fitting
    trj_E, trj_F, trj_DOFs, _ = fit.setTrjBuffs(niter=args.nstep)
    DOFnames = fit.loadDOFnames(args.dof_selection)
    
    # Create output directories if needed
    if not os.path.exists(args.plot_dir): os.makedirs(args.plot_dir)
    if not os.path.exists(args.data_dir): os.makedirs(args.data_dir)
    
    fit.run( iparallel=0, ialg=1, nstep=args.nstep, Fmax=args.fmax, dt=args.dt, damping=args.damping,   max_step=-1,  bClamp=True )

    # Save trajectory plot
    traj_path = os.path.join(args.plot_dir, f"trajectory_{os.path.splitext(os.path.basename(args.dof_selection))[0]}.png")
    fit.plot_trj_dofs(trj_DOFs, DOFnames=DOFnames, title="Optimization trajectory", save_path=traj_path)

    Eerr, E_models, _ = fit.getEs( bDOFtoTypes=True, bEs=True )

    #print(Eerr.shape)
    print(E_models.shape)
    print("sample_counts ", sample_counts, np.sum(sample_counts))
    
    istart = 0
    for i, f in enumerate(args.inputs):
        fname = os.path.join(args.dir, f)
        base = os.path.splitext(f)[0]  # Remove .xyz extension
        title = base.replace('_', ' ')  # Clean title
        iend = istart + sample_counts[i]
        
        # Get grid mapping
        Gref, seq, axis, distances, angles = fit.parse_xyz_mapping(fname)
        Es_slice = E_models[istart:iend]
        
        # Map to grid with NaN padding
        Gmodel = np.empty_like(Gref); Gmodel[:] = np.nan
        nmap = min(len(Es_slice), len(seq))
        for k in range(nmap):
            idist, iang = seq[k]
            Gmodel[idist, iang] = Es_slice[k]
        plot_path = os.path.join(args.plot_dir, f"{base}.png")
        data_base = os.path.join(args.data_dir, base)
        fit.plot_compare(
            Gref, Gmodel, angles, distances, title,
            save_prefix=plot_path,
            line=bool(args.line),
            kcal=bool(args.kcal),
            save_data_prefix=data_base,
            save_fmt=args.save_fmt
        )
        istart = iend

    if args.show == 1: plt.show()
