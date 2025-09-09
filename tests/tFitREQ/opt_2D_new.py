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
#from pyBall import atomicUtils as au

fit.plt = plt

np.set_printoptions(linewidth=300)


"""
Command-line utility for visualizing and comparing 2D energy surfaces from a single .xyz file.

Modes:
  - plot: plot only the reference energies extracted from the .xyz comments
  - model: compute model energies for the frames (no fitting) and compare to reference
  - fit: run fitting, then compute model energies and compare to reference
"""

'''
H2O-A1_H2O-D1-y.xyz   : sample out-of   plane epairs of O_3 atom in H2O with H2O probe
H2O-A1_HCN-D1-y.xyz   : sample in-plane plane epairs of O_3 atom in H2O with HCN probe
H2O-A1_HF-D1-y.xyz    : sample in-plane plane epairs of O_3 atom in H2O with HF probe
--
CH2O-A1_H2O-D1-z      : sample in-plane plane epairs of O_2 atom in CH2O with H2O probe
CH2O-A1_HCN-D1-z.xyz  : sample in-plane plane epairs of O_2 atom in CH2O with HCN  probe
CH2O-A1_HF-D1-z.xyz   : sample in-plane plane epairs of O_2 atom in CH2O with HF   probe
--
HF-A1_H2O-D1-z.xyz    : sample in-plane plane epairs of F atom in HF with H2O   probe
HF-A1_HCN-D1-z.xyz    : sample in-plane plane epairs of F atom in HF with HCN   probe
HF-A1_HF-D1-z.xyz     : sample in-plane plane epairs of F atom in HF with HF    probe
'''
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot and compare 2D energy surfaces from a single .xyz file or scan DOFs")
    parser.add_argument("--mode", choices=["plot", "model", "fit", "scan"], default="fit", help="Action: plot ref only, compare model (no fit), fit then compare, or scan DOFs")
    #parser.add_argument("-i", "--input",           default="/home/prokophapala/Desktop/CARBSIS/wb97m-split/H2O-A1_H2O-D1-y.xyz", help="Input .xyz file (single movie)")
    parser.add_argument("-i", "--input",           default="/home/prokop/Desktop/CARBSIS/PEOPLE/Paolo/HbondFit_small_mols_2025_08_15/confs/wb97m-split/H2O-A1_H2O-D1-y.xyz", help="Input .xyz file (single movie)")
    parser.add_argument("--dof-selection",         default="dofSelection_MorseSR_H2O.dat", help="DOF selection file")
    parser.add_argument("--verbosity", type=int,   default=2,    help="Verbosity for FitREQ")
    parser.add_argument("--nstep",     type=int,   default=100, help="Fitting steps")
    parser.add_argument("--fmax",      type=float, default=1e-8, help="Target force max for fitting")
    parser.add_argument("--dt",        type=float, default=0.05, help="Integrator dt")
    parser.add_argument("--max-step",  type=float, default=0.1, help="Max step")
    parser.add_argument("--damping",   type=float, default=0.0,  help="Damping")
    # Global model params
    parser.add_argument("--kMorse",    type=float, default=1.7,  help="Global kMorse parameter")
    parser.add_argument("--Lepairs",   type=float, default=1.2,  help="Global Lepairs parameter")
    parser.add_argument("--model",     type=int,   default=7,    help="Model type: 1=LJ, 5=MorseQ_SR, 7=MorseQ_SR_boys  8=MorseQ_SR_softclamp ")
    # Weighting controls
    parser.add_argument("--n_before",        type=int,   default=100,    help="#points before min to weight")
    parser.add_argument("--weight-a",      type=float, default=1.0,    help="Weight amplitude 'a' for exp weight func")
    parser.add_argument("--weight-alpha",  type=float, default=4.0,    help="Weight sharpness 'alpha' for exp weight func")
    parser.add_argument("--emin-min",      type=float, default=-0.02,  help="Emin threshold for weighting segments")
    parser.add_argument("--save",          type=str,   default=None,   help="Path to save the plot (PNG)")
    parser.add_argument("--save-data-prefix", type=str, default=None,   help="Base path to save 2D grids and 1D lines (no extension)")
    parser.add_argument("--save-fmt", choices=["both","npz","gnuplot"], default="both", help="Format for saved data")
    parser.add_argument("--epairs",        type=int,   default=1,      help="Disable epair terms when loading XYZ")
    parser.add_argument("--show",          type=int,   default=1,      help="Do not show the figure")
    parser.add_argument("--line",          type=int,   default=1,      help="Do not plot r_min(angle) and E_min(angle) lines")
    parser.add_argument("--out-xyz",       type=int,   default=1,      help="Output XYZ with fitted DOFs")
    # Scan arguments
    parser.add_argument("--scan_dofs",         type=int,  nargs='+', default=None,             help="List of DOF indices to scan. If None, all from dof-selection are scanned.")
    parser.add_argument("--scan_range",        type=float, nargs=3,  default=[-1.0, 1.0, 100], help="Scan range: min max n_steps")
    parser.add_argument("--soft_clamp",        type=int,             default=1,                help="Enable soft clamp during scan")
    parser.add_argument("--user_weights",      type=int,             default=1,                help="Enable user weights during scan")
    parser.add_argument("--regularization",    type=int,             default=0,                help="Enable regularization during scan")
    parser.add_argument("--clamp_thresholds",  nargs=2, type=float,  default=[4.0, 6.0],       help="Soft clamp thresholds: start max")
    parser.add_argument("--kcal", type=int, default=1, help="Use kcal instead of eV")

    args = parser.parse_args()

    imodel = args.model

    fit.setVerbosity(args.verbosity, PrintDOFs=1, PrintfDOFs=1, PrintBeforReg=-1, PrintAfterReg=1)
    
    regularize_val = 0 if args.regularization == 0 else 1
    fit.setup( imodel=imodel, EvalJ=1, WriteJ=1, Regularize=regularize_val, SaveJustElementXYZ=-1 )

    if args.soft_clamp:
        fit.setGlobalParams( kMorse=args.kMorse, Lepairs=args.Lepairs, softClamp_start=args.clamp_thresholds[0], softClamp_max=args.clamp_thresholds[1] )
    else:
        fit.setGlobalParams( kMorse=args.kMorse, Lepairs=args.Lepairs )

    fit.loadTypes()
    fit.loadDOFSelection(fname=args.dof_selection)
    fit.loadXYZ( args.input, bAddEpairs=(args.epairs == 1) )
    fit.getBuffs()
    
    if args.user_weights:
        Erefs, x0s = fit.read_xyz_data(args.input)
        n_before = args.n_before
        weight_func = lambda E: fit.exp_weight_func(E, a=args.weight_a, alpha=args.weight_alpha)
        weights0, lens = fit.split_and_weight_curves( Erefs, x0s, n_before_min=n_before, weight_func=weight_func, EminMin=args.emin_min )
        fit.setWeights( weights0 )

    if args.mode == "scan":
        DOFnames = fit.loadDOFnames(args.dof_selection)
        scan_dofs = args.scan_dofs
        if scan_dofs is None:
            scan_dofs = list(range(fit.nDOFs))
        xs = np.linspace(args.scan_range[0], args.scan_range[1], int(args.scan_range[2]))
        for iDOF in scan_dofs:
            Es, Fs = fit.scanParam(iDOF, xs)
            plt.figure()
            plt.plot(xs, Es)
            plt.xlabel("DOF value")
            plt.ylabel("Energy")
            plt.title(f"DOF Scan: {DOFnames[iDOF]}")
            plt.grid(True)
            if args.save:
                plt.savefig(f"{args.save}_DOF_{iDOF}.png")
    elif args.mode == "plot":
        Gref, seq, axis, distances, angles = fit.parse_xyz_mapping(args.input)
        title = os.path.basename(args.input)
        fit.plot_compare(Gref, None, angles, distances, title,
                         save_prefix=args.save,
                         show=args.show==1,
                         line=args.line == 1,
                         kcal=args.kcal,
                         save_data_prefix=args.save_data_prefix,
                         save_fmt=args.save_fmt)
    else: # fit or model
        Gref, seq, axis, distances, angles = fit.parse_xyz_mapping(args.input)
        title = os.path.basename(args.input)
        run_params = dict(nstep=args.nstep, Fmax=args.fmax, dt=args.dt, max_step=args.max_step, damping=args.damping )
        # Allocate trajectory buffers if we will fit, so the run() fills them
        if args.mode == "fit":  trj_E, trj_F, trj_DOFs, _ = fit.setTrjBuffs(niter=args.nstep)
        Gmodel = fit.compute_model_grid( args.input, seq, Gref.shape, do_fit=(args.mode == "fit"), bAddEpairs=(args.epairs == 1), run_params=run_params, bOutXYZ=args.out_xyz )
        if args.mode == "fit":
            DOFnames = fit.loadDOFnames(args.dof_selection)
            fit.plot_trj_dofs(trj_DOFs, DOFnames=DOFnames, title="Optimization trajectory")
        # Local compare plotting with requested style and optional lines
        save_prefix = None
        if args.save: save_prefix = args.save[:-4] if args.save.endswith('.png') else args.save
        fit.plot_compare(Gref, Gmodel, angles, distances, title,
                         save_prefix=save_prefix,
                         show=args.show==1,
                         line=args.line == 1,
                         kcal=args.kcal,
                         save_data_prefix=args.save_data_prefix,
                         save_fmt=args.save_fmt)
    if args.show == 1:  plt.show()