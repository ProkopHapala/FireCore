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
    parser = argparse.ArgumentParser(description="Plot and compare 2D energy surfaces from a single .xyz file")

    parser.add_argument("--mode", choices=["plot", "model", "fit"], default="plot", help="Action: plot ref only, compare model (no fit), or fit then compare")
    parser.add_argument("-i", "--input",           default="/home/prokop/Desktop/CARBSIS/PEOPLE/Paolo/HbondFit_small_mols_2025_08_15/confs/wb97m-split/H2O-A1_H2O-D1-y.xyz", help="Input .xyz file (single movie)")
    #parser.add_argument("-i", "--input",           default="/home/prokop/Desktop/CARBSIS/PEOPLE/Paolo/HbondFit_small_mols_2025_08_15/confs/wb97m-splitm/HF-A1_HF-D1-z.xyz", help="Input .xyz file (single movie)")
    #parser.add_argument("--dof-selection",         default="dofSelection_MorseSR.dat", help="DOF selection file (defaults to Morse/LJ based on --lj)")
    parser.add_argument("--dof-selection",         default="dofSelection_MorseSR_H2O.dat", help="DOF selection file (defaults to Morse/LJ based on --lj)")
    parser.add_argument("--verbosity", type=int,   default=2,    help="Verbosity for FitREQ")
    parser.add_argument("--nstep",     type=int,   default=100,  help="Fitting steps")
    parser.add_argument("--fmax",      type=float, default=1e-8, help="Target force max for fitting")
    parser.add_argument("--dt",        type=float, default=0.01, help="Integrator dt")
    parser.add_argument("--max-step",  type=float, default=0.05, help="Max step")
    parser.add_argument("--damping",   type=float, default=0.0,  help="Damping")
    parser.add_argument("--lj",        action="store_true",      help="Use LJ instead of Morse presets")
    parser.add_argument("--save",      type=str, default=None,     help="Path to save the plot (PNG)")
    parser.add_argument("--no-epairs", action="store_true",      help="Disable epair terms when loading XYZ")
    parser.add_argument("--no-show",   action="store_true",      help="Do not show the figure")
    parser.add_argument("--line",      action="store_true",      help="Also plot r_min(angle) and E_min(angle) lines")
    args = parser.parse_args()

    bMorse_local = not args.lj
    Gref, seq, axis, distances, angles = fit.parse_xyz_mapping(args.input)
    title = os.path.basename(args.input)

    if args.mode == "plot":
        fit.plot_compare(Gref, None, angles, distances, title, save_prefix=args.save, show=not args.no_show, line=args.line)
        exit()

    # model or fit: compute model grid and compare
    run_params = dict(nstep=args.nstep, Fmax=args.fmax, dt=args.dt, max_step=args.max_step, damping=args.damping, bClamp=False)
    Gmodel = fit.compute_model_grid(
        args.input,
        seq,
        Gref.shape,
        do_fit=(args.mode == "fit"),
        verbosity=args.verbosity,
        bMorse=bMorse_local,
        dof_selection=args.dof_selection,
        bAddEpairs=(not args.no_epairs),
        run_params=run_params,
    )

    # Local compare plotting with requested style and optional lines
    save_prefix = None
    if args.save:
        save_prefix = args.save[:-4] if args.save.endswith('.png') else args.save
    fit.plot_compare(Gref, Gmodel, angles, distances, title, save_prefix=save_prefix, show=not args.no_show, line=args.line)
