#!/usr/bin/python3 -u

import sys
import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt
import tkinter
import time
import argparse
from pyBall import FitREQ_PN as fit

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fit model parameters to reproduce 2D energy surfaces from multiple .xyz files")
    # setVerbosity
    parser.add_argument("--verbosity",          type=int, default=3, help="Verbosity for FitREQ")
    parser.add_argument("--idebug",             type=int, default=0, help="Flag for debug")
    parser.add_argument("--printdofs",          type=int, default=1, help="Set PrintDOFs")
    parser.add_argument("--printfdofs",         type=int, default=1, help="Set PrintfDOFs")
    parser.add_argument("--printbeforreg",      type=int, default=1, help="Set PrintBeforReg")
    parser.add_argument("--printafterreg",      type=int, default=1, help="Set PrintAfterReg")
    parser.add_argument("--printoverrepulsive", type=int, default=0, help="Set PrintOverRepulsive")
    # setModel
    parser.add_argument("--ivdw",    type=int,   default=4,   help="Set ivdW (0=no vdW, 1=LJ, 2=LJr8, 3=LJr9, 4=Morse, 5=Buck")
    parser.add_argument("--icoul",   type=int,   default=1,   help="Set iCoul (0=no Coul, 1=point charges, ...")
    parser.add_argument("--ihbond",  type=int,   default=0,   help="Set iHBond (0=no HBond correction, 1=H1 correction, 2=H2 correction, 3=H1 and H2 correction")
    parser.add_argument("--epairs",  type=int,   default=1,   help="Enable/disable Epair terms (it applies both when loading XYZ and when calculating interactions)")
    parser.add_argument("--iepairs", type=int,   default=0,   help="Set iEpairs (0=no interaction, 1=SR, 2=SR2, 3=SR3, 4=SR4)")
    parser.add_argument("--lepairs", type=float, default=0.5, help="Set distance of Epairs from host atom")
    parser.add_argument("--bpn",     type=int,   default=1,   help="Set boolean for PN functions")
    parser.add_argument("--kmorse",  type=float, default=1.8, help="Set curvature of Morse potential")
    parser.add_argument("--sr4cut",  type=float, default=1.0, help="Length of Epair interaction in SR4")
    parser.add_argument("--sr4n",    type=int,   default=2,   help="Order of polynomial of Epair interaction in SR4")
    parser.add_argument("--sr4m",    type=int,   default=2,   help="Order of polynomial of Epair interaction in SR4")
    parser.add_argument("--svdw",    type=float, default=1.0, help="Strength of vdW interaction")
    parser.add_argument("--scoul",   type=float, default=1.0, help="Strength of Coulomb interaction")
    parser.add_argument("--shcorr",  type=float, default=1.0, help="Strength of H-bond interaction")
    parser.add_argument("--sepairs", type=float, default=1.0, help="Strength of Epair interaction")
    # setInput
    parser.add_argument("--inputs",        type=str, default=None,                    help="List of input .xyz files", nargs='*')
    parser.add_argument("--inputs_dir",    type=str, default='data',                  help="Directory where input files are located")
    parser.add_argument("--dof_selection", type=str, default=None,                    help="DOF selection file")
    parser.add_argument("--fetypes",       type=str, default='data/ElementTypes.dat', help="File with Element types")
    parser.add_argument("--fatypes",       type=str, default='data/AtomTypes.dat',    help="File with Atom types")
    # setPenalty
    parser.add_argument("--clamp",           type=int,   default=1,      help="Hardly restrain the values of parameters during optimization")
    parser.add_argument("--regularize",      type=int,   default=0,      help="Apply harmonic potentials to keep parameters within a specified range")
    parser.add_argument("--regcountweight",  type=int,   default=0,      help="Apply further normalization on regularizing potentials")
    parser.add_argument("--softclamp",       type=int,   default=0,      help="Apply soft clamping to specific samples")
    parser.add_argument("--softclamp_start", type=float, default=0.0867, help="Minimum threshold in eV where soft clamping starts to kick in")
    parser.add_argument("--softclamp_max",   type=float, default=0.1735, help="Maximum value in eV of the penalty function")
    parser.add_argument("--user_weights",    type=int,   default=1,      help="Enable user weights during scan")
    parser.add_argument("--weight_a",        type=float, default=1.0,    help="Weight amplitude 'a' for exp weight func")
    parser.add_argument("--weight_alpha",    type=float, default=4.0,    help="Weight sharpness 'alpha' for exp weight func")
    # setOptimization
    parser.add_argument("--nstep",     type=int,   default=2000, help="Fitting steps")
    parser.add_argument("--ialg",      type=int,   default=2,    help="Optimization algorithm, 0=gradient descent, 1=damped dynamics, 2=gradient descent Barzilai-Borwein short step, 3=gradient descent Barzilai-Borwein long step")
    parser.add_argument("--fmax",      type=float, default=1e-8,  help="Target force max for fitting")
    parser.add_argument("--dt",        type=float, default=0.02,  help="Integrator dt")
    parser.add_argument("--max_step",  type=float, default=0.05,  help="Max step")
    parser.add_argument("--damping",   type=float, default=0.01,  help="Damping")
    parser.add_argument("--iparallel", type=int,   default=0,     help="Use OpenMP for penalty evaluation or not (0=serial, 1=OpenMP)")
    # setOutput
    parser.add_argument("--plot_dir",           type=str, default=None,      help="Directory to save all output plots")
    parser.add_argument("--out_dir",            type=str, default=None,      help="Directory to save all output data (2D maps, 1D lines)")
    parser.add_argument("--save_fmt",           type=str, default="gnuplot", help="Format for saved data", choices=["both","npz","gnuplot"])
    parser.add_argument("--save",               type=str, default=None,      help="Base name for saving plots (without extension)")
    parser.add_argument("--show",               type=int, default=0,         help="show the figure")
    parser.add_argument("--line",               type=int, default=1,         help="plot r_min(angle) and E_min(angle) lines")
    parser.add_argument("--outxyz",             type=int, default=0,         help="Output XYZ with Epairs")
    parser.add_argument("--savejustelementxyz", type=int, default=0,         help="Print only element names (not types) in output XYZ with Epairs")
    parser.add_argument("--outxyz_fname",       type=str, default=None,      help="Filename for output XYZ with Epairs")
    parser.add_argument("--kcal",               type=int, default=1,         help="Use kcal instead of eV in the output energies")
    args = parser.parse_args()
    
    # setVerbosity
    fit.setVerbosity(verbosity=args.verbosity, idebug=args.idebug, PrintDOFs=args.printdofs, PrintfDOFs=args.printfdofs,
                     PrintBeforReg=args.printbeforreg, PrintAfterReg=args.printafterreg, PrintOverRepulsive=args.printoverrepulsive)
    
    # setModel
    fit.setModel(ivdW=args.ivdw, iCoul=args.icoul, iHbond=args.ihbond, Epairs=args.epairs, iEpairs=args.iepairs, kMorse=args.kmorse,
                 Lepairs=args.lepairs, bPN=(args.bpn == 1), svdW=args.svdw, sCoul=args.scoul, sHcorr=args.shcorr, sEpairs=args.sepairs,
                 SR4cut=args.sr4cut, SR4n=args.sr4n, SR4m=args.sr4m)

    # setInput
    fit.loadTypes(fEtypes=args.fetypes, fAtypes=args.fatypes)
    fit.loadDOFSelection(fname=args.dof_selection)
    sample_counts = []
    total_loaded = 0
    all_Erefs = []
    all_pairs = []
    for i, f in enumerate(args.inputs):
        fname = os.path.join(args.inputs_dir, f)
        print(f"Loading {fname}")
        bAppend = (i > 0)
        n_total = fit.loadXYZ(fname, bAddEpairs=(args.epairs == 1), bOutXYZ=(args.outxyz == 1), bSaveJustElementXYZ=(args.savejustelementxyz == 1),
                              OutXYZ_fname=args.outxyz_fname, bEvalOnlyCorrections=False, bAppend=bAppend)
        n_delta = n_total - total_loaded
        if n_delta < 0:
            print(f"Warning: loadXYZ returned decreasing total count ({n_total} < {total_loaded}); forcing delta=0")
            n_delta = 0
        sample_counts.append(n_delta)
        total_loaded = n_total
        if args.user_weights:
            Erefs, pairs = fit.read_xyz_data_new(fname)
            all_Erefs.append(Erefs)
            all_pairs.append(pairs)
    fit.getBuffs()
    
    # setPenalty
    fit.setPenalty( Clamp=args.clamp, Regularize=args.regularize, AddRegError=args.regularize, RegCountWeight=args.regcountweight,
                    SoftClamp=args.softclamp, softClamp_start=args.softclamp_start, softClamp_max=args.softclamp_max )
    if args.user_weights:
        Erefs_all = np.concatenate(all_Erefs)
        pairs_all = np.concatenate(all_pairs)
        weight_func = lambda E: fit.exp_weight_func_new(E, a=args.weight_a, alpha=args.weight_alpha)
        weights0 = fit.split_and_weight_curves_new(Erefs_all, pairs_all, weight_func=weight_func)
        fit.setWeights( weights0 )

    # Load trajectory and DOFs
    trj_E, trj_F, trj_DOFs, _ = fit.setTrjBuffs(niter=args.nstep)
    DOFnames = fit.loadDOFnames(args.dof_selection)
    
    # Create output directoryies if needed
    if not os.path.exists(args.out_dir): os.makedirs(args.out_dir)

    # Run fitting
    fit.run_PN(iparallel=args.iparallel, ialg=args.ialg, nstep=args.nstep, Fmax=args.fmax, dt=args.dt, damping=args.damping, max_step=args.max_step)

    # Save trajectory plot
    if len(DOFnames) > 0:
        fit.save_trj_dofs(trj_DOFs, DOFnames=DOFnames, folder=args.out_dir)

    # Get result fitting
    Eerr, E_models, _ = fit.getEs(bDOFtoTypes=True, bEs=True)

    # Count NaNs
    if E_models is None:
        print("fit.getEs() returned Es=None")
    else:
        nan_total = int(np.isnan(E_models).sum())
        finite = E_models[np.isfinite(E_models)]
        if finite.size:
            print(f"E_models shape={E_models.shape} nan={nan_total} min={finite.min():.6f} max={finite.max():.6f}")
        else:
            print(f"E_models shape={E_models.shape} all values are NaN")
    print("sample_counts ", sample_counts, np.sum(sample_counts))

    # Write output
    istart = 0
    for i, f in enumerate(args.inputs):
        fname = os.path.join(args.inputs_dir, f)
        base = os.path.splitext(f)[0]  # Remove .xyz extension
        iend = istart + sample_counts[i]
        # Get grid mapping
        Gref, seq, axis, distances, angles = fit.parse_xyz_mapping(fname)
        print(f"parse_xyz_mapping[{base}] seq_len={len(seq)} grid_shape={Gref.shape} axis={axis}")
        Es_slice = E_models[istart:iend]
        # Map to grid with NaN padding
        Gmodel = np.empty_like(Gref); Gmodel[:] = np.nan
        nmap = min(len(Es_slice), len(seq))
        for k in range(nmap):
            idist, iang = seq[k]
            Gmodel[idist, iang] = Es_slice[k]
        if Es_slice.size:
            nan_slice = int(np.isnan(Es_slice).sum())
            finite_slice = Es_slice[np.isfinite(Es_slice)]
            if finite_slice.size:
                print(f"Es_slice[{base}] size={Es_slice.size} nan={nan_slice} min={finite_slice.min():.6f} max={finite_slice.max():.6f}")
            else:
                print(f"Es_slice[{base}] size={Es_slice.size} all NaN")
        else:
            print(f"Es_slice[{base}] is empty")
        finite_model = Gmodel[np.isfinite(Gmodel)]
        if finite_model.size:
            print(f"Gmodel[{base}] finite min={finite_model.min():.6f} max={finite_model.max():.6f}")
        else:
            print(f"Gmodel[{base}] has no finite values (shape={Gmodel.shape})")
        finite_ref = Gref[np.isfinite(Gref)]
        if finite_ref.size:
            print(f"Gref[{base}] finite min={finite_ref.min():.6f} max={finite_ref.max():.6f}")
        else:
            print(f"Gref[{base}] has no finite values (shape={Gref.shape})")
        data_base = os.path.join(args.out_dir, base)
        fit.save_data(Gref, Gmodel, angles, distances, save_data_prefix=data_base, save_fmt=args.save_fmt, kcal=bool(args.kcal), line=bool(args.line),)
        istart = iend

    # Optional, plotting
    if args.show == 1: plt.show()
