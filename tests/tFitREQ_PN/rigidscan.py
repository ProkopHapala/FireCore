#!/usr/bin/python3 -u

import sys
import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt
import tkinter
import time
import argparse
sys.path.append("/home/niko/work/HBOND/FireCore/")
from pyBall import FitREQ_PN as fit

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fit model parameters to reproduce 2D energy surfaces from multiple .xyz files")
    # setVerbosity
    parser.add_argument("--verbosity",          type=int, default=2, help="Verbosity for FitREQ")
    parser.add_argument("--idebug",             type=int, default=0, help="Flag for debug")
    parser.add_argument("--printdofs",          type=int, default=1, help="Set PrintDOFs")
    parser.add_argument("--printfdofs",         type=int, default=1, help="Set PrintfDOFs")
    parser.add_argument("--printbeforreg",      type=int, default=1, help="Set PrintBeforReg")
    parser.add_argument("--printafterreg",      type=int, default=1, help="Set PrintAfterReg")
    parser.add_argument("--printoverrepulsive", type=int, default=1, help="Set PrintOverRepulsive")
    # setModel
    parser.add_argument("--ivdw",    type=int,   default=1,   help="Set ivdW (0=no vdW, 1=LJ, 2=LJr8, 3=LJr9, 4=Morse, 5=Buck")
    parser.add_argument("--icoul",   type=int,   default=1,   help="Set iCoul (0=no Coul, 1=point charges, 2=soft clamping, 10-14=Boys clamping, 10=exact erf/r, 11=cubic C1, 12=quintic C2, 13=quartic even C1, 14=sextic even C2)")
    parser.add_argument("--ihbond",  type=int,   default=0,   help="Set iHBond (0=no HBond correction, 1=H1 correction, 2=H2 correction, 3=H1 and H2 correction")
    parser.add_argument("--epairs",  type=int,   default=1,   help="Enable/disable Epair terms (it applies both when loading XYZ and when calculating interactions)")
    parser.add_argument("--iepairs", type=int,   default=0,   help="Set iEpairs (0=no interaction, 1=SR interaction, 2=SR2 interaction)")
    parser.add_argument("--lepairs", type=float, default=1.0, help="Set distance of Epairs from host atom")
    parser.add_argument("--bpn",     type=int,   default=1,   help="Set boolean for PN functions")
    parser.add_argument("--kmorse",  type=float, default=1.8, help="Set curvature of Morse potential")
    # setInput
    parser.add_argument("--inputs",        type=str, default=None,   help="List of input .xyz files (single movie). If empty, uses default set", nargs='*')
    parser.add_argument("--inputs_dir",    type=str, default='data', help="Directory where input files are located")
    parser.add_argument("--dof_selection", type=str, default=None,   help="DOF selection file")
    # setPenalty
    parser.add_argument("--clamp",           type=int,   default=1,     help="Hardly restrain the values of parameters during optimization")
    parser.add_argument("--regularize",      type=int,   default=0,     help="Apply harmonic potentials to keep parameters within a specified range")
    parser.add_argument("--regcountweight",  type=int,   default=0,     help="Apply further normalization on regularizing potentials")
    parser.add_argument("--softclamp",       type=int,   default=0,     help="Apply soft clamping to specific samples")
    parser.add_argument("--softclamp_start", type=float, default=4.0,   help="Minimum threshold where soft clamping starts to kick in")
    parser.add_argument("--softclamp_max",   type=float, default=6.0,   help="Maximum value of the penalty function")
    parser.add_argument("--user_weights",    type=int,   default=1,     help="Enable user weights during scan")
    parser.add_argument("--n_before",        type=int,   default=1000,  help="#points before min to weight")
    parser.add_argument("--weight_a",        type=float, default=1.0,   help="Weight amplitude 'a' for exp weight func")
    parser.add_argument("--weight_alpha",    type=float, default=4.0,   help="Weight sharpness 'alpha' for exp weight func")
    parser.add_argument("--emin_min",        type=float, default=-0.02, help="Emin threshold for weighting segments")
    parser.add_argument("--emin0",           type=float, default=0.1,   help="Adds to denominator of the exponential function")
    # setOptimization
    parser.add_argument("--nstep", type=int, default=100, help="Fitting steps")
    # setOutput
    parser.add_argument("--out_dir",            type=str, default=None,      help="Directory to save all output data (2D maps, 1D lines)")
    parser.add_argument("--save_fmt",           type=str, default="gnuplot", help="Format for saved data", choices=["both","npz","gnuplot"])
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
    fit.setModel(ivdW=args.ivdw, iCoul=args.icoul, iHbond=args.ihbond, Epairs=args.epairs, iEpairs=args.iepairs, kMorse=args.kmorse, Lepairs=args.lepairs, bPN=(args.bpn == 1))

    # setInput
    fit.loadTypes(fEtypes="data/ElementTypes.dat", fAtypes="data/AtomTypes.dat")
    fit.loadDOFSelection(fname=args.dof_selection)
    sample_counts = []
    total_loaded = 0
    all_Erefs = []
    all_x0s = []
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
            Erefs, x0s = fit.read_xyz_data(fname)
            all_Erefs.append(Erefs)
            all_x0s.append(x0s)
    fit.getBuffs()
    
    # setPenalty
    fit.setPenalty( Clamp=args.clamp, Regularize=args.regularize, AddRegError=args.regularize, RegCountWeight=args.regcountweight,
                    SoftClamp=args.softclamp, softClamp_start=args.softclamp_start, softClamp_max=args.softclamp_max )
    if args.user_weights:
        Erefs_all = np.concatenate(all_Erefs)
        x0s_all   = np.concatenate(all_x0s)
        weight_func = lambda E: fit.exp_weight_func(E, a=args.weight_a, alpha=args.weight_alpha, Emin0=args.emin0)
        weights0, lens = fit.split_and_weight_curves(Erefs_all, x0s_all,  n_before_min=args.n_before, weight_func=weight_func, EminMin=args.emin_min)
        fit.setWeights( weights0 )

    # --- Fitting
    trj_E, trj_F, trj_DOFs, _ = fit.setTrjBuffs(niter=args.nstep)
    DOFnames = fit.loadDOFnames(args.dof_selection)
    
    # Create output directoryies if needed
    if not os.path.exists(args.out_dir): os.makedirs(args.out_dir)
    
    Err = fit.getError( iparallel=0 )
    print(Err)
    quit()

    arrays = {"E_tot": E_tot, "E_coul": E_coul, "E_vdw": E_vdw, "E_hcorr": E_hcorr, "E_epairs": E_epairs}

    def check_array(arr, name):
        if arr is None:
            print(f"fit.getEs_components() returned {name}=None")
        else:
            nan_total = int(np.isnan(arr).sum())
            finite = arr[np.isfinite(arr)]
            if finite.size:
                print(f"{name} shape={arr.shape} nan={nan_total} min={finite.min():.6f} max={finite.max():.6f}")
            else:
                print(f"{name} shape={arr.shape} all values are NaN")

    def report_slice(arr, name, base):
        if arr.size:
            nan_slice = int(np.isnan(arr).sum())
            finite = arr[np.isfinite(arr)]
            if finite.size:
                print(f"{name}[{base}] size={arr.size} nan={nan_slice} min={finite.min():.6f} max={finite.max():.6f}")
            else:
                print(f"{name}[{base}] size={arr.size} all NaN")
        else:
            print(f"{name}[{base}] is empty")

    def report_finite(G, name, base):
        finite = G[np.isfinite(G)]
        if finite.size:
            print(f"{name}[{base}] finite min={finite.min():.6f} max={finite.max():.6f}")
        else:
            print(f"{name}[{base}] has no finite values (shape={G.shape})")
            
    for name, arr in arrays.items():
        check_array(arr, name)

    print("sample_counts ", sample_counts, np.sum(sample_counts))
    
    istart = 0
    for i, f in enumerate(args.inputs):
        fname = os.path.join(args.inputs_dir, f)
        base = os.path.splitext(f)[0]  # Remove .xyz extension
        #title = base.replace('_', ' ')  # Clean title
        iend = istart + sample_counts[i]
        
        # Get grid mapping
        Gref, seq, axis, distances, angles = fit.parse_xyz_mapping(fname)
        print(f"parse_xyz_mapping[{base}] seq_len={len(seq)} grid_shape={Gref.shape} axis={axis}")

        E_comp_slice = {}
        G_comp = {}
        for name, arr in arrays.items():
            E_comp_slice[name] = arr[istart:iend]
            G_comp[name] = np.full_like(Gref, np.nan)
            nmap = min(len(E_comp_slice[name]), len(seq))
            G = G_comp[name]
            for k in range(nmap):
                idist, iang = seq[k]
                G[idist, iang] = E_comp_slice[name][k]
        for name, arr in E_comp_slice.items():
            report_slice(arr, f"{name}_slice", base)
        for name, G in G_comp.items():
            report_finite(G, name, base)        
        report_finite(Gref, "Gref", base)

        G_comp_coul_vdw = np.full_like(Gref, np.nan)
        mask = np.isfinite(G_comp["E_coul"]) & np.isfinite(G_comp["E_vdw"])
        G_comp_coul_vdw[mask] = G_comp["E_coul"][mask] + G_comp["E_vdw"][mask]

        G_comp_hcorr_epairs = np.full_like(Gref, np.nan)
        mask = np.isfinite(G_comp["E_hcorr"]) & np.isfinite(G_comp["E_epairs"])
        G_comp_hcorr_epairs[mask] = G_comp["E_hcorr"][mask] + G_comp["E_epairs"][mask]

        Gref_nocoul_novdw = np.full_like(Gref, np.nan)
        mask = np.isfinite(Gref) & np.isfinite(G_comp_coul_vdw)
        Gref_nocoul_novdw[mask] = Gref[mask] - G_comp_coul_vdw[mask]

        Gdiff_nocoul_novdw = np.full_like(Gref, np.nan)
        mask = np.isfinite(Gref_nocoul_novdw) & np.isfinite(G_comp_hcorr_epairs)
        Gdiff_nocoul_novdw[mask] = Gref_nocoul_novdw[mask] - G_comp_hcorr_epairs[mask]

        Gdiff_tot = np.full_like(Gref, np.nan)
        mask = np.isfinite(Gref) & np.isfinite(G_comp["E_tot"])
        Gdiff_tot[mask] = Gref[mask] - G_comp["E_tot"][mask]
              
        grids = {"G_ref": Gref, "G_ref_nocoul_novdw": Gref_nocoul_novdw, "G_comp_tot": G_comp["E_tot"], "G_comp_coul": G_comp["E_coul"], "G_comp_vdw": G_comp["E_vdw"], "G_comp_hcorr": G_comp["E_hcorr"], "G_comp_epairs": G_comp["E_epairs"], "G_comp_coul_vdw": G_comp_coul_vdw, "G_comp_hcorr_epairs": G_comp_hcorr_epairs, "G_diff_nocoul_novdw": Gdiff_nocoul_novdw, "G_diff_tot": Gdiff_tot}

        data_base = os.path.join(args.out_dir, base)
        fit.save_data(Gref, G_comp["E_tot"], angles, distances, save_data_prefix=data_base, save_fmt=args.save_fmt, kcal=bool(args.kcal), line=bool(args.line),)
        for name, G in grids.items():
            fit.save_data_single(G, angles, distances, save_data_prefix=data_base, save_fmt=args.save_fmt, kcal=bool(args.kcal), suffix=name)
        istart = iend

