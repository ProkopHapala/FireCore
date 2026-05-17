#!/usr/bin/python3 -u
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
fit.plt = plt

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot and compare 2D energy surfaces from a single .xyz file or scan DOFs")
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
    # setScan
    parser.add_argument("--mode",       type=str,   default="scan",           help="Action: plot ref only, compare model (no fit), fit then compare, or scan DOFs", choices=["plot","model","fit","scan"])
    parser.add_argument("--scan_dofs",  type=int,   default=None,             help="List of DOF indices to scan. If None, all from dof-selection are scanned.",  nargs='+')
    parser.add_argument("--scan_range", type=float, default=[-1.0, 1.0, 100], help="Scan range: min max n_steps", nargs=3)
    # setOutput
    parser.add_argument("--scan_outfile",     type=str, default="DOF_scan", help="Base name for scan output data files")
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
        n_total = fit.loadXYZ(fname, bAddEpairs=(args.epairs == 1), bOutXYZ=False, bEvalOnlyCorrections=False, bAppend=bAppend)
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
    fit.setPenalty( Clamp=0, Regularize=args.regularize, AddRegError=args.regularize, RegCountWeight=args.regcountweight,
                    SoftClamp=args.softclamp, softClamp_start=args.softclamp_start, softClamp_max=args.softclamp_max )
    if args.user_weights:
        Erefs_all = np.concatenate(all_Erefs)
        x0s_all   = np.concatenate(all_x0s)
        weight_func = lambda E: fit.exp_weight_func(E, a=args.weight_a, alpha=args.weight_alpha, Emin0=args.emin0)
        weights0, lens = fit.split_and_weight_curves(Erefs_all, x0s_all,  n_before_min=args.n_before, weight_func=weight_func, EminMin=args.emin_min)
        fit.setWeights( weights0 )

    # --- Scanning
    DOFnames = fit.loadDOFnames(args.dof_selection)
    xs = np.linspace(args.scan_range[0], args.scan_range[1], int(args.scan_range[2])) 
    fit.plotDOFscans(args.scan_dofs, xs, DOFnames, title="DOF scan 1D", bFs=True, bEvalSamples=True, bPrint=False, outfile=args.scan_outfile)
