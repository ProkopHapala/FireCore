#!/usr/bin/env python3
"""
Plot 2D energy maps for multiple XYZ scan files, with
  row 1: 2D imshow + Rmin overlay + global min marker
  row 2: E_min(angle) + fixed-distance slice
  row 3: molecule geometry at the global minimum

Usage:
  python tests/tFitREQ/plot_ref.py -i scan.xyz --min-lines --kcal
  python tests/tFitREQ/plot_ref.py -i a.xyz b.xyz --dir /path/ --min-lines --kcal
"""
import argparse
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from pyBall.FitREQutils import (_build_frame_grid, plot_molecule, plot_system_panel, plot_energy_panel)
from pyBall.FitREQutils import parse_xyz_with_headers

'''

HOW TO USE IT:

python3 plot_ref.py -i CH2O-A1_H2O-D1-y.xyz CH2O-A1_H2O-D1-z.xyz --dir add_epairs_full --min-lines --kcal --ncols 2 -o test_plot.png

python3 plot_ref.py --dir add_epairs_full --min-lines --kcal --ncols 4 --nrow 2 -o CH2O_plots.png -i \
CH2O-A1_HF-D1-y.xyz  CH2O-A1_HF-D1-z.xyz  \
CH2O-A1_HCN-D1-y.xyz CH2O-A1_HCN-D1-z.xyz \
CH2O-A1_H2O-D1-y.xyz CH2O-A1_H2O-D1-z.xyz \
CH2O-A1_NH3-D1-y.xyz CH2O-A1_NH3-D1-z.xyz

python3 plot_ref.py --dir add_epairs_full --min-lines --kcal --ncols 4 --nrow 2 -o H2O_plots.png -i \
H2O-A1_HF-D1-y.xyz  H2O-A1_HF-D1-z.xyz  \
H2O-A1_HCN-D1-y.xyz H2O-A1_HCN-D1-z.xyz \
H2O-A1_H2O-D1-y.xyz H2O-A1_H2O-D1-z.xyz \
H2O-A1_NH3-D1-y.xyz H2O-A1_NH3-D1-z.xyz

python3 plot_ref.py --dir add_epairs_full --min-lines --kcal --ncols 4 --nrow 2 -o HCN_plots.png -i \
HCN-A1_HF-D1-y.xyz  HCN-A1_HF-D1-z.xyz  \
HCN-A1_HCN-D1-y.xyz HCN-A1_HCN-D1-z.xyz \
HCN-A1_H2O-D1-y.xyz HCN-A1_H2O-D1-z.xyz \
HCN-A1_NH3-D1-y.xyz HCN-A1_NH3-D1-z.xyz

python3 plot_ref.py --dir add_epairs_full --min-lines --kcal --ncols 4 --nrow 2 -o NH3_plots.png -i \
NH3-A1_HF-D1-y.xyz  NH3-A1_HF-D1-z.xyz  \
NH3-A1_HCN-D1-y.xyz NH3-A1_HCN-D1-z.xyz \
NH3-A1_H2O-D1-y.xyz NH3-A1_H2O-D1-z.xyz \
NH3-A1_NH3-D1-y.xyz NH3-A1_NH3-D1-z.xyz

'''

def main():
    ap = argparse.ArgumentParser(  description="Plot 2D energy maps + min curves + molecule geometry")
    ap.add_argument("-i", "--input", nargs="+", help="XYZ file(s)")
    ap.add_argument("--dir", default=None, help="Directory of XYZ files")
    ap.add_argument("--glob", default=None, help="Glob pattern in --dir")
    ap.add_argument("-o", "--output", default=None, help="Save figure (PNG)")
    ap.add_argument("--kcal", action="store_true", help="Plot energies in kcal/mol")
    ap.add_argument("--sym", action="store_true",  help="Symmetric color scale")
    ap.add_argument("--min-lines", action="store_true",help="Overlay R_min + energy panel + molecule geometry")
    ap.add_argument("--ncols", type=int, default=None, help="Number of columns")
    ap.add_argument("--nrow", type=int, default=None, help="Number of rows (files stacked per column)")
    ap.add_argument("--debug", action="store_true")
    args = ap.parse_args()

    if not args.input and not args.glob:
        ap.print_help()
        sys.exit(1)

    paths = []
    if args.dir and args.glob:
        import glob
        paths = sorted(glob.glob(os.path.join(args.dir, args.glob)))
    elif args.dir:
        for name in args.input:
            p = os.path.join(args.dir, name)
            if not p.endswith(".xyz"): p += ".xyz"
            if os.path.isfile(p): paths.append(p)
            else: print("Warning: not found:", p)
    else:
        for name in args.input:
            if os.path.isfile(name): paths.append(name)
            else: print("Warning: not found:", name)
    if not paths:
        print("No input files found.")
        sys.exit(1)

    nfiles = len(paths)
    nrows_sys = 3 if args.min_lines else 1
    
    # Column-major layout if --nrow is specified
    if args.nrow:
        nrow = args.nrow
        ncols = args.ncols or ((nfiles + nrow - 1) // nrow)
        nr = nrows_sys * nrow
        d2 = nrow
    else:
        ncols = args.ncols or min(nfiles, 4)
        nr = nrows_sys * ((nfiles + ncols - 1) // ncols)
        d2 = (nfiles + ncols - 1) // ncols
    fig_height = {1: 3.5, 3: 3.5 + 2.6 * 2}[nrows_sys] * d2
    fig, axs = plt.subplots(nr, ncols, figsize=(4.0 * ncols, fig_height),squeeze=False)
    use_sym = args.sym or args.min_lines

    print("Plotting %d system(s) ..." % nfiles)
    for i, fp in enumerate(paths):
        if args.nrow:
            # Column-major: stack files vertically in each column
            c = i // nrow
            r0 = (i % nrow) * nrows_sys
        else:
            # Row-major: traditional left-to-right, top-to-bottom
            r0 = (i // ncols) * nrows_sys
            c = i % ncols
        ax_img = axs[r0, c]
        label = os.path.splitext(os.path.basename(fp))[0]

        try:
            V, rv, A, shift, frame_idx, Ps_raw, Ts_raw = _build_frame_grid(fp)
        except Exception as e:
            print("  ERROR: %s: %s" % (fp, e))
            ax_img.set_axis_off()
            if args.min_lines:
                axs[r0 + 1, c].set_axis_off()
                axs[r0 + 2, c].set_axis_off()
            continue

        if not np.any(np.isfinite(V)):
            print("  WARNING: all-NaN: %s" % label)
            ax_img.set_axis_off()
            if args.min_lines:
                axs[r0 + 1, c].set_axis_off()
                axs[r0 + 2, c].set_axis_off()
            continue

        glob_min = plot_system_panel(V, rv, A, ax_img, label, args.kcal,use_sym, overlay_rmin=args.min_lines)

        if args.min_lines:
            ax_en = axs[r0 + 1, c]
            plot_energy_panel(V, rv, A, ax_en, args.kcal)

            ax_mol = axs[r0 + 2, c]
            if glob_min is not None:
                iy_g, ix_g = glob_min
                iframe = frame_idx[iy_g, ix_g]
                if iframe >= 0 and iframe < len(Ps_raw):
                    apos_mol = Ps_raw[iframe]
                    enames_mol = (Ts_raw[iframe] if Ts_raw is not None and iframe < len(Ts_raw)else ["?"] * len(apos_mol))
                    n0_mol = None
                    _, _, _, N0s = parse_xyz_with_headers(fp)
                    if iframe < len(N0s) and N0s[iframe] > 0:n0_mol = int(N0s[iframe])
                    plot_molecule(ax_mol, enames_mol, apos_mol, n0=n0_mol, title="Frame %d" % iframe)
                else:
                    ax_mol.set_axis_off()
            else:
                ax_mol.set_axis_off()

    for i in range(nfiles * nrows_sys, nr * ncols):
        r, c = i // ncols, i % ncols
        axs[r, c].set_visible(False)
    fig.tight_layout()

    if args.output:
        fig.savefig(args.output, dpi=200, bbox_inches='tight')
        print("Saved:", args.output)
    else:
        plt.show()


if __name__ == "__main__":
    main()
