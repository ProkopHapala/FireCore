#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import argparse
import re
from pathlib import Path
import sys
# Ensure project root is on sys.path for `pyBall` imports when run as a script
_ROOT = Path(__file__).resolve().parents[2]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))
try:
    from pyBall.atomicUtils import scan_xyz
except Exception:
    scan_xyz = None

from pyBall.FitREQutils import (
    plot_list, plot_imshow, plot_polar, plot_profiles,
    read_scan_atomicutils, parse_xyz_blocks, compute_ra_vec,
    detect_rows_by_r, reshape_to_grid, compute_shift_from_grid,
)

if __name__ == '__main__':
    '''
    expected comment line is like:
    # n0 2 Etot -82594.964479 x0 02.16 z 161 HBr-A1_HCl-D1

    How to Run:
    with imshow:
        python split_scan_imshow.py --show --list ./HHalogens/toplot.txt --cmap bwr --sym --transpose
    polar:
        python split_scan_imshow.py --show --list ./HHalogens/toplot.txt --cmap bwr --sym  --polar --rmax 7.0 --half right --debug --transpose
    '''
    ap = argparse.ArgumentParser(description='Split packed 2D dimer scan and plot with imshow')
    ap.add_argument('--xyz',    type=str, help='input packed .xyz file', default="./HHalogens/HBr-A1_HBr-D1.xyz")
    ap.add_argument('--list',   type=str, help='panel list file (N M on first line, rows separated by ---)')
    ap.add_argument('--natoms', type=int, default=None, help='atoms per block (default: auto/infer)')
    ap.add_argument('--kcal',   action='store_true', help='plot energies in kcal/mol')
    ap.add_argument('--polar',  action='store_true', help='use polar plots (contourf) instead of rectangular imshow')
    ap.add_argument('--lines',  action='store_true', help='plot 1D profiles (radial at selected angles, angular at global r_min, and per-angle minima)')
    ap.add_argument('--rmax',   type=float, default=None, help='max radius for polar plots (clip outer radius to focus on short-range)')
    ap.add_argument('--half',   type=str, default='right', choices=['right','left'], help='half-circle to show in polar plots: right (-90..+90) or left (90..270)')
    ap.add_argument('--unsigned-angle', action='store_true', help='use unsigned angle in [0,180] deg (default uses signed angle)')
    ap.add_argument('--sym',    action='store_true', help='per-image symmetric color scale around E_far=0 (vmin=-|Emin|, vmax=+|Emin|)')
    ap.add_argument('--emin',   type=float, default=None, help='override lower bound for color scale (in displayed units); if only --emin is given (and no --emax), a symmetric scale [-emin, +emin] is used')
    ap.add_argument('--emax',   type=float, default=None, help='per-image range width above its global minimum: vmin=Emin_image, vmax=Emin_image+emax (in displayed units)')
    ap.add_argument('--cmap',   type=str, default='bwr', help='matplotlib colormap')
    ap.add_argument('--no-cbar',action='store_true', help='disable colorbar')
    ap.add_argument('--rtick-step', type=int, default=5, help='label every N-th radial sample on imshow y-axis (default: 5)')
    ap.add_argument('--transpose', action='store_true', help='swap rows and columns of the subplot grid')
    ap.add_argument('--show',   action='store_true', help='show plot (otherwise just save)')
    ap.add_argument('--out',    type=str, default=None, help='output image file (.png)')
    ap.add_argument('--debug',  action='store_true')
    args = ap.parse_args()

    if args.list:
        fig = plot_list(
            args.list,
            emin=args.emin,
            emax=args.emax,
            sym=args.sym,
            kcal=args.kcal,
            cmap=args.cmap,
            bColorbar=(not args.no_cbar),
            natoms=args.natoms,
            debug=args.debug,
            unsigned_angle=args.unsigned_angle,
            transpose=args.transpose,
            polar=args.polar,
            rmax=args.rmax,
            half=args.half,
            lines=args.lines,
            rtick_step=args.rtick_step,
        )
        if args.out is None:
            base = Path(args.list).with_suffix("")
            if args.lines:
                args.out = f"{base}_lines.png"
            elif args.polar:
                args.out = f"{base}_polar.png"
            else:
                args.out = f"{base}_grid.png"
        plt.savefig(args.out, dpi=160)
        if args.show:
            plt.show()
    else:
        # Single file path
        Es, Ps = read_scan_atomicutils(args.xyz)
        if Es.size == 0:
            print(f"WARNING: atomicUtils.scan_xyz() failed to parse {args.xyz} => fallback to local parse_xyz_blocks()")
            _, _, Ps_local = parse_xyz_blocks(args.xyz, natoms=args.natoms)
            Es = np.full((Ps_local.shape[0],), np.nan)
            Ps = Ps_local
        if args.debug:
            print(f"Parsed {len(Es)} blocks, natoms={Ps.shape[1]}")
        r, a = compute_ra_vec(Ps, signed=(not args.unsigned_angle))
        rows, step = detect_rows_by_r(r)
        if args.debug:
            print(f"Detected {len(rows)} rows; typical dr={step:.6f}")
            for k, (s, e) in enumerate(rows[:5]):
                rr = r[s:e]
                print(f" row {k}: n={e-s}, r[{s}:{e}] ~ {rr[0]:.3f}..{rr[-1]:.3f}, a~{np.nanmean(a[s:e]):.2f}")
        Vraw, R, A, rv = reshape_to_grid(Es, r, a, rows)
        shift = compute_shift_from_grid(Vraw)
        V = Vraw - shift
        title = Path(args.xyz).name
        if args.debug:
            try:
                rv_f = rv[np.isfinite(rv)] if rv is not None else np.array([])
                A_f  = A[np.isfinite(A)]   if A  is not None else np.array([])
                print(f"  rv[len={rv_f.size}] min,max = {np.nanmin(rv_f):.3f}, {np.nanmax(rv_f):.3f}   rv[:] = ", np.array2string(rv_f, precision=3))
                print(f"  A[len={A_f.size}] min,max  = {np.nanmin(A_f):.1f}, {np.nanmax(A_f):.1f}   A[:]  = ", np.array2string(A_f, precision=1))
            except Exception:
                pass
        # Compute vmin/vmax based on requested options (precedence: sym > emax > emin)
        fac = 23.060548 if args.kcal else 1.0
        if args.sym:
            emin_img = float(np.nanmin(V)) * fac
            vmag = abs(emin_img)
            vmin_plot = -vmag
            vmax_plot = +vmag
        elif args.emax is not None and args.emax > 0:
            vmin_plot = float(np.nanmin(V)) * fac
            vmax_plot = vmin_plot + args.emax
        else:
            vmin_plot = args.emin  # may be None
            vmax_plot = None
        if args.lines:
            plot_profiles(V, rv, A, R=R, rmax=args.rmax, kcal=args.kcal, title=title, vmin=vmin_plot, vmax=vmax_plot)
        elif args.polar:
            plot_polar(V, rv, A, emin=vmin_plot, vmax=vmax_plot, title=title, kcal=args.kcal, cmap=args.cmap, bColorbar=(not args.no_cbar), rmax=args.rmax, half=args.half, R=R)
        else:
            plot_imshow(V, rv, A, emin=vmin_plot, vmax=vmax_plot, title=title, kcal=args.kcal, cmap=args.cmap, bColorbar=(not args.no_cbar), rtick_step=args.rtick_step)
        if args.out is None:
            base = Path(args.xyz).with_suffix("")
            if args.lines:
                args.out = f"{base}_lines.png"
            elif args.polar:
                args.out = f"{base}_polar.png"
            else:
                args.out = f"{base}_imshow.png"
        plt.savefig(args.out, dpi=160)
        if args.show:
            plt.show()
