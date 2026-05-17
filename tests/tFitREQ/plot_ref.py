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
from matplotlib import collections as mc

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from pyBall import elements
from pyBall.FitREQutils import (compute_panel_data, plot_imshow,
                                extract_min_curves, parse_xyz_blocks,
                                read_scan_atomicutils, compute_ra_vec,
                                detect_rows_by_r, reshape_to_grid,
                                compute_shift_from_grid)
from pyBall.atomicUtils import findBondsNP


def _build_frame_grid(fp):
    """Parse XYZ, build V grid + frame-index grid + atom data."""
    Es, Ps = read_scan_atomicutils(fp)
    Ts = None
    if Es.size == 0:
        Es, Ts, Ps = parse_xyz_blocks(fp)
    else:
        # read_scan_atomicutils works but doesn't give types — get them separately
        _, Ts2, _ = parse_xyz_blocks(fp)
        Ts = Ts2 if len(Ts2) == len(Ps) else None

    r, a = compute_ra_vec(Ps, signed=True)
    rows, _ = detect_rows_by_r(r)
    ny = len(rows)
    nx = max(e - s for s, e in rows)
    frame_idx = np.full((ny, nx), -1, dtype=int)
    for iy, (s, e) in enumerate(rows):
        n = e - s
        frame_idx[iy, :n] = np.arange(s, e)

    Vraw, R, A, rv = reshape_to_grid(Es, r, a, rows)
    shift = compute_shift_from_grid(Vraw)
    V = Vraw - shift
    return V, rv, A, shift, frame_idx, Ps, Ts


def _element_color(elem):
    """Hex color for an element; 'E' (dummy) -> purple."""
    if elem == 'E':
        return '#9932CC'  # purple
    try:
        return elements.ELEMENT_DICT[elem][8]
    except (KeyError, IndexError):
        return '#808080'


def _element_size(elem):
    """Marker size for an element; 'E' (dummy) -> small."""
    if elem == 'E':
        return 25
    try:
        return elements.ELEMENT_DICT[elem][6] * 30
    except (KeyError, IndexError):
        return 30


def plot_molecule(ax, enames, apos, n0=None, title=""):
    """Draw molecular geometry on ax.
    
    Colors by element (E=purple). Labels = 1-based atom indices.
    n0 splits donor/acceptor groups (labels shown per fragment).
    """
    if enames is None or apos is None or len(apos) == 0:
        ax.text(0.5, 0.5, "No geometry", ha='center', va='center',
                transform=ax.transAxes)
        ax.set_axis_off()
        return

    elem = [e.split('_')[0] for e in enames]
    nat = len(apos)
    if nat == 0:
        ax.set_axis_off()
        return

    # Colors and sizes by element (E = electron pair = purple)
    colors = [_element_color(e) for e in elem]
    sizes  = [_element_size(e) for e in elem]

    ax.scatter(apos[:, 0], apos[:, 1], c=colors, s=sizes, zorder=5,
               edgecolors='black', linewidths=0.5)

    # Labels: 1-based index, numbered per fragment if n0 given
    if n0 is not None and 0 < n0 < nat:
        for i in range(n0):
            ax.annotate(str(i + 1), apos[i, :2], textcoords="offset points",
                        xytext=(0, 6), fontsize=6, ha='center', color='blue',
                        fontweight='bold')
        for i in range(n0, nat):
            ax.annotate(str(i + 1 - n0) + "'", apos[i, :2],
                        textcoords="offset points", xytext=(0, 6),
                        fontsize=6, ha='center', color='red',
                        fontweight='bold')
        # Separator line
        mid_x = (apos[:n0, 0].max() + apos[n0:, 0].min()) / 2
        ax.axvline(x=mid_x, color='gray', ls=':', lw=1, alpha=0.5)
    else:
        for i in range(nat):
            ax.annotate(str(i + 1), apos[i, :2], textcoords="offset points",
                        xytext=(0, 6), fontsize=6, ha='center',
                        color='black', fontweight='bold')

    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])
    if title:
        ax.set_title(title, fontsize=8)
    if nat > 0:
        ptp = np.ptp(apos[:, :2], axis=0)
        half = max(ptp.max() / 2, 1.5) + 0.7
        center = apos[:, :2].mean(axis=0)
        ax.set_xlim(center[0] - half, center[0] + half)
        ax.set_ylim(center[1] - half, center[1] + half)


def plot_system_panel(V, rv, A, ax, label, kcal, sym, overlay_rmin=False):
    """Row 1: 2D imshow + Rmin overlay + global min marker."""
    if not np.any(np.isfinite(V)):
        ax.set_axis_off()
        return None

    im = plot_imshow(V, rv, A, title=label, kcal=kcal, ax=ax,
                     bColorbar=True, rtick_step=5, bSym=False)
    if sym and im is not None:
        conv = 23.060548 if kcal else 1.0
        vmin_ref = np.nanmin(V) * conv
        if vmin_ref < 0:
            im.set_clim(vmin_ref, -vmin_ref)
        else:
            vabs = np.nanmax(np.abs(V)) * conv
            im.set_clim(-vabs, vabs)
        cb = im.colorbar
        if cb is not None:
            cb.update_normal(im)

    glob_min = None
    if np.any(np.isfinite(V)):
        iy_g, ix_g = np.unravel_index(np.nanargmin(V), V.shape)
        glob_min = (iy_g, ix_g)
        yr = rv[np.isfinite(rv)]
        if yr.size >= 2 and len(A) > iy_g:
            y0, y1 = np.nanmin(yr), np.nanmax(yr)
            if y1 > y0:
                y_mark = y0 + (ix_g + 0.5) * (y1 - y0) / V.shape[1]
                ax.plot(A[iy_g], y_mark, '+', color='white', ms=12, mew=2)

    if overlay_rmin:
        ny, nx = V.shape
        pix = np.full(ny, np.nan)
        for iy in range(ny):
            if np.isfinite(V[iy, :]).any():
                pix[iy] = np.nanargmin(V[iy, :])
        o = np.argsort(A)
        As, ps = A[o], pix[o]
        if np.isfinite(ps).any() and np.nanmin(As) != np.nanmax(As):
            yr = rv[np.isfinite(rv)]
            if yr.size >= 2:
                y0, y1 = np.nanmin(yr), np.nanmax(yr)
                if y1 > y0:
                    ax.plot(As, y0 + (ps + 0.5) * (y1 - y0) / nx,
                            'k-', lw=1.5, alpha=0.8)
    return glob_min


def plot_energy_panel(V, rv, A, ax, kcal):
    """Row 2: E_min(angle) + fixed-r slice."""
    conv = 23.060548 if kcal else 1.0
    unit = "kcal/mol" if kcal else "eV"
    _, emin = extract_min_curves(A, rv, V.T)
    o = np.argsort(A)
    ax.plot(A[o], emin[o] * conv, 'k-', lw=1.5, label="E$_{min}$")
    if np.any(np.isfinite(V)):
        iy_g, ix_g = np.unravel_index(np.nanargmin(V), V.shape)
        r_glob = rv[ix_g] if ix_g < len(rv) else np.nan
        sl = V[:, ix_g] * conv
        so = np.argsort(A)
        lbl = "E @ r=%.2f" % r_glob if np.isfinite(r_glob) else ""
        ax.plot(A[so], sl[so], 'r--', lw=1.5, label=lbl)
    ax.set_xlabel("Angle (deg)")
    ax.set_ylabel("E (%s)" % unit)
    ax.legend(fontsize=7)
    ax.grid(alpha=0.2)


def main():
    ap = argparse.ArgumentParser(
        description="Plot 2D energy maps + min curves + molecule geometry")
    ap.add_argument("-i", "--input", nargs="+", help="XYZ file(s)")
    ap.add_argument("--dir", default=None, help="Directory of XYZ files")
    ap.add_argument("--glob", default=None, help="Glob pattern in --dir")
    ap.add_argument("-o", "--output", default=None, help="Save figure (PNG)")
    ap.add_argument("--kcal", action="store_true",
                    help="Plot energies in kcal/mol")
    ap.add_argument("--sym", action="store_true",
                    help="Symmetric color scale")
    ap.add_argument("--min-lines", action="store_true",
                    help="Overlay R_min + energy panel + molecule geometry")
    ap.add_argument("--ncols", type=int, default=None,
                    help="Number of columns")
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
    ncols = args.ncols or min(nfiles, 4)
    nrows_sys = 3 if args.min_lines else 1
    nr = nrows_sys * ((nfiles + ncols - 1) // ncols)
    d2 = (nfiles + ncols - 1) // ncols
    fig_height = {1: 3.5, 3: 3.5 + 2.6 * 2}[nrows_sys] * d2
    fig, axs = plt.subplots(nr, ncols, figsize=(4.0 * ncols, fig_height),
                            squeeze=False)
    use_sym = args.sym or args.min_lines

    print("Plotting %d system(s) ..." % nfiles)
    for idx, fp in enumerate(paths):
        r0 = (idx // ncols) * nrows_sys
        c = idx % ncols
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

        glob_min = plot_system_panel(V, rv, A, ax_img, label, args.kcal,
                                     use_sym, overlay_rmin=args.min_lines)

        if args.min_lines:
            ax_en = axs[r0 + 1, c]
            plot_energy_panel(V, rv, A, ax_en, args.kcal)

            ax_mol = axs[r0 + 2, c]
            if glob_min is not None:
                iy_g, ix_g = glob_min
                iframe = frame_idx[iy_g, ix_g]
                if iframe >= 0 and iframe < len(Ps_raw):
                    apos_mol = Ps_raw[iframe]
                    enames_mol = (Ts_raw[iframe] if Ts_raw is not None
                                  and iframe < len(Ts_raw)
                                  else ["?"] * len(apos_mol))
                    n0_mol = None
                    try:
                        from pyBall.FitREQutils import parse_xyz_with_headers
                        _, _, _, N0s = parse_xyz_with_headers(fp)
                        if iframe < len(N0s) and N0s[iframe] > 0:
                            n0_mol = int(N0s[iframe])
                    except Exception:
                        pass
                    plot_molecule(ax_mol, enames_mol, apos_mol,
                                  n0=n0_mol, title="Frame %d" % iframe)
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
