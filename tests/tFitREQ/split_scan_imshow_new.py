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

# ----------------------------
# Parsing packed XYZ scan file
# ----------------------------

def parse_xyz_blocks(fname, natoms=None):
    """
    Parse a packed XYZ where blocks appear as:
      [optional natoms line] \n
      # E = <energy> eV \n
      <natoms lines of 'Sym x y z'>\n
    Some files repeat a single integer line (natoms) between blocks. We ignore
    such lines and rely on the comment + natoms atom lines pattern.

    Returns arrays: Es [n], types [n, natoms], pos [n, natoms, 3]
    """
    Es   = []
    Ts   = []
    Ps   = []

    # Support both old '# E = ...' and new 'Etot ...' formats
    energy_re = re.compile(r"(?:#\s*E\s*=\s*|\bEtot\s+)([+-]?(?:\d*\.\d+|\d+)(?:[Ee][+-]?\d+)?)")
    with open(fname, 'r') as f: lines = f.readlines()
    i = 0
    nline = len(lines)

    # If natoms not specified, try to infer from the first block
    def try_peek_natoms(k):
        # Find the next non-empty, non-comment line that is integer
        while k < nline and lines[k].strip() == "": k += 1
        if k < nline and lines[k].strip().isdigit(): return int(lines[k].strip()), k+1
        return None, k

    # Advance to first comment line with energy
    while i < nline and not lines[i].lstrip().startswith('#'): i += 1

    # Try infer natoms if not provided
    if natoms is None:
        # Look back one line for natoms if possible, otherwise forward search
        n_guess = None
        if i-1 >= 0 and lines[i-1].strip().isdigit():
            n_guess = int(lines[i-1].strip())
        else:
            n_guess, _ = try_peek_natoms(i+1)
        if n_guess is None:  n_guess = 4  # Fallback: expect 4 atoms per the H–halogen dimers case
        natoms = n_guess

    while i < nline:
        s = lines[i].lstrip()
        if not s.startswith('#'):
            i += 1
            continue
        m = energy_re.search(s)
        if not m:
            i += 1
            continue
        E = float(m.group(1))
        i += 1
        # Skip any natoms integer lines or empty lines between comment and atoms
        taken = 0
        types = []
        pos   = []
        while i < nline and taken < natoms:
            t = lines[i].strip()
            if t == "" or t.isdigit():
                i += 1
                continue
            parts = t.split()
            if len(parts) < 4:
                # Not an atom line; skip
                i += 1
                continue
            sym = parts[0]
            try:
                x, y, z = map(float, parts[1:4])
            except ValueError:
                # Malformed numeric line, skip
                i += 1
                continue
            types.append(sym)
            pos.append((x, y, z))
            taken += 1
            i += 1
        if taken == natoms:
            Es.append(E)
            Ts.append(types)
            Ps.append(pos)
        # else: incomplete block at EOF; stop

        # Optionally skip trailing natoms line after atom block
        while i < nline and (lines[i].strip() == "" or lines[i].strip().isdigit()):
            i += 1

    Es = np.array(Es)
    Ps = np.array(Ps, dtype=float)  # shape [n, natoms, 3]
    return Es, Ts, Ps

def read_scan_atomicutils(fname):
    """
    Read a packed XYZ using pyBall.atomicUtils.scan_xyz() and extract
    energies from the comment lines. Returns (Es [n], Ps [n,nat,3]).
    """
    # Support both old '# E = ...' and new 'Etot ...' formats
    energy_re = re.compile(r"(?:#\s*E\s*=\s*|\bEtot\s+)([+-]?(?:\d*\.\d+|\d+)(?:[Ee][+-]?\d+)?)")

    def _cb(block, id=None, comment=None):
        apos, _es = block  # _es is element names list from load_xyz
        E = np.nan
        if comment:
            m = energy_re.search(comment)
            if m:
                E = float(m.group(1))
        return (E, np.array(apos, dtype=float))

    if scan_xyz is None:
        return np.empty((0,)), np.empty((0, 0, 3))
    res = scan_xyz(fname, callback=_cb)
    if len(res) == 0:
        return np.empty((0,)), np.empty((0, 0, 3))
    Es, Plist = zip(*res)
    Es = np.array(Es, dtype=float)
    nat = Plist[0].shape[0]
    Ps = np.stack(Plist, axis=0).reshape(len(Plist), nat, 3)
    return Es, Ps

# ----------------------------
# Geometry-derived scan params
# ----------------------------

def derive_ra_from_block(P):
    """
    Given a 4-atom H–halogen dimer block with ordering:
      mol A: atoms 0,1
      mol B: atoms 2,3
    compute:
      r = |COM_B - COM_A|
      a = angle(axis_B, COM_B - COM_A) in degrees
    This should be stable across angles and dimers. If natoms != 4, it still
    uses first half vs second half atoms as the two fragments.
    """
    nat = P.shape[0]
    h = nat // 2
    A = P[:h, :]
    B = P[h:, :]
    cA = A.mean(axis=0)
    cB = B.mean(axis=0)
    R  = cB - cA
    r  = np.linalg.norm(R)
    if r < 1e-12:
        a = np.nan
    else:
        # Axis of B: use first two atoms of B (index h and h+1); if more, use PCA? here simple vector
        vB = B[-1] - B[0]
        nb = np.linalg.norm(vB)
        if nb < 1e-12:
            a = np.nan
        else:
            cosang = np.dot(vB, R) / (nb * r)
            cosang = np.clip(cosang, -1.0, 1.0)
            a = np.degrees(np.arccos(cosang))
    return r, a


def compute_ra_vec(P, signed=True):
    """Compute r and a from the moving molecule's first atom in the XZ plane.
    Convention:
    - Pivot is the first atom of molecule A (index 0).
    - Moving point is the first atom of molecule B (index h).
    - r = sqrt((x_B0-x_A0)^2 + (z_B0-z_A0)^2)
    - a = atan2(z_B0-x_A0, x_B0-x_A0) in degrees; signed in [-180,180], or |a| if unsigned.
    """
    n, nat, _ = P.shape
    h = nat // 2
    A0 = P[:, 0, :]
    B0 = P[:, h + 0, :]
    V  = B0 - A0
    x = V[:, 0]
    z = V[:, 2]
    r = np.sqrt(x*x + z*z)
    a = np.degrees(np.arctan2(z, x))
    if not signed:
        a = np.abs(a)
    return r, a

# ----------------------------
# Header-derived r/a (optional)
# ----------------------------

def parse_headers_ra(fname):
    """
    Parse comment headers to extract per-block Etot (energy), radius (x0) and angle (z).
    Supports both old energy-only headers and new enriched headers.
    Returns (Eh, Rh, Ah) as float arrays; values may be NaN if missing.
    """
    Eh = []
    Rh = []
    Ah = []
    reE = re.compile(r"(?:#\s*E\s*=\s*|\bEtot\s+)([+-]?(?:\d*\.\d+|\d+)(?:[Ee][+-]?\d+)?)")
    reR = re.compile(r"\bx0\s+([+-]?(?:\d*\.\d+|\d+)(?:[Ee][+-]?\d+)?)")
    reA = re.compile(r"\bz\s+([+-]?(?:\d*\.\d+|\d+)(?:[Ee][+-]?\d+)?)")
    with open(fname, 'r') as f:
        for ln in f:
            s = ln.lstrip()
            if not s.startswith('#'):
                continue
            mE = reE.search(s)
            mR = reR.search(s)
            mA = reA.search(s)
            Eh.append(float(mE.group(1)) if mE else np.nan)
            Rh.append(float(mR.group(1)) if mR else np.nan)
            Ah.append(float(mA.group(1)) if mA else np.nan)
    return np.array(Eh, dtype=float), np.array(Rh, dtype=float), np.array(Ah, dtype=float)

# ----------------------------
# Row detection and reshaping
# ----------------------------

def detect_rows_by_r(r, tol=None):
    """
    Find indices where a new angular scanline starts. The faster index is r.
    We detect row starts where dr = r[i] - r[i-1] is a negative jump larger
    in magnitude than half the typical step.
    Returns list of (start, end) slices covering the whole sequence.
    """
    dr = np.diff(r)
    pos = dr[dr > 0]
    if len(pos) == 0:
        step = np.median(np.abs(dr)) if len(dr) > 0 else 1.0
    else:
        step = np.median(pos)
    if tol is None:
        tol = 0.5 * step
    # Row breaks where dr < -tol
    brk = np.where(dr < -tol)[0] + 1
    splits = np.concatenate(([0], brk, [len(r)]))
    rows = [(int(splits[i]), int(splits[i+1])) for i in range(len(splits)-1)]
    return rows, step


def reshape_to_grid(vals, r, a, rows):
    """
    Pad rows to a rectangular grid. Returns:
      V [ny, nx], R [ny, nx], A [ny] (row-mean angle), rv [nx] (first row r)
    Missing points are filled with NaN.
    """
    ny = len(rows)
    nx = max(e - s for s, e in rows)
    V = np.full((ny, nx), np.nan)
    R = np.full((ny, nx), np.nan)
    A = np.full((ny,), np.nan)
    rv = np.full((nx,), np.nan)
    for iy, (s, e) in enumerate(rows):
        n = e - s
        V[iy, :n] = vals[s:e]
        R[iy, :n] = r[s:e]
        # Use angle measured at the maximum radius within the row to minimize
        # sensitivity to atom-picking errors. Fallback to mean if needed.
        try:
            if n > 0:
                irel = int(np.nanargmax(r[s:e]))
                A[iy] = a[s + irel]
            else:
                A[iy] = np.nan
        except Exception:
            A[iy] = np.nanmean(a[s:e])
        if iy == 0:
            rv[:n] = r[s:e]
    return V, R, A, rv

# ----------------------------
# Referencing / shifting energy
# ----------------------------

def compute_ref_shift(Es, r, rows):
    """
    Reference shift: choose, for each row, the energy at its maximum r, then
    return the minimum among those as the asymptotic reference.
    """
    refs = []
    for s, e in rows:
        if e > s:
            irel = np.argmax(r[s:e])
            refs.append(Es[s + irel])
    if len(refs) == 0:
        return 0.0
    return float(np.nanmin(refs))

def compute_shift_from_grid(V):
    """Simpler reference: take nan-min of the last column (max distance) across rows."""
    if V.size == 0:
        return 0.0
    col = V[:, -1]
    if not np.any(np.isfinite(col)):
        # Avoid all-NaN warnings and return zero shift so V stays NaN (diagnostic)
        return 0.0
    return float(np.nanmin(col))

# ----------------------------
# Plotting
# ----------------------------

def plot_imshow(V, rv, A, emin=None, vmax=None, title=None, cmap='bwr', kcal=False, ax=None, bColorbar=True, rtick_step=5):
    fac = 23.060548 if kcal else 1.0
    Z = V * fac
    # Build extent from finite rv/A; avoid identical y-limits
    extent = None
    x_label = 'angle a [deg]'
    y_label = 'r'
    # Invert axes compared to previous: x <- A (angles), y <- rv (distance)
    xr = A[np.isfinite(A)]   if A  is not None else np.array([])
    yr = rv[np.isfinite(rv)] if rv is not None else np.array([])
    if xr.size >= 1 and yr.size >= 2:
        # Use min/max to avoid inverted ranges if rows are not angle-sorted
        x0 = np.nanmin(xr)
        x1 = np.nanmax(xr)
        y0 = np.nanmin(yr)
        y1 = np.nanmax(yr)
        if np.isfinite(x0) and np.isfinite(x1) and np.isfinite(y0) and np.isfinite(y1):
            # pad if needed to avoid singular transform
            if abs(x1 - x0) < 1e-9:
                padx = 0.5 if x0 != 0 else 1.0
                x0, x1 = x0 - padx, x0 + padx
            if abs(y1 - y0) < 1e-9:
                pady = 0.5 if y0 != 0 else 1.0
                y0, y1 = y0 - pady, y0 + pady
            extent = [x0, x1, y0, y1]
    # If everything is NaN, render an empty panel with a note
    if not np.any(np.isfinite(Z)):
        if ax is None:
            fig = plt.figure(figsize=(7, 5)); ax = plt.gca()
        ax.set_axis_off()
        ax.text(0.5, 0.5, 'No finite data', ha='center', va='center', transform=ax.transAxes)
        if title: ax.set_title(title)
        return None
    # Color scale handling
    vmin = None
    if vmax is None and emin is not None:
        # Interpret `emin` as symmetric magnitude
        vmag = abs(emin)
        vmin = -vmag
        vmax = +vmag
    else:
        vmin = emin if emin is not None else np.nanmin(Z)
    if vmax is None:
        vmax = np.nanmax(Z)
    if ax is None:
        fig = plt.figure(figsize=(7, 5))
        ax = plt.gca()
    # Transpose so that x (width) corresponds to angles (A), y (height) to distances (rv)
    im = ax.imshow(Z.T, origin='lower', aspect='auto', extent=extent, vmin=vmin, vmax=vmax, cmap=cmap)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    # Custom y ticks for non-uniform rv sampling: label every Nth sample at row centers
    try:
        if (rtick_step is not None) and (rtick_step > 0) and (rv is not None) and (rv.size > 0):
            M = int(rv.size)
            mask = np.isfinite(rv)
            valid = np.nonzero(mask)[0]
            if valid.size > 0:
                sel = valid[::int(rtick_step)]
                if extent is None:
                    y0, y1 = 0.0, float(M)
                else:
                    y0, y1 = float(extent[2]), float(extent[3])
                dy = (y1 - y0) / float(M)
                yticks = y0 + (sel + 0.5) * dy
                ylabels = [f"{rv[i]:.2f}" for i in sel]
                ax.set_yticks(yticks)
                ax.set_yticklabels(ylabels)
    except Exception:
        pass
    if bColorbar:
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('E [kcal/mol]' if kcal else 'E [eV]')
    if title:
        ax.set_title(title)
    return im

def plot_polar(V, rv, A, emin=None, vmax=None, title=None, cmap='bwr', kcal=False, ax=None, bColorbar=True, rmax=None, half='right', R=None):
    """Polar plot using angles A (deg) and distances rv. V is [ny,nx].
    half: 'right' shows -90..+90 deg, 'left' shows 90..270 deg.
    """
    fac = 23.060548 if kcal else 1.0
    # Sort rows by angle to make contours well-behaved
    order = np.argsort(A)
    A = A[order]
    V = V[order, :]
    if R is not None:
        R = R[order, :]
    Z = V * fac
    # Build theta (radians), apply shift like in reference scripts to show -90..+90 within 90..270 window
    thetas = np.radians(A)
    if half == 'right':
        thetas = thetas - np.pi
    # Build coordinate grids matching V's shape. Prefer per-row R if provided.
    ny, nx = V.shape
    if R is None:
        Rg, Tg = np.meshgrid(rv, thetas)
    else:
        Rg = R
        Tg = np.repeat(thetas[:, None], nx, axis=1)
    Zt = np.ma.masked_invalid(Z)
    # Color scale
    vmin = None
    if vmax is None and emin is not None:
        vmag = abs(emin)
        vmin = -vmag
        vmax = +vmag
    else:
        vmin = emin
    if ax is None:
        ax = plt.subplot(111, projection='polar')
    # Plot
    # For contourf, generate levels within [vmin, vmax] and clamp data to avoid outliers skewing levels
    if (vmin is not None) and (vmax is not None):
        levels = np.linspace(vmin, vmax, 100)
        Zp = np.clip(Zt, vmin, vmax)
    else:
        levels = 100
        Zp = Zt
    # If everything masked (all NaN), skip plotting to avoid matplotlib errors
    if not np.any(np.isfinite(Zp)):
        if ax is None:
            ax = plt.subplot(111, projection='polar')
        ax.text(0.5, 0.5, 'No finite data', transform=ax.transAxes, ha='center', va='center')
        if title: ax.set_title(title)
        return ax
    cs = ax.contourf(Tg, Rg, Zp, levels=levels, cmap=cmap, vmin=vmin, vmax=vmax)
    # Half circle setup: always use 90..270 window; for 'right' we shifted by -pi above
    ax.set_thetamin(90)
    ax.set_thetamax(270)
    # Radial limits: include center to show blank area below smallest sampled radius
    base = R if R is not None else rv
    if base is None or not np.any(np.isfinite(base)):
        rmax_eff = rmax if (rmax is not None) else 1.0
    else:
        rmax_eff = rmax if (rmax is not None) else float(np.nanmax(base))
    ax.set_ylim(0.0, rmax_eff)
    if bColorbar:
        cbar = plt.colorbar(cs, ax=ax)
        cbar.set_label('E [kcal/mol]' if kcal else 'E [eV]')
    if title:
        ax.set_title(title)
    return cs

def plot_profiles(V, rv, A, R=None, rmax=None, kcal=False, ax=None, title=None, vmin=None, vmax=None):
    """Plot multiple 1D profiles on a single axes for given 2D grid V[r, a].
    Plots:
    - Radial profiles at selected angles (nearest to -90 deg and 0 deg)
    - Angular profile at radius of the global energy minimum
    - Per-angle minimum energy (min over r for each angle)

    Overlays angular profiles using theta in radians in [0, 2π) so they can share
    the same x-axis with r if typical rmax ~ 6-7.
    """
    fac = 23.060548 if kcal else 1.0
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(6, 4))
    # Sort by angle for consistent ordering
    order = np.argsort(A)
    A = A[order]
    V = V[order, :]
    if R is not None:
        R = R[order, :]
    # Helper: pick nearest angle index
    def nearest_angle_idx(target_deg):
        return int(np.nanargmin(np.abs(A - target_deg)))
    # Radial profiles at -90 and 0 deg (use nearest available)
    idxs = []
    for tdeg in (-90.0, 0.0):
        try:
            idxs.append(nearest_angle_idx(tdeg))
        except Exception:
            pass
    # Deduplicate
    idxs = sorted(set(idxs))
    colors = ['tab:blue', 'tab:red']
    for j, irow in enumerate(idxs):
        rr = (R[irow, :] if R is not None else rv).astype(float)
        ee = (V[irow, :] * fac).astype(float)
        if (vmin is not None) and (vmax is not None):
            ee = np.clip(ee, vmin, vmax)
        if rmax is not None:
            m = rr <= float(rmax)
            rr = rr[m]; ee = ee[m]
        ax.plot(rr, ee, '-', lw=2, color=colors[j % len(colors)], label=f"radial @{A[irow]:.0f}°")
    # Angular profile at global energy minimum (pick column index of global min)
    try:
        iy, ix = np.unravel_index(np.nanargmin(V), V.shape)
        theta = np.radians(A)
        theta = np.mod(theta, 2*np.pi)  # wrap to [0, 2π)
        srt = np.argsort(theta)
        theta = theta[srt]
        e_ang = (V[:, ix] * fac)[srt]
        if (vmin is not None) and (vmax is not None):
            e_ang = np.clip(e_ang, vmin, vmax)
        ax.plot(theta, e_ang, '--', lw=2, color='k', label=f"angular @ r≈{(rv[ix] if rv is not None else (R[iy, ix] if R is not None else np.nan)):.2f}")
    except Exception:
        pass
    # Per-angle minima (min over r for each angle)
    try:
        e_min_ang = np.nanmin(V, axis=1) * fac
        if (vmin is not None) and (vmax is not None):
            e_min_ang = np.clip(e_min_ang, vmin, vmax)
        theta = np.radians(A)
        theta = np.mod(theta, 2*np.pi)
        srt = np.argsort(theta)
        ax.plot(theta[srt], e_min_ang[srt], '-', lw=1.5, color='tab:green', label='min over r per angle')
    except Exception:
        pass
    # Cosmetics
    if rmax is not None:
        ax.set_xlim(0.0, float(rmax))
    if (vmin is not None) and (vmax is not None):
        ax.set_ylim(float(vmin), float(vmax))
    ax.grid(True, ls=':')
    ax.set_xlabel('r [Å] / θ [rad]')
    ax.set_ylabel('E [kcal/mol]' if kcal else 'E [eV]')
    if title:
        ax.set_title(title)
    ax.legend(loc='best', fontsize=8)
    return ax

def parse_panel_list(list_path):
    """
    Parse a panel list file with format:
      Nrows Ncols\n
      --- (row separator)\n
      filename.xyz or #comment to skip.
    Returns (nrows, ncols, entries) where entries is a flat list of either
    absolute file paths or None for skipped panels.
    Paths are resolved relative to the list file directory.
    """
    p = Path(list_path)
    base = p.parent
    with open(p, 'r') as f:
        lines = [ln.strip() for ln in f if ln.strip() != '']
    # First non-empty line must have two ints
    rc = lines[0].split()
    nrows, ncols = int(rc[0]), int(rc[1])
    entries = []
    r = 0
    c = 0
    for ln in lines[1:]:
        if ln.startswith('---'):
            # start next row
            if c != 0 and c != ncols:
                # pad incomplete row
                while c < ncols:
                    entries.append(None); c += 1
            r += 1
            c = 0
            continue
        if ln.startswith('#'):
            entries.append(None)
        else:
            entries.append(str((base / ln).resolve()))
        c += 1
        if c == ncols:
            c = 0
    # Ensure size nrows*ncols
    if len(entries) < nrows * ncols:
        entries += [None] * (nrows * ncols - len(entries))
    else:
        entries = entries[: nrows * ncols]
    return nrows, ncols, entries

# ----------------------------
# Main
# ----------------------------

def compute_panel_data(xyz, natoms=None, debug=False, unsigned_angle=False):
    # Header-derived values (may be NaN where missing)
    Eh, Rh, Ah = parse_headers_ra(xyz)
    Es, Ps = read_scan_atomicutils(xyz)
    if Es.size == 0:
        print(f"WARNING: atomicUtils.scan_xyz() failed to parse {xyz} => fallback to local parse_xyz_blocks()")
        _, _, Ps_local = parse_xyz_blocks(xyz, natoms=natoms)
        Es = np.full((Ps_local.shape[0],), np.nan)
        Ps = Ps_local
    if debug:
        print(f"Parsed {len(Es)} blocks, natoms={Ps.shape[1]} from {Path(xyz).name}")
    # Prefer energies from header if available
    if Eh.size == Es.size and np.any(np.isfinite(Eh)):
        Es = Eh
    # Geometry r/a, then override from headers where provided
    r, a = compute_ra_vec(Ps, signed=(not unsigned_angle))
    if Rh.size == r.size and np.any(np.isfinite(Rh)):
        r = np.where(np.isfinite(Rh), Rh, r)
    if Ah.size == a.size and np.any(np.isfinite(Ah)):
        a_h = np.abs(Ah) if unsigned_angle else Ah
        a = np.where(np.isfinite(a_h), a_h, a)
    rows, step = detect_rows_by_r(r)
    if debug:
        print(f"Detected {len(rows)} rows; typical dr={step:.6f}")
        try:
            #print(f"  r[min,max] {np.nanmin(r):.3f}, {np.nanmax(r):.3f} r[:]: ", np.array2string(r, precision=3) )
            #print(f"  a[min,max] {np.nanmin(a):.1f}, {np.nanmax(a):.1f} a[:]: ", np.array2string(a, precision=1))
            print(f"  r[min,max] {np.nanmin(r):.3f}, {np.nanmax(r):.3f}")
            print(f"  a[min,max] {np.nanmin(a):.1f}, {np.nanmax(a):.1f}")
            print(f"  Es finite={np.isfinite(Es).sum()} / {Es.size}")
            print(f"  Rh finite={np.isfinite(Rh).sum()} Ah finite={np.isfinite(Ah).sum()}")
        except Exception:
            pass
    Vraw, R, A, rv = reshape_to_grid(Es, r, a, rows)
    shift = compute_shift_from_grid(Vraw)
    V = Vraw - shift
    if debug:
        try:
            ny, nx = V.shape
            rv_f = rv[np.isfinite(rv)] if rv is not None else np.array([])
            A_f  = A[np.isfinite(A)]   if A  is not None else np.array([])
            print(f"  Grid shape: ny={ny}, nx={nx}")
            if rv_f.size:
                print(f"  rv[len={rv_f.size}] min,max = {np.nanmin(rv_f):.3f}, {np.nanmax(rv_f):.3f}   rv[:] = ", np.array2string(rv_f, precision=3))
            if A_f.size:
                print(f"  A[len={A_f.size}] min,max  = {np.nanmin(A_f):.1f}, {np.nanmax(A_f):.1f}   A[:]  = ", np.array2string(A_f, precision=1))
            print(f"  V finite={np.isfinite(V).sum()} / {V.size}   last_col finite={np.isfinite(V[:, -1]).sum()} / {V.shape[0]}")
        except Exception:
            pass
    return V, rv, A, shift, R

def plot_list(list_path, emin=None, emax=None, sym=False, kcal=False, cmap='bwr', bColorbar=True, natoms=None, debug=False, unsigned_angle=False, transpose=False, polar=False, rmax=None, half='right', lines=False, rtick_step=5):
    nrows, ncols, entries = parse_panel_list(list_path)
    # For line mode we want normal Cartesian axes; ignore polar projection
    subplot_kw = {'projection':'polar'} if (polar and not lines) else None
    if transpose:
        # swap rows/cols for layout
        fig, axs = plt.subplots(ncols, nrows, figsize=(3.0*nrows, 2.6*ncols), squeeze=False, subplot_kw=subplot_kw)
    else:
        fig, axs = plt.subplots(nrows, ncols, figsize=(3.0*ncols, 2.6*nrows), squeeze=False, subplot_kw=subplot_kw)
    # Draw each panel independently, with its own colorbar and autoscaled range
    for idx, fp in enumerate(entries):
        r = idx // ncols
        c = idx % ncols
        ax = axs[c, r] if transpose else axs[r, c]
        if fp is None:
            ax.set_axis_off()
            continue
        try:
            V, rv, A, shift, R = compute_panel_data(fp, natoms=natoms, debug=debug, unsigned_angle=unsigned_angle)
        except Exception as e:
            print(f"ERROR processing {fp}: {e}")
            ax.set_axis_off(); continue
        title = Path(fp).name
        # Determine vmin/vmax per image according to options precedence: sym > emax > emin/autoscale
        fac = 23.060548 if kcal else 1.0
        if sym:
            # Guard: if all-NaN, skip this panel
            if not np.any(np.isfinite(V)):
                print(f"WARNING: all-NaN V for {title}; skipping panel")
                ax.set_axis_off(); continue
            emin_img = float(np.nanmin(V)) * fac
            vmag = abs(emin_img)
            vmin_plot = -vmag
            vmax_plot = +vmag
        elif emax is not None and emax > 0:
            if not np.any(np.isfinite(V)):
                print(f"WARNING: all-NaN V for {title}; skipping panel")
                ax.set_axis_off(); continue
            vmin_plot = float(np.nanmin(V)) * fac
            vmax_plot = vmin_plot + emax
        else:
            vmin_plot = emin  # may be None
            vmax_plot = None
        if lines:
            plot_profiles(V, rv, A, R=R, rmax=rmax, kcal=kcal, ax=ax, title=title, vmin=vmin_plot, vmax=vmax_plot)
        elif polar:
            plot_polar(V, rv, A, emin=vmin_plot, vmax=vmax_plot, title=title, cmap=cmap, kcal=kcal, ax=ax, bColorbar=bColorbar, rmax=rmax, half=half, R=R)
        else:
            plot_imshow(V, rv, A, emin=vmin_plot, vmax=vmax_plot, title=title, cmap=cmap, kcal=kcal, ax=ax, bColorbar=bColorbar, rtick_step=rtick_step)
    fig.tight_layout()
    return fig

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