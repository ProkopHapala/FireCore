"""
Pure-Python plotting utilities for FitREQ energy maps.
No C/OpenCL dependency — only numpy + matplotlib + re + os.

Copied verbatim from pyBall/FitREQ.py (pure-Python functions only).
"""
from pathlib import Path
import numpy as np
import os
import re
import matplotlib.pyplot as plt

from pyBall.atomicUtils import scan_xyz

ev2kcal = 23.060547831


# ─────────────────────────────────────────────────────────
# Data I/O and slicing (pure Python, from FitREQ.py)
# ─────────────────────────────────────────────────────────

def read_xyz_data(fname="input_all.xyz"):
    """Extract Erefs and x0s from a concatenated XYZ file comment lines."""
    Erefs = []
    x0s = []
    try:
        with open(fname,'r') as f:
            for line in f:
                ws = line.split()
                if len(ws)>=10 and ws[0]=="#" and ws[1]=="n0":
                    Erefs.append(float(ws[4]))
                    x0s.append(float(ws[6]))
    except Exception as e:
        print(f"Error reading {fname}: {e}")
        return np.array([]), np.array([])
    return np.array(Erefs), np.array(x0s)


def mark_molecule_blocks(lines):
    """
    Parse comment lines from an all.xyz / concatenated XYZ file
    to determine molecule block boundaries and how many points per angle.
    """
    i  = 0
    marks    = []
    angle_data = []
    blocks={}
    for l in lines:
        ws = l.split()
        if len(ws)<10: continue
        name= ws[-1] if ws[-1]!='eV' else ws[-2]
        i+=1
        if not name in blocks: blocks.setdefault(name,[]).append(i)
        else:                  blocks[    name].append(i)
    for name,is_ in blocks.items():
        marks.append((is_[0]-1,is_[-1]))
        # count consecutive angle groups
        preva=None
        nxs=[]
        nx=0
        for ii in is_:
            ws = lines[ii-1].split()
            a = int(float(ws[8]))
            if preva is None: preva=a
            if a==preva: nx+=1
            else:
                nxs.append(nx)
                nx=1
                preva=a
        nxs.append(nx)
        angle_data.append(nxs)
    # correct order by first index
    marks2=sorted(marks)
    angd2=[angle_data[marks.index(m)] for m in marks2]
    return marks2, angd2


def slice_and_reshape(Es, marks, angle_data):
    """
    Reshape flat energy array E into panels, one per molecule block.
    Each panel: (n_angles, n_distances) with NaN for missing points.
    """
    Eplots = []
    for (i0,i1),nxs in zip(marks,angle_data):
        nx = np.max(nxs)
        ny = len(nxs)
        Eplot = np.empty((ny,nx)); Eplot[:]=np.nan
        ii = i0
        for iy,nxi in enumerate(nxs):
            Eplot[iy,0:nxi] = Es[ii:ii+nxi]
            ii+=nxi
        Eplots.append(Eplot)
    return Eplots


def concatenate_xyz_files(directories=None, base_path='./', fname="all.xyz", output_file="all.xyz", mode='w'):
    """
    Concatenate per-molecule all.xyz files into one file,
    tagging comment lines with source directory name.
    Returns marks: list of (i0,i1) per molecule block.
    """
    marks=[]
    if directories is None:
        directories = find_all_dirs(base_path)
    if isinstance(directories,str):
        directories=[directories]
    if mode == 'w':
        fout=open(output_file, 'w')
    else:
        fout=open(output_file, 'a')
    i=0
    for dir_ in directories:
        #print( "dir_=", dir_, fname)
        path1 = base_path + "/" + dir_ + "/" + fname
        path2 = dir_       + "/" + fname
        path = path1 if os.path.exists(path1) else path2
        i0=i
        with open(path, 'r') as f:
            for line in f:
                if line[0] == '#':
                    line =  line.replace("\n","") + " " + dir_ + "\n"
                i+=1
                fout.write(line)
        i1=i
        marks.append((i0,i1))
    fout.close()
    return marks


def concatenate_xyz_files_flat(names=None, base_path='./', output_file="all.xyz", mode='w'):
    """
    Concatenate .xyz files from a flat directory (no subdirs).
    """
    marks=[]
    if mode == 'w':
        fout=open(output_file, 'w')
    else:
        fout=open(output_file, 'a')
    i=0
    for name in names:
        path = name if os.path.exists(name) else base_path+"/"+name
        i0=i
        with open(path, 'r') as f:
            for line in f:
                if line[0] == '#':
                    line = line.replace("\n","") + " " + name + "\n"
                i+=1
                fout.write(line)
        i1=i
        marks.append((i0,i1))
    fout.close()
    return marks


def find_all_dirs(base_path):
    """Return all subdirectories containing an all.xyz file."""
    dirs = [d for d in os.listdir(base_path)
            if os.path.isdir(os.path.join(base_path, d))]
    dirs_with_xyz = []
    for d in dirs:
        if os.path.exists(os.path.join(base_path, d, "all.xyz")):
            dirs_with_xyz.append(d)
    return dirs_with_xyz


def read_file_comments(fname, comment_sign='#'):
    """Read only comment lines from an XYZ file."""
    comments=[]
    with open(fname,'r') as f:
        for line in f:
            if line[0]==comment_sign:
                comments.append(line.replace("\n",""))
    return comments


def extract_comments_and_types(fname, comment_sign='#'):
    """Extract atom type names and comment lines from XYZ."""
    type_names=[]
    comments=[]
    with open(fname,'r') as f:
        il=-1
        for line in f:
            if il==-1:
                try:
                    n=int(line.split()[0])
                    il=n
                except:
                    continue
            elif il==n:
                ws=line.split()
                if len(ws) >=5:
                    comments.append( ws[4] + " " + ws[5] + " " + ws[6] )
                il-=1
            elif il >0:
                ws=line.split()
                if ws[0] not in type_names:
                    type_names.append(ws[0])
                il-=1
        return type_names, comments


# ─────────────────────────────────────────────────────────
# Grid operations (pure numpy)
# ─────────────────────────────────────────────────────────

def shift_grid(G):
    """Shift grid by subtracting asymptotic baseline (min at largest distance)."""
    # Clean garbage
    G[np.abs(G) > 1e10] = np.nan
    last = G[-1, :]
    finite = np.isfinite(last)
    if finite.any():
        ref = np.nanmin(last[finite])
    else:
        ref = 0.0
    GS = G - ref
    mloc = np.nanmin(GS)
    return GS, ref, mloc


def extract_min_curves(angles, distances, G, rmax=None):
    """
    For each angle column, find distance and energy at minimum.
    G: (n_distances, n_angles) grid.
    """
    nA = len(angles)
    rmin = np.full(nA, np.nan)
    emin = np.full(nA, np.nan)
    for j in range(nA):
        col = G[:, j]
        finite = np.isfinite(col)
        if finite.any():
            idx = np.nanargmin(col)
            rmin[j] = distances[idx]
            emin[j] = col[idx]
            if rmax is not None and rmin[j] > rmax:
                rmin[j] = np.nan
                emin[j] = np.nan
    return rmin, emin


def compute_min_lines_from_panel(Epanel, Xpanel, angles, rmax=None, do_shift=True):
    """Extract rmin(angle) and emin(angle) from panel-shaped data."""
    G = Epanel.T
    ny, nx = G.shape
    distances = np.full(nx, np.nan)
    for j in range(nx):
        col = Xpanel[:, j] if Xpanel.ndim == 2 else Xpanel
        finite = np.isfinite(col)
        if finite.any():
            distances[j] = col[finite][0]
    if do_shift:
        GS, _, _ = shift_grid(G)
    else:
        GS = G
    return extract_min_curves(angles, distances, GS, rmax)


def _distances_from_Xpanel(Xpanel):
    """Extract 1D distance array from a panel."""
    ny, nx = Xpanel.shape
    distances = np.full(nx, np.nan)
    for j in range(nx):
        col = Xpanel[:, j]
        finite = np.isfinite(col)
        if finite.any():
            distances[j] = col[finite][0]
    return distances


# ─────────────────────────────────────────────────────────
# Data export helpers
# ─────────────────────────────────────────────────────────

def save_grid_npz(angles, distances, grid, filepath):
    np.savez(filepath, angles=angles, distances=distances, grid=grid)

def save_grid_gnuplot(angles, distances, grid, filepath):
    with open(filepath, 'w') as f:
        f.write("# angle distance value\n")
        ny, nx = grid.shape
        for i in range(ny):
            for j in range(nx):
                v = grid[i, j]
                if not np.isnan(v):
                    f.write(f"{angles[j]:8.1f} {distances[i]:8.3f} {v:15.8e}\n")

def save_min_lines_npz(angles, rmin, emin, filepath):
    np.savez(filepath, angles=angles, rmin=rmin, emin=emin)

def save_min_lines_gnuplot(angles, rmin, emin, filepath):
    with open(filepath, 'w') as f:
        f.write("# angle rmin emin\n")
        for i in range(len(angles)):
            if not np.isnan(rmin[i]):
                f.write(f"{angles[i]:8.1f} {rmin[i]:8.3f} {emin[i]:15.8e}\n")

# ─────────────────────────────────────────────────────────
# Plotting functions
# ─────────────────────────────────────────────────────────

def plot_Epanels(Eplots, ref_dirs, bColorbar=True, Emin=-5.0, bKcal=False):
    E_units = 1.0
    if bKcal:
        E_units = ev2kcal
    nmols = len(Eplots)
    if nmols != len(ref_dirs):
        print(f"Error: len(Eplots={nmols}) != len(ref_dirs={len(ref_dirs)})")
        return
    fig, axs = plt.subplots(1, nmols, figsize=(20, 3))
    for i in range(nmols):
        im = axs[i].imshow(Eplots[i].T * E_units, aspect='auto', origin='lower',
                           vmin=Emin, vmax=-Emin, cmap='bwr')
        axs[i].set_title(f"Ref: {ref_dirs[i]}")
        axs[i].set_ylabel('Reference Energies')
        if bColorbar:
            plt.colorbar(im, ax=axs[i])
    plt.tight_layout()
    return fig


def plot_Epanels_diff(Emodels, Erefs, ref_dirs, bColorbar=True, Emin=-5.0, bKcal=False):
    E_units = 1.0
    if bKcal:
        E_units = ev2kcal
    nmols = len(Erefs)
    if nmols == 1:
        fig, axs = plt.subplots(3, 1, figsize=(6, 9))
        for i in range(nmols):
            im0 = axs[0].imshow(Erefs[i].T * E_units, aspect='auto', origin='lower',
                                vmin=Emin, vmax=-Emin, cmap='bwr')
            axs[0].set_title(f"Ref: {ref_dirs[i]}")
            if bColorbar: plt.colorbar(im0, ax=axs[0])

            im1 = axs[1].imshow(Emodels[i].T * E_units, aspect='auto', origin='lower',
                                vmin=Emin, vmax=-Emin, cmap='bwr')
            axs[1].set_title(f"Model: {ref_dirs[i]}")
            if bColorbar: plt.colorbar(im1, ax=axs[1])

            Ediff = (Emodels[i] - Erefs[i]) * E_units
            dmax = np.nanmax(np.abs(Ediff))
            im2 = axs[2].imshow(Ediff.T, aspect='auto', origin='lower',
                                vmin=-dmax, vmax=dmax, cmap='bwr')
            axs[2].set_title(f"Diff: {ref_dirs[i]}")
            if bColorbar: plt.colorbar(im2, ax=axs[2])
    else:
        fig, axs = plt.subplots(3, nmols, figsize=(6 * nmols, 9))
        for i in range(nmols):
            im0 = axs[0, i].imshow(Erefs[i].T * E_units, aspect='auto', origin='lower',
                                   vmin=Emin, vmax=-Emin, cmap='bwr')
            axs[0, i].set_title(f"Ref: {ref_dirs[i]}")
            if bColorbar: plt.colorbar(im0, ax=axs[0, i])

            im1 = axs[1, i].imshow(Emodels[i].T * E_units, aspect='auto', origin='lower',
                                   vmin=Emin, vmax=-Emin, cmap='bwr')
            axs[1, i].set_title(f"Model: {ref_dirs[i]}")
            if bColorbar: plt.colorbar(im1, ax=axs[1, i])

            Ediff = (Emodels[i] - Erefs[i]) * E_units
            dmax = np.nanmax(np.abs(Ediff))
            im2 = axs[2, i].imshow(Ediff.T, aspect='auto', origin='lower',
                                   vmin=-dmax, vmax=dmax, cmap='bwr')
            axs[2, i].set_title(f"Diff: {ref_dirs[i]}")
            if bColorbar: plt.colorbar(im2, ax=axs[2, i])
    plt.tight_layout()
    return fig


def plot2Dlong_diff(Erefs, Es, lens):
    nmols = len(lens)
    fig, axs = plt.subplots(2, nmols, figsize=(6 * nmols, 8))
    if nmols == 1:
        axs = axs.reshape(2, 1)
    for i in range(nmols):
        im0 = axs[0, i].imshow(Erefs[i].T, aspect='auto', origin='lower', cmap='bwr')
        axs[0, i].set_title(f"Ref {i}")
        plt.colorbar(im0, ax=axs[0, i])
        im1 = axs[1, i].imshow((Es[i] - Erefs[i]).T, aspect='auto', origin='lower', cmap='bwr')
        axs[1, i].set_title(f"Diff {i}")
        plt.colorbar(im1, ax=axs[1, i])
    plt.tight_layout()
    return fig


def plot_Epanels_diff_separate(Emodels, Erefs, ref_dirs, save_prefix=None,
                               bColorbar=True, Emin=-5.0, bKcal=False, bClose=False):
    E_units = 1.0
    if bKcal:
        E_units = ev2kcal
    nmols = len(Erefs)
    figs = []
    for i in range(nmols):
        fig, axs = plt.subplots(1, 3, figsize=(18, 4))
        im0 = axs[0].imshow(Erefs[i].T * E_units, aspect='auto', origin='lower',
                            vmin=Emin, vmax=-Emin, cmap='bwr')
        axs[0].set_title(f"Ref: {ref_dirs[i]}")
        if bColorbar: plt.colorbar(im0, ax=axs[0])

        im1 = axs[1].imshow(Emodels[i].T * E_units, aspect='auto', origin='lower',
                            vmin=Emin, vmax=-Emin, cmap='bwr')
        axs[1].set_title(f"Model: {ref_dirs[i]}")
        if bColorbar: plt.colorbar(im1, ax=axs[1])

        Ediff = (Emodels[i] - Erefs[i]) * E_units
        dmax = np.nanmax(np.abs(Ediff))
        im2 = axs[2].imshow(Ediff.T, aspect='auto', origin='lower',
                            vmin=-dmax, vmax=dmax, cmap='bwr')
        axs[2].set_title(f"Diff: {ref_dirs[i]}")
        if bColorbar: plt.colorbar(im2, ax=axs[2])
        plt.tight_layout()
        if save_prefix:
            plt.savefig(f"{save_prefix}_{i}.png")
        if bClose:
            plt.close(fig)
        figs.append(fig)
    return figs


def plot_energy_2d_from_xyz(xyz_path, distances=None, angles=None, title=None,
                            cmap='bwr', vmin=None, vmax=None, save_path=None):
    """Read XYZ, build 2D grid, shift by baseline, and plot."""
    if distances is None:
        distances = np.arange(1.40, 20.05, 0.05)
    if angles is None:
        angles = list(range(-90, 100, 10))

    dist_to_ix = {d: idx for idx, d in enumerate(distances)}
    ang_to_iy = {a: idx for idx, a in enumerate(angles)}

    G = np.full((len(distances), len(angles)), np.nan)
    seq = []
    direction_found = None

    with open(xyz_path, 'r') as f:
        while True:
            line = f.readline()
            if not line: break
            line = line.strip()
            if not line: continue
            try:
                natoms = int(line.split()[0])
            except ValueError:
                continue
            comment = f.readline()
            if not comment: break
            for _ in range(natoms):
                f.readline()

            etot, x0, angle, direction = _parse_comment(comment)
            if direction is not None and direction_found is None:
                direction_found = direction
            if etot is not None and x0 is not None and angle is not None:
                i_dist = dist_to_ix.get(x0)
                i_ang = ang_to_iy.get(angle)
                if i_dist is not None and i_ang is not None:
                    G[i_dist, i_ang] = etot
                    seq.append((i_dist, i_ang))

    # Baseline shift
    last = G[-1, :]
    finite = np.isfinite(last)
    if finite.any():
        ref_val = np.nanmin(last[finite])
    else:
        ref_val = 0.0
    GS = G - ref_val
    mloc = np.nanmin(GS)

    if vmin is None:
        vmin = np.nanmin(GS)
        if vmin > 0: vmin = 0
    if vmax is None:
        vmax = -vmin

    fig, ax = plt.subplots(figsize=(8, 5))
    extent = [angles[0], angles[-1], distances[0], distances[-1]]
    im = ax.imshow(GS, origin='lower', aspect='auto', cmap=cmap,
                   vmin=vmin, vmax=vmax, interpolation='nearest',
                   extent=extent)
    plt.colorbar(im, ax=ax, label="Energy [a.u.]")
    ax.set_xlabel(f"angle ({direction_found or '?'}) [deg]")
    ax.set_ylabel("distance x0 [A]")
    if title:
        ax.set_title(title)
    else:
        ax.set_title(os.path.basename(xyz_path))
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=200)
    return GS


def _parse_comment(comment):
    """Parse comment line for Etot, x0, angle, and axis direction."""
    etot = None; x0 = None; angle = None; direction = None
    toks = comment.strip().split()
    for i, t in enumerate(toks):
        if t == 'Etot' and i + 1 < len(toks):
            try: etot = float(toks[i + 1])
            except ValueError: pass
        if t == 'x0' and i + 1 < len(toks):
            try: x0 = float(toks[i + 1])
            except ValueError: pass
        if t in ('y', 'Y', 'z', 'Z') and i + 1 < len(toks):
            try:
                angle = int(float(toks[i + 1]))
                direction = t.lower()
            except ValueError: pass
    return etot, x0, angle, direction


def parse_xyz_mapping(xyz_path, distances=None, angles=None):
    """Parse XYZ file to get reference grid and frame-to-grid mapping."""
    if distances is None:
        distances = np.arange(1.40, 20.05, 0.05)
    if angles is None:
        angles = list(range(-90, 100, 10))

    dist_to_ix = {d: idx for idx, d in enumerate(distances)}
    ang_to_iy = {a: idx for idx, a in enumerate(angles)}

    G = np.full((len(distances), len(angles)), np.nan)
    seq = []
    axis = None

    with open(xyz_path, 'r') as f:
        while True:
            line = f.readline()
            if not line: break
            line = line.strip()
            if not line: continue
            try:
                natoms = int(line.split()[0])
            except ValueError:
                continue
            comment = f.readline()
            if not comment: break
            for _ in range(natoms):
                f.readline()

            etot, x0, angle, direction = _parse_comment(comment)
            if direction is not None and axis is None:
                axis = direction
            if etot is not None and x0 is not None and angle is not None:
                x0r = round(x0 / 0.05) * 0.05
                ar = round(angle / 10) * 10
                i_dist = dist_to_ix.get(x0r)
                i_ang = ang_to_iy.get(ar)
                if i_dist is not None and i_ang is not None:
                    G[i_dist, i_ang] = etot
                    seq.append((i_dist, i_ang))
    return G, seq, axis, distances, angles



def shift_grid(G):
    """Shift grid by subtracting asymptotic baseline (min at largest distance)."""
    G[np.abs(G) > 1e10] = np.nan
    last = G[-1, :]
    finite = np.isfinite(last)
    if finite.any():
        ref = np.nanmin(last[finite])
    else:
        ref = 0.0
    GS = G - ref
    mloc = np.nanmin(GS)
    return GS, ref, mloc


def extract_min_curves(angles, distances, G, rmax=None):
    """For each angle column, find distance and energy at minimum."""
    nA = len(angles)
    rmin = np.full(nA, np.nan)
    emin = np.full(nA, np.nan)
    for j in range(nA):
        col = G[:, j]
        finite = np.isfinite(col)
        if finite.any():
            idx = np.nanargmin(col)
            rmin[j] = distances[idx]
            emin[j] = col[idx]
            if rmax is not None and rmin[j] > rmax:
                rmin[j] = np.nan
                emin[j] = np.nan
    return rmin, emin


def plot_compare(Gref, Gmodel, angles, distances, title, save_prefix=None,
                vmin=None, vmax=None, line=False, kcal=False,
                save_data_prefix=None, save_fmt="both"):
    """Plot reference/model/difference 2D maps, optionally with min lines."""
    conv = ev2kcal if kcal else 1.0
    unit = "kcal" if kcal else "eV"

    GRS, refR, mlocR = shift_grid(Gref)
    GMS, refM, mlocM = shift_grid(Gmodel)

    GRS *= conv
    GMS *= conv

    if vmin is None:
        vmin = np.nanmin([GRS, GMS])
        if vmin > 0: vmin = 0
    if vmax is None:
        vmax = -vmin

    if line:
        EpanelR = GRS.T
        EpanelM = GMS.T
        ny, nx = EpanelR.shape
        Xpanel = np.tile(distances, (ny, 1))  # same distances for all angles
        plot_min_lines_pair(EpanelR, EpanelM, Xpanel, angles, title=title,
                           to_kcal=False, save_fmt=save_fmt)

    fig, axs = plt.subplots(3, 1, figsize=(8, 10))
    extent = [angles[0], angles[-1], distances[0], distances[-1]]

    im0 = axs[0].imshow(GRS, origin='lower', aspect='auto', cmap='bwr',
                        vmin=vmin, vmax=vmax, interpolation='nearest', extent=extent)
    axs[0].set_title(f"Reference ({unit})")
    plt.colorbar(im0, ax=axs[0])

    im1 = axs[1].imshow(GMS, origin='lower', aspect='auto', cmap='bwr',
                        vmin=vmin, vmax=vmax, interpolation='nearest', extent=extent)
    axs[1].set_title(f"Model ({unit})")
    plt.colorbar(im1, ax=axs[1])

    D = GMS - GRS
    dmax = max(-np.nanmin(D), np.nanmax(D))
    im2 = axs[2].imshow(D, origin='lower', aspect='auto', cmap='bwr',
                        vmin=-dmax, vmax=dmax, interpolation='nearest', extent=extent)
    axs[2].set_title(f"Difference ({unit})")
    plt.colorbar(im2, ax=axs[2])

    for ax in axs:
        ax.set_xlabel("Angle (deg)")
        ax.set_ylabel("Distance (A)")

    plt.suptitle(title)
    plt.tight_layout()

    if save_prefix:
        plt.savefig(f"{save_prefix}.png", dpi=200, bbox_inches='tight')
        if save_data_prefix:
            save_data(Gref, Gmodel, angles, distances, save_data_prefix, save_fmt, kcal)

    return fig


def plot_compare_combined(Gref, Gmodel, angles, distances, title, save_path=None,
                          kcal=False, params_text=None):
    """3x2 layout: 2D maps + min lines + text panel."""
    conv = ev2kcal if kcal else 1.0
    unit = "kcal" if kcal else "eV"

    # Shift grids
    GRS, refR, mlocR = shift_grid(Gref)
    GMS, _, _ = shift_grid(Gmodel)

    GRS *= conv
    GMS *= conv

    # Difference
    D = GMS - GRS
    E_MARGIN = 1.0 / conv if not kcal else 1.0
    E_SPAN = 2.0 / conv if not kcal else 2.0
    dmax = max(-np.nanmin(D), np.nanmax(D))

    # Min lines
    EpanelR = GRS.T
    EpanelM = GMS.T
    ny, nx = EpanelR.shape
    Xpanel = np.tile(distances, (ny, 1))
    rR, eR = compute_min_lines_from_panel(EpanelR, Xpanel, angles, do_shift=False)
    rM, eM = compute_min_lines_from_panel(EpanelM, Xpanel, angles, do_shift=False)

    fig = plt.figure(figsize=(14, 12))
    gs = fig.add_gridspec(3, 2, width_ratios=[1, 1])

    extent = [angles[0], angles[-1], distances[0], distances[-1]]

    vmin = min(np.nanmin(GRS), np.nanmin(GMS))
    if vmin > 0: vmin = 0
    vmax = -vmin

    # Ref 2D
    ax00 = fig.add_subplot(gs[0, 0])
    im00 = ax00.imshow(GRS, origin='lower', aspect='auto', cmap='bwr',
                       vmin=vmin, vmax=vmax, interpolation='nearest', extent=extent)
    ax00.set_title(f"Reference ({unit})")
    ax00.set_ylabel("Distance (A)")
    plt.colorbar(im00, ax=ax00)

    # r_min line
    ax01 = fig.add_subplot(gs[0, 1])
    ax01.plot(angles, rR, 'k-', lw=1.5, label="Ref")
    ax01.plot(angles, rM, 'r-', lw=1.5, label="Model")
    ax01.set_title("R$_{min}$(angle)")
    ax01.set_xlabel("Angle (deg)")
    ax01.set_ylabel("R$_{min}$ (A)")
    ax01.legend()
    ax01.grid(alpha=0.2)

    # Model 2D
    ax10 = fig.add_subplot(gs[1, 0])
    im10 = ax10.imshow(GMS, origin='lower', aspect='auto', cmap='bwr',
                       vmin=vmin, vmax=vmax, interpolation='nearest', extent=extent)
    ax10.set_title(f"Model ({unit})")
    ax10.set_ylabel("Distance (A)")
    plt.colorbar(im10, ax=ax10)

    # E_min line
    ax11 = fig.add_subplot(gs[1, 1])
    ax11.plot(angles, eR, 'k-', lw=1.5, label="Ref")
    ax11.plot(angles, eM, 'r-', lw=1.5, label="Model")
    ax11.set_title(f"E$_{{min}}$(angle) ({unit})")
    ax11.set_xlabel("Angle (deg)")
    ax11.set_ylabel(f"E$_{{min}}$ ({unit})")
    ax11.legend()
    ax11.grid(alpha=0.2)

    # Diff 2D
    ax20 = fig.add_subplot(gs[2, 0])
    im20 = ax20.imshow(D, origin='lower', aspect='auto', cmap='bwr',
                       vmin=-dmax, vmax=dmax, interpolation='nearest', extent=extent)
    ax20.set_title(f"Difference ({unit})")
    ax20.set_xlabel("Angle (deg)")
    ax20.set_ylabel("Distance (A)")
    plt.colorbar(im20, ax=ax20)

    # Text panel
    ax21 = fig.add_subplot(gs[2, 1])
    ax21.axis('off')
    if params_text:
        ax21.text(0.05, 0.95, params_text, transform=ax21.transAxes,
                 fontsize=10, verticalalignment='top', fontfamily='monospace')

    plt.suptitle(title)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    return fig


def compute_min_lines_from_panel(Epanel, Xpanel, angles, rmax=None, do_shift=True):
    """Extract rmin(angle) and emin(angle) from panel-shaped data."""
    G = Epanel.T
    ny, nx = G.shape
    distances = np.full(nx, np.nan)
    for j in range(nx):
        col = Xpanel[:, j] if Xpanel.ndim == 2 else Xpanel
        finite = np.isfinite(col)
        if finite.any():
            distances[j] = col[finite][0]
    if do_shift:
        GS, _, _ = shift_grid(G)
    else:
        GS = G
    return extract_min_curves(angles, distances, GS, rmax)


def plot_min_lines_pair(Epanel_ref, Epanel_mod, Xpanel, angles, title=None,
                        save_path=None, to_kcal=False, ms=2, lw=0.5,
                        save_data_prefix=None, save_fmt="both"):
    """Plot rmin(angle) and emin(angle) comparing reference vs model."""
    conv = ev2kcal if to_kcal else 1.0
    rR, eR = compute_min_lines_from_panel(Epanel_ref, Xpanel, angles)
    rM, eM = compute_min_lines_from_panel(Epanel_mod, Xpanel, angles)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
    ax1.plot(angles, rR, 'k.-', ms=ms, lw=lw, label="Ref")
    ax1.plot(angles, rM, 'r.-', ms=ms, lw=lw, label="Model")
    ax1.set_ylim(1.5, 3.0)
    ax1.set_xlabel("Angle (deg)")
    ax1.set_ylabel("R$_{min}$ (A)")
    ax1.legend()
    ax1.grid(alpha=0.2)

    ax2.plot(angles, eR * conv, 'k.-', ms=ms, lw=lw, label="Ref")
    ax2.plot(angles, eM * conv, 'r.-', ms=ms, lw=lw, label="Model")
    ax2.set_xlabel("Angle (deg)")
    ax2.set_ylabel("E$_{min}$")
    ax2.legend()
    ax2.grid(alpha=0.2)

    if title:
        fig.suptitle(title)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=200)
    if save_data_prefix:
        save_min_lines_npz(angles, rR, eR, f"{save_data_prefix}_ref.npz")
        save_min_lines_npz(angles, rM, eM, f"{save_data_prefix}_model.npz")
        if save_fmt in ("gnuplot", "both"):
            save_min_lines_gnuplot(angles, rR, eR, f"{save_data_prefix}_ref.dat")
            save_min_lines_gnuplot(angles, rM, eM, f"{save_data_prefix}_model.dat")
    return fig


def save_data(Gref, Gmodel, angles, distances, save_data_prefix=None, save_fmt="both", kcal=False, line=False):
    """Save grid data to NPZ and/or gnuplot format."""
    if save_data_prefix is None:
        return
    conv = ev2kcal if kcal else 1.0
    GRS = Gref * conv
    GMS = Gmodel * conv
    D = GMS - GRS

    save_grid_npz(angles, distances, GRS, f"{save_data_prefix}_ref.npz")
    save_grid_npz(angles, distances, GMS, f"{save_data_prefix}_model.npz")
    save_grid_npz(angles, distances, D, f"{save_data_prefix}_diff.npz")

    if save_fmt in ("gnuplot", "both"):
        save_grid_gnuplot(angles, distances, GRS, f"{save_data_prefix}_ref.dat")
        save_grid_gnuplot(angles, distances, GMS, f"{save_data_prefix}_model.dat")
        save_grid_gnuplot(angles, distances, D, f"{save_data_prefix}_diff.dat")

    if line:
        for prefix, G in [("_ref", GRS), ("_model", GMS)]:
            Epanel = G.T
            ny, nx = Epanel.shape
            Xpanel = np.tile(distances, (ny, 1))
            rvals, evals = compute_min_lines_from_panel(Epanel, Xpanel, angles, do_shift=False)
            fn = f"{save_data_prefix}{prefix}"
            save_min_lines_npz(angles, rvals, evals, f"{fn}.npz")
            if save_fmt in ("gnuplot", "both"):
                save_min_lines_gnuplot(angles, rvals, evals, f"{fn}.dat")



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
    energy_re = re.compile(r"(?:#\s*E\s*=\s*|\bE_?tot\s+)([+-]?(?:\d*\.\d+|\d+)(?:[Ee][+-]?\d+)?)")
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

def parse_xyz_with_headers(fname, natoms=None):
    """
    Like parse_xyz_blocks, but also returns header-derived info like n0.
    """
    energy_re = re.compile(r"(?:#\s*E\s*=\s*|\bE_?tot\s+)([+-]?(?:\d*\.\d+|\d+)(?:[Ee][+-]?\d+)?)")
    n0_re = re.compile(r"n0\s+(\d+)")

    Es, Ts, Ps, N0s = [], [], [], []
    with open(fname, 'r') as f: lines = f.readlines()
    i = 0
    nline = len(lines)

    # Try infer natoms if not provided
    if natoms is None:
        while i < nline and not lines[i].lstrip().startswith('#'): i += 1
        n_guess = None
        if i-1 >= 0 and lines[i-1].strip().isdigit(): n_guess = int(lines[i-1].strip())
        natoms = n_guess if n_guess is not None else 4
        i = 0 # reset

    while i < nline:
        s = lines[i].lstrip()
        if not s.startswith('#'):
            i += 1
            continue
        mE = energy_re.search(s)
        if not mE:
            i += 1
            continue
        E = float(mE.group(1))

        mN = n0_re.search(s)
        n0 = int(mN.group(1)) if mN else natoms // 2

        i += 1
        taken, types, pos = 0, [], []
        while i < nline and taken < natoms:
            t = lines[i].strip()
            if t == "" or t.isdigit():
                i += 1
                continue
            parts = t.split()
            if len(parts) >= 4:
                types.append(parts[0])
                pos.append([float(x) for x in parts[1:4]])
                taken += 1
            i += 1
        if taken == natoms:
            Es.append(E); Ts.append(types); Ps.append(pos); N0s.append(n0)

    return np.array(Es), Ts, np.array(Ps, dtype=float), np.array(N0s)

def read_scan_atomicutils(fname):
    """
    Read a packed XYZ using pyBall.atomicUtils.scan_xyz() and extract
    energies from the comment lines. Returns (Es [n], Ps [n,nat,3]).
    """
    # Support both old '# E = ...' and new 'Etot ...' formats
    energy_re = re.compile(r"(?:#\s*E\s*=\s*|\bE_?tot\s+)([+-]?(?:\d*\.\d+|\d+)(?:[Ee][+-]?\d+)?)")

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


def compute_ra_vec(P, h=None, signed=True):
    """Compute r and a from the moving molecule's first atom.
    Convention:
    - Pivot is the first atom of molecule A (index 0).
    - Moving point is the first atom of molecule B (index h).
    - r = |B0 - A0|
    - a = angle in degrees (detects primary plane of scan)
    """
    n, nat, _ = P.shape
    if h is None:
        h = nat // 2
    A0 = P[:, 0, :]
    B0 = P[:, h, :]
    V  = B0 - A0

    r = np.linalg.norm(V, axis=1)

    # Detect which plane to use for angle (the one with most variation)
    vx = V[:, 0]
    vy = V[:, 1]
    vz = V[:, 2]

    std_x = np.std(vx)
    std_y = np.std(vy)
    std_z = np.std(vz)

    # Default to XZ if it looks like the scan plane
    if std_y < min(std_x, std_z) * 0.1:
        a = np.degrees(np.arctan2(vz, vx))
    elif std_z < min(std_x, std_y) * 0.1:
        a = np.degrees(np.arctan2(vy, vx))
    else:
        # Fallback to whatever has most variation vs the largest component
        # This is a bit heuristic but better than hardcoded XZ
        a = np.degrees(np.arctan2(vz, np.sqrt(vx**2 + vy**2)))

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

def plot_imshow(V, rv, A, emin=None, vmax=None, title=None, cmap='bwr', kcal=False, ax=None, bColorbar=True, rtick_step=5, bSym=False, bByMin=False):
    fac = 23.060548 if kcal else 1.0

    print(f"plot_imshow title({title}) V.shape", V.shape)
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
    if vmax is None:
        if bSym:
            if bByMin:
                vmin = np.nanmin(V)
                vmax = -vmin
            else:
                vmax = np.nanmax(np.abs(V))
                vmin = -vmax
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
