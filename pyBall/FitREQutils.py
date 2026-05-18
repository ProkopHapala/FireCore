"""
Utilities for FitREQ energy map analysis, XYZ file processing, and electron pair handling.

Data I/O and file operations:
- read_xyz_data: Extract energy and distance data from concatenated XYZ comment lines
- mark_molecule_blocks: Identify molecule boundaries in multi-molecule XYZ files
- slice_and_reshape: Slice and reshape data arrays based on molecule block marks
- concatenate_xyz_files: Concatenate XYZ files from multiple directories
- concatenate_xyz_files_flat: Concatenate XYZ files from a flat list
- find_all_dirs: Find all subdirectories in a base path
- read_file_comments: Extract comment lines from a file
- extract_comments_and_types: Extract comments and atom type counts from XYZ files
- parse_xyz_blocks: Parse XYZ file into blocks of (es, apos, qs, comment)
- parse_xyz_with_headers: Parse XYZ file with header (n0, Etot, x0, angle) extraction
- read_scan_atomicutils: Read scan data using atomicUtils.load_xyz_movie
- parse_panel_list: Parse list file for batch plotting

Grid manipulation and reshaping:
- shift_grid: Shift grid values to make minimum zero
- reshape_to_grid: Reshape 1D arrays into 2D grid based on row detection
- compute_ref_shift: Compute reference shift from energy and distance arrays
- compute_shift_from_grid: Compute shift from 2D energy grid

Energy analysis and extraction:
- extract_min_curves: Extract minimum energy curves along distance axis
- compute_min_lines_from_panel: Compute minimum energy lines from panel data
- _distances_from_Xpanel: Extract distances from Xpanel array

Plotting utilities:
- plot_Epanels: Plot energy panels from multiple reference directories
- plot_Epanels_diff: Plot energy panel differences between model and reference
- plot2Dlong_diff: Plot 2D energy difference as long format
- plot_Epanels_diff_separate: Plot separate energy panel differences
- plot_energy_2d_from_xyz: Plot 2D energy map directly from XYZ file
- plot_compare: Compare reference and model energy grids
- plot_compare_combined: Combined comparison plot with multiple subplots
- plot_min_lines_pair: Plot minimum energy lines for reference and model
- plot_imshow: Plot 2D energy grid as imshow
- plot_polar: Plot 2D energy grid in polar coordinates
- plot_profiles: Plot energy profiles along distance axis
- plot_list: Batch plotting from a list file
- plot_molecule: Plot molecule geometry with atom colors and sizes
- plot_system_panel: Plot system energy panel with Rmin overlay
- plot_energy_panel: Plot energy minimum curve panel

XYZ parsing and geometry analysis:
- _parse_comment: Parse XYZ comment line for energy and geometry info
- parse_xyz_mapping: Parse XYZ file to extract distance/angle mapping
- derive_ra_from_block: Derive distance and angle from atomic positions
- compute_ra_vec: Compute distance and angle vectors between fragments
- parse_headers_ra: Parse comment headers to extract energy, distance, angle
- detect_rows_by_r: Detect row boundaries based on distance values

Frame grid building:
- compute_panel_data: Compute panel data from XYZ file for plotting
- _build_frame_grid: Build frame grid and extract energy/geometry data

Element helpers:
- _element_color: Get color for atom element
- _element_size: Get size for atom element

Frame reordering:
- reorder_frames_by_angle: Reorder XYZ frames by angle then distance

Electron pair processing:
- AtomType: Data class for atom type properties
- read_atom_types: Read AtomTypes.dat file
- parse_n0_from_comment: Parse n0 (fragment split) from comment
- update_n0_in_comment: Update n0 value in comment string
- element_from_typename: Extract base element name from type name
- count_epairs_from_type: Get number of electron pairs for atom type
- has_sigma_hole_type: Check if atom type is a sigma-hole donor
- _z_to_qs: Get atomic charge from atomic number
- _z_to_rs: Get atomic radius from atomic number
- build_fragment: Build AtomicSystem for one fragment with epairs/sigma holes
- process_frame: Process single XYZ frame to add epairs/sigma holes
- process_one_file: Process entire XYZ file to add epairs/sigma holes
- find_data_path: Find AtomTypes.dat data directory

File I/O (saving):
- save_grid_npz: Save energy grid to NumPy format
- save_grid_gnuplot: Save energy grid to gnuplot format
- save_min_lines_npz: Save minimum energy curves to NumPy format
- save_min_lines_gnuplot: Save minimum energy curves to gnuplot format
- save_data: Save comparison data in multiple formats
"""
from pathlib import Path
import numpy as np
import os
import re
import matplotlib.pyplot as plt

from pyBall.atomicUtils import scan_xyz, load_xyz_movie, normalize, writeToXYZ
from pyBall import elements
from pyBall.AtomicSystem import AtomicSystem
import pyBall.atomicUtils as au

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
    """Slice and reshape data arrays based on molecule block marks.
    
    Takes energy arrays and block marks to slice them into per-molecule
    chunks, then reshapes based on angle data for 2D grid representation.
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
    """Concatenate XYZ files from multiple directories into single output.
    
    Finds specified file in each subdirectory and concatenates them.
    Useful for combining results from multiple calculations.
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
    """Concatenate XYZ files from a flat list of filenames.
    
    Simpler version for when files are in a single directory.

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
    """Find all subdirectories in a base path.
    
    Return all subdirectories containing an all.xyz file."""
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
    """Shift grid values to make minimum zero by subtracting asymptotic baseline.
    
    Subtracts the minimum value at the largest distance from the entire grid
    to set the baseline to zero, useful for energy landscape visualization.
    Returns shifted grid, reference value, and minimum location.
    """
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
    """Extract minimum energy curves along distance axis.
    
    For each angle, finds the distance and energy at the minimum.
    Returns arrays of rmin and emin for each angle, optionally filtered by rmax.
    G: (n_distances, n_angles) energy grid.
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
    """Compute minimum energy lines from panel data.
    
    Extracts minimum energy curves from 2D energy panel data,
    optionally shifting to zero baseline and filtering by rmax.
    Returns rmin, emin arrays.
    """
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
    """Extract distances from Xpanel array.
    
    Xpanel contains distance values for each angle/distance point.
    Returns unique distance values.
    """
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
    """Save energy grid to NumPy format.
    
    Stores angles, distances, and energy grid in .npz file
    for later loading and analysis.
    """
    np.savez(filepath, angles=angles, distances=distances, grid=grid)

def save_grid_gnuplot(angles, distances, grid, filepath):
    """Save energy grid to gnuplot format.
    
    Writes grid data in format suitable for gnuplot plotting.
    """
    with open(filepath, 'w') as f:
        f.write("# angle distance value\n")
        ny, nx = grid.shape
        for i in range(ny):
            for j in range(nx):
                v = grid[i, j]
                if not np.isnan(v):
                    f.write(f"{angles[j]:8.1f} {distances[i]:8.3f} {v:15.8e}\n")

def save_min_lines_npz(angles, rmin, emin, filepath):
    """Save minimum energy curves to NumPy format.
    
    Stores angles, Rmin distances, and minimum energies.
    """
    np.savez(filepath, angles=angles, rmin=rmin, emin=emin)

def save_min_lines_gnuplot(angles, rmin, emin, filepath):
    """Save minimum energy curves to gnuplot format.
    
    Writes min curve data in format suitable for gnuplot.
    """
    with open(filepath, 'w') as f:
        f.write("# angle rmin emin\n")
        for i in range(len(angles)):
            if not np.isnan(rmin[i]):
                f.write(f"{angles[i]:8.1f} {rmin[i]:8.3f} {emin[i]:15.8e}\n")

# ─────────────────────────────────────────────────────────
# Plotting functions
# ─────────────────────────────────────────────────────────

def plot_Epanels(Eplots, ref_dirs, bColorbar=True, Emin=-5.0, bKcal=False):
    """Plot energy panels from multiple reference directories.
    
    Creates a row of subplots showing energy panels for each reference directory.
    Useful for comparing energy landscapes across different calculations.
    """
    E_units = 1.0
    if bKcal:
        E_units = ev2kcal
    nmols = len(Eplots)
    if nmols != len(ref_dirs):
        print(f"Error: len(Eplots={nmols}) != len(ref_dirs={len(ref_dirs)})")
        return
    fig, axs = plt.subplots(1, nmols, figsize=(20, 3))
    for i in range(nmols):
        im = axs[i].imshow(Eplots[i].T * E_units, aspect='auto', origin='lower', vmin=Emin, vmax=-Emin, cmap='bwr')
        axs[i].set_title(f"Ref: {ref_dirs[i]}")
        axs[i].set_ylabel('Reference Energies')
        if bColorbar:
            plt.colorbar(im, ax=axs[i])
    plt.tight_layout()
    return fig


def plot2Dlong_diff(Erefs, Es, lens):
    """Plot 2D energy difference as long format.
    
    Plots reference vs model energies as scatter plot for comparison.
    """
    plt.figure(figsize=(10, 6))
    for i, (eref, emod, l) in enumerate(zip(Erefs, Es, lens)):
        plt.plot(eref, emod, 'o', label=f'len={l}')
    plt.xlabel('Reference Energy')
    plt.ylabel('Model Energy')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    return plt.gcf()


def plot_Epanels_diff(Emodels, Erefs, ref_dirs, bColorbar=True, Emin=-5.0, bKcal=False):
    """Plot energy panel differences between model and reference.
    
    For single molecule: shows ref, model, and diff in 3 rows.
    For multiple molecules: shows diff in single row.
    """
    E_units = 1.0
    if bKcal:
        E_units = ev2kcal
    nmols = len(Erefs)
    if nmols == 1:
        fig, axs = plt.subplots(3, 1, figsize=(6, 9))
        for i in range(nmols):
            im0 = axs[0].imshow(Erefs[i].T * E_units, aspect='auto', origin='lower',  vmin=Emin, vmax=-Emin, cmap='bwr')
            axs[0].set_title(f"Ref: {ref_dirs[i]}")
            if bColorbar: plt.colorbar(im0, ax=axs[0])

            im1 = axs[1].imshow(Emodels[i].T * E_units, aspect='auto', origin='lower', vmin=Emin, vmax=-Emin, cmap='bwr')
            axs[1].set_title(f"Model: {ref_dirs[i]}")
            if bColorbar: plt.colorbar(im1, ax=axs[1])

            Ediff = (Emodels[i] - Erefs[i]) * E_units
            dmax = np.nanmax(np.abs(Ediff))
            im2 = axs[2].imshow(Ediff.T, aspect='auto', origin='lower', vmin=-dmax, vmax=dmax, cmap='bwr')
            axs[2].set_title(f"Diff: {ref_dirs[i]}")
            if bColorbar: plt.colorbar(im2, ax=axs[2])
    else:
        fig, axs = plt.subplots(3, nmols, figsize=(6 * nmols, 9))
        for i in range(nmols):
            im0 = axs[0, i].imshow(Erefs[i].T * E_units, aspect='auto', origin='lower',  vmin=Emin, vmax=-Emin, cmap='bwr')
            axs[0, i].set_title(f"Ref: {ref_dirs[i]}")
            if bColorbar: plt.colorbar(im0, ax=axs[0, i])

            im1 = axs[1, i].imshow(Emodels[i].T * E_units, aspect='auto', origin='lower', vmin=Emin, vmax=-Emin, cmap='bwr')
            axs[1, i].set_title(f"Model: {ref_dirs[i]}")
            if bColorbar: plt.colorbar(im1, ax=axs[1, i])

            Ediff = (Emodels[i] - Erefs[i]) * E_units
            dmax = np.nanmax(np.abs(Ediff))
            im2 = axs[2, i].imshow(Ediff.T, aspect='auto', origin='lower', vmin=-dmax, vmax=dmax, cmap='bwr')
            axs[2, i].set_title(f"Diff: {ref_dirs[i]}")
            if bColorbar: plt.colorbar(im2, ax=axs[2, i])
    plt.tight_layout()
    return fig


def plot2Dlong_diff(Erefs, Es, lens):
    """Plot 2D energy difference as long format with ref and diff panels.
    
    Creates 2-row plot showing reference energies and differences.
    """
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


def plot_Epanels_diff_separate(Emodels, Erefs, ref_dirs, save_prefix=None, bColorbar=True, Emin=-5.0, bKcal=False, bClose=False):
    """Plot separate energy panel differences for each reference directory.
    
    Creates individual files for each molecule showing ref, model, and diff.
    Useful for detailed per-molecule analysis.
    """
    E_units = 1.0
    if bKcal:
        E_units = ev2kcal
    nmols = len(Erefs)
    figs = []
    for i in range(nmols):
        fig, axs = plt.subplots(1, 3, figsize=(18, 4))
        im0 = axs[0].imshow(Erefs[i].T * E_units, aspect='auto', origin='lower', vmin=Emin, vmax=-Emin, cmap='bwr')
        axs[0].set_title(f"Ref: {ref_dirs[i]}")
        if bColorbar: plt.colorbar(im0, ax=axs[0])

        im1 = axs[1].imshow(Emodels[i].T * E_units, aspect='auto', origin='lower',  vmin=Emin, vmax=-Emin, cmap='bwr')
        axs[1].set_title(f"Model: {ref_dirs[i]}")
        if bColorbar: plt.colorbar(im1, ax=axs[1])

        Ediff = (Emodels[i] - Erefs[i]) * E_units
        dmax = np.nanmax(np.abs(Ediff))
        im2 = axs[2].imshow(Ediff.T, aspect='auto', origin='lower',  vmin=-dmax, vmax=dmax, cmap='bwr')
        axs[2].set_title(f"Diff: {ref_dirs[i]}")
        if bColorbar: plt.colorbar(im2, ax=axs[2])
        plt.tight_layout()
        if save_prefix:
            plt.savefig(f"{save_prefix}_{i}.png")
        if bClose:
            plt.close(fig)
        figs.append(fig)
    return figs


def plot_energy_2d_from_xyz(xyz_path, distances=None, angles=None, title=None,  cmap='bwr', vmin=None, vmax=None, save_path=None):
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
    im = ax.imshow(GS, origin='lower', aspect='auto', cmap=cmap, vmin=vmin, vmax=vmax, interpolation='nearest',  extent=extent)
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
    """Parse XYZ comment line for energy and geometry info.
    
    Extracts Etot, x0 (distance), and angle (y or z) from comment string.
    """
    # comment format: "# n0 N Etot E x0 R y A" or "# n0 N Etot E x0 R z A"
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
    """Parse XYZ file to extract distance/angle mapping.
    
    Reads scan data and extracts the mapping between frames and their
    geometric parameters (distance and angle).
    """
    Es, Ps = read_scan_atomicutils(xyz_path)
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


def plot_compare(Gref, Gmodel, angles, distances, title, save_prefix=None,vmin=None, vmax=None, line=False, kcal=False,  save_data_prefix=None, save_fmt="both"):
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
        plot_min_lines_pair(EpanelR, EpanelM, Xpanel, angles, title=title,to_kcal=False, save_fmt=save_fmt)

    fig, axs = plt.subplots(3, 1, figsize=(8, 10))
    extent = [angles[0], angles[-1], distances[0], distances[-1]]

    im0 = axs[0].imshow(GRS, origin='lower', aspect='auto', cmap='bwr', vmin=vmin, vmax=vmax, interpolation='nearest', extent=extent)
    axs[0].set_title(f"Reference ({unit})")
    plt.colorbar(im0, ax=axs[0])

    im1 = axs[1].imshow(GMS, origin='lower', aspect='auto', cmap='bwr', vmin=vmin, vmax=vmax, interpolation='nearest', extent=extent)
    axs[1].set_title(f"Model ({unit})")
    plt.colorbar(im1, ax=axs[1])

    D = GMS - GRS
    dmax = max(-np.nanmin(D), np.nanmax(D))
    im2 = axs[2].imshow(D, origin='lower', aspect='auto', cmap='bwr',  vmin=-dmax, vmax=dmax, interpolation='nearest', extent=extent)
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


def plot_compare_combined(Gref, Gmodel, angles, distances, title, save_path=None, kcal=False, params_text=None):
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
    im00 = ax00.imshow(GRS, origin='lower', aspect='auto', cmap='bwr',  vmin=vmin, vmax=vmax, interpolation='nearest', extent=extent)
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
    im10 = ax10.imshow(GMS, origin='lower', aspect='auto', cmap='bwr',  vmin=vmin, vmax=vmax, interpolation='nearest', extent=extent)
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
    im20 = ax20.imshow(D, origin='lower', aspect='auto', cmap='bwr',  vmin=-dmax, vmax=dmax, interpolation='nearest', extent=extent)
    ax20.set_title(f"Difference ({unit})")
    ax20.set_xlabel("Angle (deg)")
    ax20.set_ylabel("Distance (A)")
    plt.colorbar(im20, ax=ax20)

    # Text panel
    ax21 = fig.add_subplot(gs[2, 1])
    ax21.axis('off')
    if params_text:
        ax21.text(0.05, 0.95, params_text, transform=ax21.transAxes, fontsize=10, verticalalignment='top', fontfamily='monospace')

    plt.suptitle(title)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    return fig


def compute_min_lines_from_panel(Epanel, Xpanel, angles, rmax=None, do_shift=True):
    """Compute minimum energy lines from panel data.
    
    Extracts minimum energy curves from 2D energy panel data,
    optionally shifting to zero baseline and filtering by rmax.
    Returns rmin, emin arrays.
    """
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


def plot_min_lines_pair(Epanel_ref, Epanel_mod, Xpanel, angles, title=None, save_path=None, to_kcal=False, ms=2, lw=0.5, save_data_prefix=None, save_fmt="both"):
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
    Parse comment headers to extract per-block Etot (energy), radius (x0) and angle (y or z).
    Supports both old energy-only headers and new enriched headers.
    The angle field can be 'y' or 'z'; both are tried and the non-zero value is preferred.
    Returns (Eh, Rh, Ah) as float arrays; values may be NaN if missing.
    """
    Eh = []
    Rh = []
    Ah = []
    _num = r'([+-]?(?:\d*\.\d+|\d+)(?:[Ee][+-]?\d+)?)'
    reE = re.compile(r'(?:#\s*E\s*=\s*|\bEtot\s+)' + _num)
    reR = re.compile(r'\bx0\s+' + _num)
    reAy = re.compile(r'\by\s+' + _num)
    reAz = re.compile(r'\bz\s+' + _num)
    with open(fname, 'r') as f:
        for ln in f:
            s = ln.lstrip()
            if not s.startswith('#'):
                continue
            mE = reE.search(s)
            mR = reR.search(s)
            mAy = reAy.search(s)
            mAz = reAz.search(s)
            Eh.append(float(mE.group(1)) if mE else np.nan)
            Rh.append(float(mR.group(1)) if mR else np.nan)
            # prefer the non-zero angle field; if both present pick the one != 0
            ay = float(mAy.group(1)) if mAy else np.nan
            az = float(mAz.group(1)) if mAz else np.nan
            if np.isfinite(ay) and ay != 0:
                Ah.append(ay)
            elif np.isfinite(az):
                Ah.append(az)
            elif np.isfinite(ay):
                Ah.append(ay)
            else:
                Ah.append(np.nan)
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
    """Compute reference shift from energy and distance arrays.
    
    Calculates baseline shift from the last distance point in each row.
    Used to shift energy grids to zero baseline.

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
    """Plot 2D energy grid as imshow.
    
    Creates a 2D heatmap of energy values with distance and angle axes.
    Supports various styling options like colormap, units, and symmetry.
    """
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
    """Plot 2D energy grid in polar coordinates.
    
    Creates polar plot with angle as theta and distance as radius.
    Useful for visualizing directional energy landscapes.
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

def plot_polar_symmetric(V, rv, A, title=None, cmap='seismic', kcal=False, ax=None, bColorbar=True, rmax=6.0, R=None, vmin=None, vmax=None, geometry=None, n0=None, plane='xy'):
    """Plot 2D energy grid in polar coordinates with symmetric color scale and fixed Rcut.
    
    Creates polar plot with:
    - Symmetric color scale (vmax = -vmin) for balanced visualization
    - Real distances and angles from grid
    - Fixed radial limit at Rcut = 6.0A to focus on relevant interaction range
    - Optional geometry overlay for first fragment (atoms and electron pairs)
    
    Args:
        V: 2D energy grid V[angle, distance]
        rv: Distance values (1D array)
        A: Angle values in degrees (1D array)
        title: Plot title
        cmap: Colormap name (default: 'bwr')
        kcal: If True, convert energies to kcal/mol
        ax: Matplotlib axes (if None, creates new polar subplot)
        bColorbar: If True, add colorbar
        rmax: Maximum radial distance in Angstrom (default: 6.0)
        R: Optional per-row distance grid (if None, uses rv)
        vmin: Minimum color limit (if None, computed from data)
        vmax: Maximum color limit (if None, computed from data)
        geometry: Optional tuple (enames, apos) to plot fragment A geometry
        n0: Optional split index for geometry (fragment A atoms are indices < n0)
        plane: Which plane to project geometry onto ('xy' or 'xz', default: 'xy')
    """
    fac = 23.060548 if kcal else 1.0
    # Sort rows by angle to make contours well-behaved
    order = np.argsort(A)
    A = A[order]
    V = V[order, :]
    if R is not None:
        R = R[order, :]
    Z = V * fac
    # Build theta (radians), shift by -pi to show -90..+90 within 90..270 window
    thetas = np.radians(A) - np.pi
    # Build coordinate grids matching V's shape. Use per-row R if available to avoid NaN padding in rv
    ny, nx = V.shape
    if R is not None:
        Rg = R
        Tg = np.repeat(thetas[:, None], nx, axis=1)
    else:
        Rg, Tg = np.meshgrid(rv, thetas)
    # Replace NaN in coordinate grids with large value (outside plot range) to avoid pcolormesh error
    Rg = np.nan_to_num(Rg, nan=rmax * 2)
    Tg = np.nan_to_num(Tg, nan=0.0)
    # Symmetric color scale: if vmin/vmax not provided, compute from data
    if vmin is None or vmax is None:
        zmax = np.nanmax(np.abs(Z))
        vmin = -zmax
        vmax = +zmax
    if ax is None:
        ax = plt.subplot(111, projection='polar')
    # Plot with symmetric color scale - clip values and replace NaN to show oversaturation instead of white spots
    Zp = np.clip(np.nan_to_num(Z, nan=vmin), vmin, vmax)
    # If everything masked (all NaN), skip plotting to avoid matplotlib errors
    if not np.any(np.isfinite(Zp)):
        ax.text(0.5, 0.5, 'No finite data', transform=ax.transAxes, ha='center', va='center')
        if title: ax.set_title(title)
        return ax
    # Use pcolormesh instead of contourf to handle irregular grids and avoid white spots
    cs = ax.pcolormesh(Tg, Rg, Zp, cmap=cmap, vmin=vmin, vmax=vmax, shading='auto')
    # Half circle setup: 90..270 window
    ax.set_thetamin(90)
    ax.set_thetamax(270)
    # Fixed radial limit at Rcut
    ax.set_ylim(0.0, rmax)
    # Add radial grid lines for distance reference (1, 2, 3, 4, 5 Angstrom)
    ax.grid(True, linestyle='--', linewidth=0.5, alpha=0.5)
    # Plot geometry if provided
    if geometry is not None and len(geometry) == 2:
        enames, apos = geometry
        # Get fragment A atoms (indices < n0 if n0 provided, else all)
        if n0 is not None and n0 > 0:
            apos_A = apos[:n0]
            enames_A = enames[:n0]
        else:
            apos_A = apos
            enames_A = enames
        
        # Convert to polar coordinates based on specified plane
        if plane == 'xz':
            # Use xz plane for z-scans
            r_geo = np.sqrt(apos_A[:, 0]**2 + apos_A[:, 2]**2)
            theta_geo = np.arctan2(apos_A[:, 2], apos_A[:, 0]) - np.pi
        else:
            # Default to xy plane for y-scans
            r_geo = np.sqrt(apos_A[:, 0]**2 + apos_A[:, 1]**2)
            theta_geo = np.arctan2(apos_A[:, 1], apos_A[:, 0]) - np.pi  # Shift by -pi to match plot
        
        # Plot atoms with colors by element
        for i, (ename, r, theta) in enumerate(zip(enames_A, r_geo, theta_geo)):
            elem = ename.split('_')[0]
            color = _element_color(elem)
            size = _element_size(elem) if elem != 'E' else 15  # Smaller for electron pairs
            ax.scatter(theta, r, c=color, s=size, zorder=10, alpha=0.8, edgecolors='black', linewidth=0.5)
    
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
    """Compute panel data from XYZ file for plotting.
    
    Parses XYZ file, extracts energies and geometry, and builds
    2D grid data suitable for panel plotting.
    """
    Es, Ps = read_scan_atomicutils(xyz)
    Eh, Rh, Ah = parse_headers_ra(xyz)
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
    """Batch plotting from a list file.
    
    Reads a list file and creates a grid of plots for multiple XYZ files.
    Supports various plot types (imshow, polar, profiles) and styling.
    """
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

    # Extract n0 (fragment split) and header-derived r/a from comment lines
    Eh, Rh, Ah = parse_headers_ra(fp)
    _, _, _, N0s = parse_xyz_with_headers(fp)
    h = int(N0s[0]) if len(N0s) > 0 and N0s[0] > 0 else None

    # Use header energy values if available and Es is all NaN (fallback from read_scan_atomicutils)
    if np.sum(np.isfinite(Es)) == 0 and np.sum(np.isfinite(Eh)) > 0:
        Es = Eh[:len(Es)] if len(Eh) >= len(Es) else Eh

    # Use header values for r (x0) and a (y or z) - mandatory, no fallback
    n = len(Es)
    if np.sum(np.isfinite(Rh)) != n: raise ValueError(f"Header distance (x0) values not available for all {n} samples")
    if np.sum(np.isfinite(Ah)) != n: raise ValueError(f"Header angle (y or z) values not available for all {n} samples")
    r = Rh[:n]
    a = Ah[:n]

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
    """Get color for atom element.
    
    Returns matplotlib color for element based on CPK coloring scheme.
    """
    # CPK coloring: O red, N blue, C gray, H white, etc.
    if elem == 'E':
        return '#9932CC'  # purple
    try:
        return elements.ELEMENT_DICT[elem][8]
    except (KeyError, IndexError):
        return '#808080'


def _element_size(elem):
    """Get size for atom element.
    
    Returns marker size for atom based on element type.
    Dummy atoms (E) get smaller size.
    """
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
        ax.text(0.5, 0.5, "No geometry", ha='center', va='center',  transform=ax.transAxes)
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

    ax.scatter(apos[:, 0], apos[:, 1], c=colors, s=sizes, zorder=5)
    # Labels: 1-based index, numbered per fragment if n0 given
    if n0 is not None and 0 < n0 < nat:
        for i in range(n0):
            ax.annotate(str(i + 1), apos[i, :2], textcoords="offset points", xytext=(0, 6), fontsize=6, ha='center', color='blue',fontweight='bold')
        for i in range(n0, nat):
            ax.annotate(str(i + 1 - n0) + "'", apos[i, :2],textcoords="offset points", xytext=(0, 6), fontsize=6, ha='center', color='red',fontweight='bold')
        # Separator line
        mid_x = (apos[:n0, 0].max() + apos[n0:, 0].min()) / 2
        ax.axvline(x=mid_x, color='gray', ls=':', lw=1, alpha=0.5)
    else:
        for i in range(nat):
            ax.annotate(str(i + 1), apos[i, :2], textcoords="offset points",xytext=(0, 6), fontsize=6, ha='center',color='black', fontweight='bold')
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


def plot_system_panel(V, rv, A, ax, label, kcal, sym, overlay_rmin=False, cmap='seismic'):
    """Row 1: 2D imshow + Rmin overlay + global min marker."""
    if not np.any(np.isfinite(V)):
        ax.set_axis_off()
        return None
    im = plot_imshow(V, rv, A, title=label, cmap=cmap, kcal=kcal, ax=ax, bColorbar=True, rtick_step=5, bSym=False)
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
                    ax.plot(As, y0 + (ps + 0.5) * (y1 - y0) / nx,  'k-', lw=1.5, alpha=0.8)
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

def reorder_frames_by_angle(inpath, outpath=None):
    """Reorder XYZ frames by angle (y or z field) then distance (x0 field).
    
    Reads the file, parses angle and distance from comment headers,
    sorts by angle ascending, then distance ascending, and writes back.
    """
    _num = r'([+-]?(?:\d*\.\d+|\d+)(?:[Ee][+-]?\d+)?)'
    reE = re.compile(r'(?:#\s*E\s*=\s*|\bEtot\s+)' + _num)
    reR = re.compile(r'\bx0\s+' + _num)
    reAy = re.compile(r'\by\s+' + _num)
    reAz = re.compile(r'\bz\s+' + _num)
    reN0 = re.compile(r'\bn0\s+(\d+)')
    
    frames = []
    with open(inpath, 'r') as f:
        lines = f.readlines()
    
    i = 0
    nline = len(lines)
    while i < nline:
        # Read natoms line
        if not lines[i].strip().isdigit():
            i += 1
            continue
        natoms = int(lines[i].strip())
        i += 1
        if i >= nline: break
        
        # Read comment line
        if not lines[i].lstrip().startswith('#'):
            i += 1
            continue
        comment = lines[i].strip()
        mE = reE.search(comment)
        if not mE:
            i += 1
            continue
        
        # Parse angle and distance
        mR = reR.search(comment)
        mAy = reAy.search(comment)
        mAz = reAz.search(comment)
        
        r = float(mR.group(1)) if mR else 0.0
        ay = float(mAy.group(1)) if mAy else np.nan
        az = float(mAz.group(1)) if mAz else np.nan
        # Prefer non-zero angle, default to 0
        a = ay if (np.isfinite(ay) and ay != 0) else (az if np.isfinite(az) else 0.0)
        
        i += 1
        
        # Read atom lines
        atom_lines = []
        taken = 0
        while i < nline and taken < natoms:
            line = lines[i].strip()
            if line and not line.isdigit():
                atom_lines.append(line)
                taken += 1
            i += 1
        
        if taken == natoms:
            frames.append({
                'comment': comment,
                'natoms': natoms,
                'atoms': atom_lines,
                'angle': a,
                'distance': r
            })
    
    # Sort by angle, then distance
    frames.sort(key=lambda f: (f['angle'], f['distance']))
    
    # Write back
    if outpath is None:
        outpath = inpath
    with open(outpath, 'w') as f:
        for frame in frames:
            f.write(str(frame['natoms']) + '\n')
            f.write(frame['comment'] + '\n')
            for atom_line in frame['atoms']:
                f.write(atom_line + '\n')
    
    print(f"Reordered {len(frames)} frames: {inpath} -> {outpath}")
    return len(frames)

# ── Data structures for AtomTypes.dat / ElementTypes.dat ──

class AtomType:
    def __init__(self, name, parent_name="*", element_name="", epair_name="",
                 valence=0, nepair=0, npi=0, sym=0,
                 Ruff=0.0, RvdW=0.0, EvdW=0.0, Qbase=0.0, Hb=0.0):
        self.name = name
        self.parent_name = parent_name
        self.element_name = element_name
        self.epair_name = epair_name
        self.valence = valence
        self.nepair = nepair
        self.npi = npi
        self.sym = sym
        self.Ruff = Ruff
        self.RvdW = RvdW
        self.EvdW = EvdW
        self.Qbase = Qbase
        self.Hb = Hb


def read_atom_types(filepath):
    """Read AtomTypes.dat file and return dictionary of AtomType objects."""
    atom_types = {}
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:
                continue
            parts = line.split()
            if len(parts) < 5:
                continue
            name = parts[0]
            at = AtomType(name=name, parent_name=parts[1], element_name=parts[2], epair_name=parts[3])
            try:
                at.valence = int(parts[4])
                at.nepair = int(parts[5])
                at.npi = int(parts[6])
                at.sym = int(parts[7])
                at.Ruff = float(parts[8])
                at.RvdW = float(parts[9])
                at.EvdW = float(parts[10])
                at.Qbase = float(parts[11])
                at.Hb = float(parts[12])
            except (ValueError, IndexError):
                continue
            atom_types[name] = at
    return atom_types


# ── Core logic ──

def parse_n0_from_comment(comment):
    m = re.search(r"\bn0\s+(\d+)", comment)
    return int(m.group(1)) if m else None


def update_n0_in_comment(comment, new_n0):
    return re.sub(r"\bn0\s+\d+", f"n0 {new_n0}", comment)


def element_from_typename(tname):
    """Extract base element name from type name (e.g., O_2 -> O)."""
    return tname.split("_")[0] if "_" in tname else tname


def count_epairs_from_type(tname, atom_types):
    """Return how many epairs this atom type should have, or 0 if unknown."""
    at = atom_types.get(tname)
    if at is None:
        return 0
    return max(0, at.nepair)


def has_sigma_hole_type(tname, atom_types):
    """Return True if this atom type is a sigma-hole donor (has epair_name, nepair==0)."""
    at = atom_types.get(tname)
    if at is None:
        return False
    return at.nepair == 0 and at.epair_name not in ("*", "", "E")


def _z_to_qs(z):
    """Get atomic charge from atomic number."""
    try:
        return elements.ELEMENTS[z - 1][9]
    except (IndexError, KeyError):
        return 0.0


def _z_to_rs(z):
    """Get atomic radius from atomic number."""
    try:
        return elements.ELEMENTS[z - 1][7]
    except (IndexError, KeyError):
        return 1.0


def build_fragment(apos, enames, qs, atom_types, lepair, sigma_dist, do_epairs, do_sigma):
    """
    Build an AtomicSystem for one fragment, add epairs and/or sigma holes.
    Returns (sys, n_added).
    """
    elem = []
    atypes_list = []
    valid = []
    for i, e in enumerate(enames):
        en = element_from_typename(e)
        if (en in elements.ELEMENT_DICT) and (en != 'E'):
            elem.append(en)
            atypes_list.append(elements.ELEMENT_DICT[en][0])
            valid.append(i)
        else:
            # Unknown type (e.g., 'E' dummy from previous run) — keep but use Z=0
            elem.append(en)
            atypes_list.append(0)
            valid.append(i)

    atypes = np.array(atypes_list, dtype=np.int32)
    elem_list = list(elem)

    sys = AtomicSystem(
        apos=apos.copy(),
        atypes=atypes,
        enames=elem_list,
        qs=qs.copy() if qs is not None else None,
        bPreinit=False,
    )
    if sys.qs is None:
        sys.qs = np.array([_z_to_qs(z) for z in sys.atypes])
    if sys.Rs is None:
        sys.Rs = np.array([_z_to_rs(z) for z in sys.atypes])
    sys.neighs()

    n_orig = len(sys.apos)
    n_added = 0

    # ── Electron pairs ──
    if do_epairs:
        # Override VALENCE_DICT in the AtomicSystem instance so it knows
        # which atoms get epairs from AtomTypes.dat (not just O/N).
        # VALENCE_DICT is a module-level dict; we temporarily replace it.
        from pyBall.AtomicSystem import VALENCE_DICT as _orig_vd

        valence_map = {}
        for i, tname in enumerate(enames):
            at = atom_types.get(tname)
            if at is not None and at.nepair > 0:
                elem_name = elem[i]
                nb = at.valence
                nsigma = len(sys.ngs[i]) if (sys.ngs is not None and i < len(sys.ngs)) else 0
                npi = nb - nsigma
                valence_map[elem_name] = (nb, at.nepair, npi)

        # Patch VALENCE_DICT so add_electron_pairs uses our data
        import pyBall.AtomicSystem as _asmod
        _backup = dict(_asmod.VALENCE_DICT)
        for ename, (nb, nep, _) in valence_map.items():
            _asmod.VALENCE_DICT[ename] = (nb, nep)

        # Need to also set npi per atom. add_electron_pairs/dd_epair uses npi from
        # difference between valence and sigma neighbors, which we already did above.
        # The VALENCE_DICT only gives (nb, nep), npi is derived.
        # So we just need the correct (nb, nep) entries.
        sys.add_electron_pairs(distance=lepair)
        _asmod.VALENCE_DICT.clear()
        _asmod.VALENCE_DICT.update(_backup)

        n_added += len(sys.apos) - n_orig

    # ── Sigma holes ──
    if do_sigma:
        for i in range(n_orig):
            if not has_sigma_hole_type(enames[i], atom_types): continue
            neighs = (list(sys.ngs[i].keys()) if (sys.ngs is not None and i < len(sys.ngs)) else [])
            if len(neighs) != 1: continue
            j = neighs[0]
            direction = au.normalize(sys.apos[i] - sys.apos[j])
            at = atom_types.get(enames[i])
            ep_tname = (at.epair_name if (at is not None) else "")
            if not ep_tname or ep_tname in ("*", "E"):
                ep_tname = "E"
            sys.place_electron_pair(  i, direction, distance=sigma_dist, ename=ep_tname, atype=200, qs=0.0, Rs=1.0, )
            n_added += 1

    # Restore original type names
    for i in range(len(enames)):
        # Keep epair placeholders ("E") to be relabeled below
        if element_from_typename(enames[i]) == 'E':
            continue
        sys.enames[i] = enames[i]

    # Assign specific epair type names (for both pre-existing and newly added epairs)
    for iep in range(len(sys.enames)):
        # Candidate epair atoms: either originally provided as 'E', or newly added as 'E'
        if sys.enames[iep] != 'E': continue
        ng = sys.ngs[iep] if (sys.ngs is not None and iep < len(sys.ngs)) else None
        if not ng: continue
        ih = next(iter(ng.keys()))
        if ih < 0 or ih >= len(enames): continue
        at = atom_types.get(enames[ih])
        if at is None: continue
        ep_tname = at.epair_name
        if not ep_tname or ep_tname == "*":  continue
        sys.enames[iep] = ep_tname

    return sys, n_added


def derive_ra_from_block(P):
    """Derive distance and angle from atomic positions.
    
    Computes distance between fragment centers and angle of approach
    from atomic coordinates.
    """
    n = len(P)
    if n < 2:
        return None, None


def process_frame(es, apos, qs, rs, comment, atom_types, lepair, sigma_dist, do_epairs, do_sigma):
    comment = comment.strip()
    n0 = parse_n0_from_comment(comment)
    if n0 is None:
        return None, None
    natoms = len(es)
    if n0 > natoms:
        return None, None

    esA, esB = es[:n0], es[n0:]
    aposA, aposB = apos[:n0].copy(), apos[n0:].copy()
    qsA = qs[:n0].copy() if qs is not None else None
    qsB = qs[n0:].copy() if qs is not None else None

    sysA, nA = build_fragment(aposA, esA, qsA, atom_types, lepair, sigma_dist, do_epairs, do_sigma)
    sysB, nB = build_fragment(aposB, esB, qsB, atom_types, lepair, sigma_dist, do_epairs, do_sigma)

    new_es = list(sysA.enames) + list(sysB.enames)
    new_apos = np.vstack([sysA.apos, sysB.apos])
    new_qs = (np.concatenate([sysA.qs, sysB.qs]) if (sysA.qs is not None and sysB.qs is not None) else None)
    new_comment = update_n0_in_comment(comment, n0 + nA)
    return (new_es, new_apos, new_qs, None, new_comment), (nA, nB)


def process_one_file(inpath, outpath, atom_types, lepair, sigma_dist, do_epairs, do_sigma, verbose, simple_names):
    """Process entire XYZ file to add epairs/sigma holes."""
    trj = au.load_xyz_movie(inpath)
    if not trj:
        return 0, 0, 0
    if verbose:
        print(f"    {len(trj)} frames")

    fout = open(outpath, "w")
    total_frames = 0
    total_A = 0
    total_B = 0
    for es, apos, qs, rs, comment in trj:
        res = process_frame(es, apos, qs, rs, comment, atom_types, lepair, sigma_dist, do_epairs, do_sigma)
        if res is None:  continue
        (new_es, new_apos, new_qs, _, new_comment), (nA, nB) = res
        out_es = [element_from_typename(e) for e in new_es] if simple_names else new_es
        au.writeToXYZ(fout, out_es, new_apos, qs=new_qs, comment=new_comment, bHeader=True)
        total_frames += 1
        total_A += nA
        total_B += nB
        if verbose: print(f"      frame -> {len(new_es)} atoms, n0→{parse_n0_from_comment(new_comment)}, +{nA}/+{nB}")
    fout.close()
    return total_frames, total_A, total_B


def find_data_path(script_dir):
    """Find AtomTypes.dat data directory."""
    candidates = [
        os.path.join(script_dir, "..", "..", "tests", "tFitREQ_PN", "data"),
        os.path.join(script_dir, "..", "..", "tests", "tFitREQ", "data"),
        os.path.join(script_dir, "data"),
    ]
    for d in candidates:
        d = os.path.abspath(d)
        if os.path.isfile(os.path.join(d, "AtomTypes.dat")):
            return d
    return None


def plot_profile_row(fig, axes, V_ref, V_model_total, V_model_hbond, V_model_eout, rv, A, 
                    frame_idx, Ps_raw, Ts_raw, n0_first, kcal, Rmax1D):
    """Plot second row with radial slice, angular slice, E_min(angle), and geometry."""
    # Use existing axes from row 2
    ax_radial = axes[1, 0]
    ax_angular = axes[1, 1]
    ax_emin = axes[1, 2]
    ax_geom = axes[1, 3]
    
    # Find global minimum indices from Reference
    iy_g, ix_g = np.unravel_index(np.nanargmin(V_ref), V_ref.shape)
    angle_min = A[iy_g]
    r_min = rv[ix_g] if ix_g < len(rv) else np.nan
    
    # 1) Radial slice at minimum angle (Eref, Etot, Ein, Eout)
    fac = 23.060548 if kcal else 1.0
    ax_radial.plot(rv, V_ref[iy_g, :] * fac, 'k:', lw=1.5, label='Eref')
    ax_radial.plot(rv, V_model_total[iy_g, :] * fac, 'r-', lw=0.5, label='Etot')
    ax_radial.plot(rv, V_model_hbond[iy_g, :] * fac, 'b-', lw=0.5, label='Ein')
    ax_radial.plot(rv, V_model_eout[iy_g, :] * fac, 'g-', lw=0.5, label='Eout')
    ax_radial.axvline(r_min, color='gray', linestyle='--', alpha=0.5, label='r_min')
    ax_radial.set_xlabel('Distance (Å)')
    ax_radial.set_ylabel('E (kcal/mol)' if kcal else 'E (eV)')
    ax_radial.set_title(f'Radial slice @ {angle_min:.1f}°')
    ax_radial.set_xlim(1.4, Rmax1D)
    ax_radial.legend(fontsize=8)
    ax_radial.grid(alpha=0.3)
    
    # 2) Angular slice at minimum radius (Eref, Etot, Ein, Eout)
    o = np.argsort(A)
    ax_angular.plot(A[o], V_ref[:, ix_g][o] * fac, 'k:', lw=1.5, label='Eref')
    ax_angular.plot(A[o], V_model_total[:, ix_g][o] * fac, 'r-', lw=0.5, label='Etot')
    ax_angular.plot(A[o], V_model_hbond[:, ix_g][o] * fac, 'b-', lw=0.5, label='Ein')
    ax_angular.plot(A[o], V_model_eout[:, ix_g][o] * fac, 'g-', lw=0.5, label='Eout')
    ax_angular.axvline(angle_min, color='gray', linestyle='--', alpha=0.5, label='angle_min')
    ax_angular.set_xlabel('Angle (deg)')
    ax_angular.set_ylabel('E (kcal/mol)' if kcal else 'E (eV)')
    ax_angular.set_title(f'Angular slice @ r={r_min:.2f}Å')
    ax_angular.legend(fontsize=8)
    ax_angular.grid(alpha=0.3)
    
    # 3) E_min(angle) using plot_energy_panel for reference
    plot_energy_panel(V_ref, rv, A, ax_emin, kcal)
    # Overlay model curves on top
    _, emin_tot = extract_min_curves(A, rv, V_model_total.T)
    _, emin_hbond = extract_min_curves(A, rv, V_model_hbond.T)
    _, emin_eout = extract_min_curves(A, rv, V_model_eout.T)
    ax_emin.plot(A[o], emin_tot[o] * fac, 'r-', lw=0.5, label='Etot')
    ax_emin.plot(A[o], emin_hbond[o] * fac, 'b-', lw=0.5, label='Ein')
    ax_emin.plot(A[o], emin_eout[o] * fac, 'g-', lw=0.5, label='Eout')
    ax_emin.legend(fontsize=8)
    
    # 4) Geometry at global minimum using plot_molecule (2D)
    iframe = frame_idx[iy_g, ix_g]
    if iframe >= 0 and iframe < len(Ps_raw):
        apos_mol = Ps_raw[iframe]
        enames_mol = (Ts_raw[iframe] if Ts_raw is not None and iframe < len(Ts_raw) else ["?"] * len(apos_mol))
        n0_mol = n0_first
        plot_molecule(ax_geom, enames_mol, apos_mol, n0=n0_mol, title=f"Geometry @ min\nFrame {iframe}")
    else:
        ax_geom.text(0.5, 0.5, 'No geometry', transform=ax_geom.transAxes, ha='center', va='center')