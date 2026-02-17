#!/usr/bin/env python3
"""
test_ptcda.py  –  AFM simulation of PTCDA using pyBall/OCL/AFM.py

Height convention (PTCDA is flat, all atoms at z=0):
  Contact distance:  R_O + R_C ≈ 1.52 + 1.7 = 3.2 Å  (probe just touching C atom)
  CO tip length:     ~3.5 Å  (tip anchor to O-atom)
  Tip anchor height: ~6-8 Å above molecular plane at start of scan
  Probe height:      6 → 0.5 Å above mol plane (nz=28 steps at 0.2 Å)

Outputs (all in tests/tAFM/):
  Fz_raw_PTCDA.npy      (nx,ny,nz)  – Fz before PP relaxation
  Fz_relax_PTCDA.npy    (nx,ny,nz)  – Fz after  PP relaxation
  FEs_raw_PTCDA.npy     (nx,ny,nz,4)
  FEs_relax_PTCDA.npy   (nx,ny,nz,4)
  probe_heights.npy     (nz,)       – probe height above mol plane [Å]

Usage:
  python test_ptcda.py                  # LJ + fixed charges (default)
  python test_ptcda.py --morse          # Morse + fixed charges
  python test_ptcda.py --qeq            # LJ + QEq charges
  python test_ptcda.py --noplot         # skip matplotlib
  python test_ptcda.py --nxy 80 80      # finer scan grid
"""

import sys, os, argparse
import numpy as np

_root = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..'))
if _root not in sys.path:
    sys.path.insert(0, _root)

from pyBall.OCL.AFM import AFMulator

# ── CLI ───────────────────────────────────────────────────────────────────────
ap = argparse.ArgumentParser()
ap.add_argument('--morse',  action='store_true')
ap.add_argument('--qeq',    action='store_true')
ap.add_argument('--noplot', action='store_true')
ap.add_argument('--nxy',    type=int, nargs=2, default=[60, 60], metavar=('NX','NY'))
# height range: probe scans from z_start to z_start + nz*dtip above molecular plane
ap.add_argument('--z_start', type=float, default=6.0,  help='probe start height above mol plane [Å]')
ap.add_argument('--z_end',   type=float, default=0.5,  help='probe end   height above mol plane [Å]')
ap.add_argument('--dz',      type=float, default=0.2,  help='z step [Å]')
args = ap.parse_args()

nz   = int(round((args.z_start - args.z_end) / args.dz)) + 1
dtip = -args.dz
tag  = f"{'morse' if args.morse else 'lj'}_{'qeq' if args.qeq else 'fixedQ'}"
print(f"=== test_ptcda.py  tag={tag}  nxy={args.nxy}  nz={nz}  dz={args.dz}  "
      f"probe_height {args.z_start:.1f}→{args.z_end:.1f} Å ===")

# ── paths ─────────────────────────────────────────────────────────────────────
xyz_path    = os.path.join(_root, 'cpp', 'common_resources', 'xyz', 'PTCDA.xyz')
params_path = os.path.join(_root, 'cpp', 'common_resources', 'ElementTypes.dat')
assert os.path.exists(xyz_path),    f"PTCDA.xyz not found: {xyz_path}"
assert os.path.exists(params_path), f"ElementTypes.dat not found: {params_path}"

# ── build AFMulator ───────────────────────────────────────────────────────────
afm = AFMulator(use_morse=args.morse)
mol = afm.load_molecule(xyz_path)
print(f"Loaded {mol.natoms} atoms  enames set={set(mol.enames)}")
if mol.qs is not None:
    print(f"  charges from .xyz: min={mol.qs.min():.3f} max={mol.qs.max():.3f} sum={mol.qs.sum():.3f}")

afm.assign_params(params_path=params_path)

if args.qeq:
    qs = afm.solve_QEq(Q_total=0.0)
    print(f"  QEq q range=[{qs.min():.3f},{qs.max():.3f}] sum={qs.sum():.4f}")

# Grid: use larger z_top to accommodate tip going up to mol+10 Å
L = afm.setup_grid(n=(80, 80, 60), margin=3.0, z_top=12.0)
print(f"Grid L={L}  n={afm.n}  mol_z={afm.mol_z:.2f} Å (kernel space)")

afm.make_forcefield()

# ── Compute scan positions explicitly (shared between raw and relaxed) ────────
mol_z  = afm.mol_z
nx_s, ny_s = args.nxy
apos = afm.atoms_arr[:,:3]
mn, mx = apos.min(axis=0), apos.max(axis=0)
x0 = mn[0] + (mx[0]-mn[0])*0.05
y0 = mn[1] + (mx[1]-mn[1])*0.05
# TIP starts at mol_z + z_start + |dpos0_z| (so probe starts at mol_z + z_start)
z0_tip = mol_z + args.z_start + abs(float(afm.dpos0[2]))
scan_p0 = np.array([x0, y0, z0_tip], dtype=np.float32)
scan_da = np.array([(mx[0]-mn[0])*0.9/max(nx_s-1,1), 0., 0.], dtype=np.float32)
scan_db = np.array([0., (mx[1]-mn[1])*0.9/max(ny_s-1,1), 0.], dtype=np.float32)
# Probe heights above mol plane for each z step
probe_heights = mol_z + args.z_start + np.arange(nz) * dtip - mol_z  # = z_start + iz*dtip
print(f"scan_p0={scan_p0}  probe_heights[0]={probe_heights[0]:.2f}  probe_heights[-1]={probe_heights[-1]:.2f} Å")
# scan extent in absolute Ang (for plotting axis labels)
x_ext = [float(scan_p0[0] - afm.mol_shift[0]),
          float(scan_p0[0] + scan_da[0]*(nx_s-1) - afm.mol_shift[0])]
y_ext = [float(scan_p0[1] - afm.mol_shift[1]),
          float(scan_p0[1] + scan_db[1]*(ny_s-1) - afm.mol_shift[1])]

# ── Raw FF scan (before PP relaxation) ───────────────────────────────────────
print("\n--- RAW SCAN (no PP relaxation) ---")
FEs_raw, _ = afm.get_raw_FE(nxy=tuple(args.nxy), nz=nz, dtip=dtip,
                              scan_p0=scan_p0, scan_da=scan_da, scan_db=scan_db)
Fz_raw = FEs_raw[:,:,:,2]
print(f"Fz_raw  stats: min={Fz_raw.min():.4f}  max={Fz_raw.max():.4f}  mean={Fz_raw.mean():.4f}")

# ── Relaxed scan (with PP relaxation) ─────────────────────────────────────────
print("\n--- RELAXED SCAN (with PP relaxation) ---")
FEs_relax, _ = afm.run_scan(nxy=tuple(args.nxy), nz=nz, dtip=dtip,
                              scan_p0=scan_p0, scan_da=scan_da, scan_db=scan_db)
Fz_relax = FEs_relax[:,:,:,2]
print(f"Fz_relax stats: min={Fz_relax.min():.4f}  max={Fz_relax.max():.4f}  mean={Fz_relax.mean():.4f}")

# ── df per-slice: -dFz/dz  (proxy for constant-height frequency shift) ────────
# Positive df = repulsive gradient; negative df = attractive gradient (bond contrast)
dz_step = abs(dtip)
df_raw   = -np.gradient(Fz_raw,   dz_step, axis=2)   # (nx,ny,nz)
df_relax = -np.gradient(Fz_relax, dz_step, axis=2)
print(f"df_relax stats: min={df_relax.min():.5f}  max={df_relax.max():.5f}")

# ── Save 3D .npy outputs ──────────────────────────────────────────────────────
np.save(f'Fz_raw_{tag}.npy',        Fz_raw)
np.save(f'Fz_relax_{tag}.npy',      Fz_relax)
np.save(f'FEs_raw_{tag}.npy',       FEs_raw)
np.save(f'FEs_relax_{tag}.npy',     FEs_relax)
np.save(f'df_relax_{tag}.npy',      df_relax)
np.save(f'df_raw_{tag}.npy',        df_raw)
np.save(f'probe_heights_{tag}.npy', probe_heights)
print(f"\nSaved: Fz_raw_{tag}.npy       shape={Fz_raw.shape}")
print(f"Saved: Fz_relax_{tag}.npy     shape={Fz_relax.shape}")
print(f"Saved: df_relax_{tag}.npy     shape={df_relax.shape}")
print(f"Saved: probe_heights_{tag}.npy  values={np.round(probe_heights,2)}")

# ── Determine ring-center pixels (most-attractive Fz at contact height) ───────
iz_contact = int(round((args.z_start - 3.2) / args.dz))  # h ≈ 3.2 Å (R_O+R_C)
iz_contact = max(0, min(nz-1, iz_contact))
# ring centers: 3 pixels with most-negative Fz_relax at contact height
Fz_contact = Fz_relax[:,:,iz_contact]
m = 4  # boundary exclusion (pixels)
Fz_search = Fz_contact.copy()
Fz_search[:m,:] = np.inf; Fz_search[-m:,:] = np.inf
Fz_search[:,:m] = np.inf; Fz_search[:,-m:] = np.inf
flat_idx = np.argsort(Fz_search.ravel())   # ascending → most negative first
ring_pixels = []
for fi in flat_idx:
    ix, iy = divmod(int(fi), ny_s)
    if Fz_search[ix, iy] == np.inf:
        continue
    # avoid duplicates closer than 5 pixels
    if all(abs(ix-px)>5 or abs(iy-py)>5 for px,py in ring_pixels):
        ring_pixels.append((ix, iy))
    if len(ring_pixels) == 3:
        break
print(f"Ring-center pixels at h={probe_heights[iz_contact]:.1f}Å: {ring_pixels}  "
      f"Fz_contact={[f'{Fz_contact[p]:.4f}' for p in ring_pixels]}")

# ── Plots ─────────────────────────────────────────────────────────────────────
if args.noplot:
    print("=== DONE (--noplot) ===")
    sys.exit(0)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

def safe_norm(data_2d, pct=99):
    """Symmetric ±vabs TwoSlopeNorm; handles all-positive/all-negative slices."""
    vabs = max(float(np.percentile(np.abs(data_2d), pct)), 1e-6)
    return TwoSlopeNorm(vmin=-vabs, vcenter=0, vmax=vabs)

def plot_grid_Fz(Fz, heights, tag, label, fname, x_ext, y_ext, ncols=7):
    """Plot grid of Fz images; each panel has its own colorscale (per-slice)."""
    nz = len(heights)
    nrows = int(np.ceil(nz / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(2.5*ncols, 2.8*nrows))
    axes = np.array(axes).reshape(nrows, ncols)
    fig.suptitle(f"PTCDA  Fz {label}  [{tag}]  (eV/Å) [per-slice scale]", fontsize=10)
    kw = dict(origin='lower', cmap='bwr', aspect='auto',
              extent=[x_ext[0], x_ext[1], y_ext[0], y_ext[1]])
    for k in range(nz):
        r, c = divmod(k, ncols)
        ax = axes[r, c]
        norm = safe_norm(Fz[:,:,k])
        vmax = norm.vmax
        im = ax.imshow(Fz[:,:,k].T, norm=norm, **kw)
        ax.set_title(f"h={heights[k]:.1f}Å ±{vmax:.2g}", fontsize=6)
        ax.tick_params(labelsize=4)
        plt.colorbar(im, ax=ax, shrink=0.8)
    for k in range(nz, nrows*ncols):
        r, c = divmod(k, ncols)
        axes[r, c].set_visible(False)
    plt.tight_layout()
    plt.savefig(fname, dpi=90, bbox_inches='tight')
    plt.close()
    print(f"Saved {fname}")

# Figure 1: Raw Fz at all heights
plot_grid_Fz(Fz_raw, probe_heights, tag, 'BEFORE relaxation',
             f'afm_Fz_raw_{tag}.png', x_ext, y_ext)

# Figure 2: Relaxed Fz at all heights
plot_grid_Fz(Fz_relax, probe_heights, tag, 'AFTER relaxation',
             f'afm_Fz_relax_{tag}.png', x_ext, y_ext)

# Figure 3: df per-slice grid (frequency-shift proxy, sharper bond contrast)
plot_grid_Fz(df_relax, probe_heights, tag, 'df=-dFz/dz AFTER relax',
             f'afm_df_relax_{tag}.png', x_ext, y_ext)

# Figure 4: Side-by-side comparison at representative heights
# Use separate row-wise colorbars: raw row gets its own scale, relax row gets its own
sel_heights = [5.0, 4.0, 3.5, 3.0, 2.5, 2.0]
sel_iz = [int(round((args.z_start - h) / args.dz)) for h in sel_heights]
sel_iz = [iz for iz in sel_iz if 0 <= iz < nz]
nsel = len(sel_iz)
fig, axes = plt.subplots(2, nsel, figsize=(3*nsel, 6))
fig.suptitle(f"PTCDA AFM before vs after PP relaxation  [{tag}]", fontsize=11)
kw_base = dict(origin='lower', cmap='bwr', aspect='auto',
               extent=[x_ext[0], x_ext[1], y_ext[0], y_ext[1]])
for col, iz in enumerate(sel_iz):
    h = probe_heights[iz]
    # Each height gets independent scale for raw and relaxed separately
    norm_raw   = safe_norm(Fz_raw[:,:,iz])
    norm_relax = safe_norm(Fz_relax[:,:,iz])
    im0 = axes[0,col].imshow(Fz_raw[:,:,iz].T,   norm=norm_raw,   **kw_base)
    im1 = axes[1,col].imshow(Fz_relax[:,:,iz].T, norm=norm_relax, **kw_base)
    axes[0,col].set_title(f"raw   h={h:.1f}Å\n±{norm_raw.vmax:.2g}eV/Å",   fontsize=7)
    axes[1,col].set_title(f"relax h={h:.1f}Å\n±{norm_relax.vmax:.2g}eV/Å", fontsize=7)
    plt.colorbar(im0, ax=axes[0,col], shrink=0.8)
    plt.colorbar(im1, ax=axes[1,col], shrink=0.8)
    for ax in [axes[0,col], axes[1,col]]:
        ax.tick_params(labelsize=5); ax.set_xlabel('x(Å)', fontsize=6)
axes[0,0].set_ylabel('y (Å)', fontsize=7); axes[1,0].set_ylabel('y (Å)', fontsize=7)
plt.tight_layout()
cmp_fname = f'afm_comparison_{tag}.png'
plt.savefig(cmp_fname, dpi=110, bbox_inches='tight')
plt.close()
print(f"Saved {cmp_fname}")

# Figure 5: Fz approach traces – pixels over ATOMS and RING CENTERS
def nearest_scan_pixel(pos_xy_kernel, scan_p0, scan_da, scan_db, nx_s, ny_s):
    ix = int(round((pos_xy_kernel[0] - scan_p0[0]) / scan_da[0]))
    iy = int(round((pos_xy_kernel[1] - scan_p0[1]) / scan_db[1]))
    return max(0, min(nx_s-1, ix)), max(0, min(ny_s-1, iy))

# Atom-position pixels
atom_xy = afm.atoms_arr[:, :2]
centroid_xy = atom_xy.mean(axis=0)
dists = np.linalg.norm(atom_xy - centroid_xy, axis=1)
for ia, lbl in [(int(np.argmin(dists)),'center-atom'),
                (int(np.argsort(dists)[len(dists)//2]),'mid-atom'),
                (int(np.argmax(dists)),'edge-atom')]:
    ix, iy = nearest_scan_pixel(atom_xy[ia], scan_p0, scan_da, scan_db, nx_s, ny_s)
    print(f"  atom trace {lbl}: pixel({ix},{iy})  xy=({atom_xy[ia,0]:.2f},{atom_xy[ia,1]:.2f})")
atom_pixels = [nearest_scan_pixel(atom_xy[ia], scan_p0, scan_da, scan_db, nx_s, ny_s)
               for ia in [int(np.argmin(dists)),
                           int(np.argsort(dists)[len(dists)//2]),
                           int(np.argmax(dists))]]
atom_labels = ['center-atom','mid-atom','edge-atom']

# Ring-center pixels (most-negative relaxed Fz at contact height)
print(f"Ring-center pixels: {ring_pixels}")
ring_labels = [f'ring#{i+1}({p[0]},{p[1]})' for i,p in enumerate(ring_pixels)]

def plot_traces(pixels, labels, title_suffix, fname):
    fig, axes = plt.subplots(2, 3, figsize=(14, 7))
    fig.suptitle(f"PTCDA  Fz traces  [{tag}]  {title_suffix}  (top=full, bottom=zoomed ±0.5)")
    for col, ((ix,iy), lbl) in enumerate(zip(pixels, labels)):
        for row, (ax, ylim) in enumerate(zip(axes[:,col], [None, (-0.5, 0.5)])):
            ax.plot(probe_heights, Fz_raw[ix,iy,:],   color='steelblue', label='raw FF', lw=1.2)
            ax.plot(probe_heights, Fz_relax[ix,iy,:], color='tomato', lw=1.2, ls='--', label='PP relax')
            ax.axhline(0, color='k', lw=0.5)
            ax.axvline(3.2, color='gray', lw=0.8, ls=':', label='R_O+R_C=3.2Å')
            ax.set_xlabel('Probe height (Å)', fontsize=8); ax.set_ylabel('Fz (eV/Å)', fontsize=8)
            ax.set_title(lbl if row==0 else f"{lbl} [zoom]", fontsize=8)
            ax.legend(fontsize=6); ax.invert_xaxis()
            if ylim:
                ax.set_ylim(ylim)
                ax.axhspan(ylim[0], 0, alpha=0.04, color='blue')
    plt.tight_layout()
    plt.savefig(fname, dpi=110, bbox_inches='tight'); plt.close()
    print(f"Saved {fname}")

plot_traces(atom_pixels,  atom_labels,  'over ATOMS',       f'afm_traces_atoms_{tag}.png')
plot_traces(ring_pixels,  ring_labels,  'over RING CENTERS', f'afm_traces_rings_{tag}.png')

# Figure 6: LJ vs Morse comparison (if both .npy files exist)
other_tag = tag.replace('morse','lj') if 'morse' in tag else tag.replace('lj','morse')
Fz_other_path = f'Fz_relax_{other_tag}.npy'
heights_other_path = f'probe_heights_{other_tag}.npy'
import os
if os.path.exists(Fz_other_path) and os.path.exists(heights_other_path):
    Fz_other   = np.load(Fz_other_path)
    hts_other  = np.load(heights_other_path)
    lbl_this   = 'Morse' if 'morse' in tag else 'LJ'
    lbl_other  = 'LJ'    if 'morse' in tag else 'Morse'
    sel_iz_cmp = [int(round((args.z_start - h) / args.dz)) for h in [4.0, 3.5, 3.0, 2.5, 2.0, 1.6]]
    sel_iz_cmp = [iz for iz in sel_iz_cmp if 0 <= iz < nz]
    nsel_cmp   = len(sel_iz_cmp)
    fig, axes  = plt.subplots(2, nsel_cmp, figsize=(3*nsel_cmp, 6))
    fig.suptitle(f"PTCDA  Fz relaxed: {lbl_this} vs {lbl_other}  [{tag[:3]}…]  (per-panel scale)", fontsize=10)
    for col, iz in enumerate(sel_iz_cmp):
        h = probe_heights[iz]
        iz_other = min(np.argmin(np.abs(hts_other - h)), Fz_other.shape[2]-1)
        for row, (Fz_arr, lbl_row) in enumerate([(Fz_relax, lbl_this), (Fz_other, lbl_other)]):
            iz_use = iz if row==0 else iz_other
            norm = safe_norm(Fz_arr[:,:,iz_use])
            im = axes[row,col].imshow(Fz_arr[:,:,iz_use].T, norm=norm, **kw_base)
            axes[row,col].set_title(f"{lbl_row} h={h:.1f}Å\n±{norm.vmax:.2g}", fontsize=7)
            plt.colorbar(im, ax=axes[row,col], shrink=0.8)
            axes[row,col].tick_params(labelsize=4)
    plt.tight_layout()
    cmp2_fname = f'afm_LJvsMorse_{tag}.png'
    plt.savefig(cmp2_fname, dpi=110, bbox_inches='tight'); plt.close()
    print(f"Saved {cmp2_fname}")

print("=== DONE ===")
