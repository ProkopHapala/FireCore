#!/usr/bin/env python3
"""
test_ptcda.py  –  AFM simulation of PTCDA using pyBall/OCL/AFM.py
Thin top-level script: LJ/Morse force field + PP relaxation via AFMulator GPU kernels.

Usage:
  python test_ptcda.py                  # LJ + fixed charges (default)
  python test_ptcda.py --morse          # Morse + fixed charges
  python test_ptcda.py --qeq            # LJ + QEq charges
  python test_ptcda.py --noplot         # skip matplotlib
"""
import sys, os, argparse
import numpy as np

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_ROOT     = os.path.realpath(os.path.join(_THIS_DIR, '..', '..'))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from pyBall.OCL.AFM import AFMulator, compute_df, plot_grid_Fz, safe_norm

# ── CLI ───────────────────────────────────────────────────────────────────────
ap = argparse.ArgumentParser()
ap.add_argument('--morse',  action='store_true')
ap.add_argument('--qeq',    action='store_true')
ap.add_argument('--noplot', action='store_true')
ap.add_argument('--nxy',    type=int, nargs=2, default=[60, 60], metavar=('NX','NY'))
ap.add_argument('--z_start', type=float, default=6.0,  help='probe start height [A]')
ap.add_argument('--z_end',   type=float, default=0.5,  help='probe end height [A]')
ap.add_argument('--dz',      type=float, default=0.2,  help='z step [A]')
args = ap.parse_args()

nz   = int(round((args.z_start - args.z_end) / args.dz)) + 1
dtip = -args.dz
tag  = f"{'morse' if args.morse else 'lj'}_{'qeq' if args.qeq else 'fixedQ'}"
print(f"=== test_ptcda.py  tag={tag}  nxy={args.nxy}  nz={nz}  dz={args.dz}  "
      f"probe_height {args.z_start:.1f}->{args.z_end:.1f} A ===")

# ── paths ─────────────────────────────────────────────────────────────────────
xyz_path    = os.path.join(_ROOT, 'cpp', 'common_resources', 'xyz', 'PTCDA.xyz')
params_path = os.path.join(_ROOT, 'cpp', 'common_resources', 'ElementTypes.dat')
assert os.path.exists(xyz_path),    f"PTCDA.xyz not found: {xyz_path}"
assert os.path.exists(params_path), f"ElementTypes.dat not found: {params_path}"

# ── build AFMulator ──────────────────────────────────────────────────────────
afm = AFMulator(use_morse=args.morse)
mol = afm.load_molecule(xyz_path)
print(f"Loaded {mol.natoms} atoms  enames set={set(mol.enames)}")
if mol.qs is not None:
    print(f"  charges from .xyz: min={mol.qs.min():.3f} max={mol.qs.max():.3f} sum={mol.qs.sum():.3f}")

afm.assign_params(params_path=params_path)
if args.qeq:
    qs = afm.solve_QEq(Q_total=0.0)
    print(f"  QEq q range=[{qs.min():.3f},{qs.max():.3f}] sum={qs.sum():.4f}")

L = afm.setup_grid(n=(80, 80, 60), margin=3.0, z_top=12.0)
print(f"Grid L={L}  n={afm.n}  mol_z={afm.mol_z:.2f} A (kernel space)")
afm.make_forcefield()

# ── Scan positions ───────────────────────────────────────────────────────────
mol_z  = afm.mol_z
nx_s, ny_s = args.nxy
apos = afm.atoms_arr[:,:3]
mn, mx = apos.min(axis=0), apos.max(axis=0)
x0 = mn[0] + (mx[0]-mn[0])*0.05
y0 = mn[1] + (mx[1]-mn[1])*0.05
z0_tip = mol_z + args.z_start + abs(float(afm.dpos0[2]))
scan_p0 = np.array([x0, y0, z0_tip], dtype=np.float32)
scan_da = np.array([(mx[0]-mn[0])*0.9/max(nx_s-1,1), 0., 0.], dtype=np.float32)
scan_db = np.array([0., (mx[1]-mn[1])*0.9/max(ny_s-1,1), 0.], dtype=np.float32)
probe_heights = mol_z + args.z_start + np.arange(nz) * dtip - mol_z
print(f"scan_p0={scan_p0}  probe_heights[0]={probe_heights[0]:.2f}  probe_heights[-1]={probe_heights[-1]:.2f} A")
x_ext = [float(scan_p0[0] - afm.mol_shift[0]),
          float(scan_p0[0] + scan_da[0]*(nx_s-1) - afm.mol_shift[0])]
y_ext = [float(scan_p0[1] - afm.mol_shift[1]),
          float(scan_p0[1] + scan_db[1]*(ny_s-1) - afm.mol_shift[1])]

# ── Raw + Relaxed scans ─────────────────────────────────────────────────────
print("\n--- RAW SCAN ---")
FEs_raw, _ = afm.get_raw_FE(nxy=tuple(args.nxy), nz=nz, dtip=dtip,
                              scan_p0=scan_p0, scan_da=scan_da, scan_db=scan_db)
Fz_raw = FEs_raw[:,:,:,2]
print(f"Fz_raw  stats: min={Fz_raw.min():.4f}  max={Fz_raw.max():.4f}  mean={Fz_raw.mean():.4f}")

print("\n--- RELAXED SCAN ---")
FEs_relax, _ = afm.run_scan(nxy=tuple(args.nxy), nz=nz, dtip=dtip,
                              scan_p0=scan_p0, scan_da=scan_da, scan_db=scan_db)
Fz_relax = FEs_relax[:,:,:,2]
print(f"Fz_relax stats: min={Fz_relax.min():.4f}  max={Fz_relax.max():.4f}  mean={Fz_relax.mean():.4f}")

# ── df ───────────────────────────────────────────────────────────────────────
df_raw   = compute_df(Fz_raw,   args.dz)
df_relax = compute_df(Fz_relax, args.dz)
print(f"df_relax stats: min={df_relax.min():.5f}  max={df_relax.max():.5f}")

# ── Save ─────────────────────────────────────────────────────────────────────
for name, arr in [('Fz_raw',Fz_raw), ('Fz_relax',Fz_relax), ('FEs_raw',FEs_raw),
                  ('FEs_relax',FEs_relax), ('df_relax',df_relax), ('df_raw',df_raw),
                  ('probe_heights',probe_heights)]:
    np.save(os.path.join(_THIS_DIR, f'{name}_{tag}.npy'), arr)
    print(f"Saved {name}_{tag}.npy  shape={arr.shape}")

# ── Plots ────────────────────────────────────────────────────────────────────
if args.noplot:
    print("=== DONE (--noplot) ==="); sys.exit(0)

import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

plot_grid_Fz(Fz_raw,   probe_heights, f'PTCDA Fz raw [{tag}]',   f'afm_Fz_raw_{tag}.png',   x_ext=x_ext, y_ext=y_ext, save_dir=_THIS_DIR)
plot_grid_Fz(Fz_relax, probe_heights, f'PTCDA Fz relax [{tag}]', f'afm_Fz_relax_{tag}.png', x_ext=x_ext, y_ext=y_ext, save_dir=_THIS_DIR)
plot_grid_Fz(df_relax, probe_heights, f'PTCDA df relax [{tag}]',  f'afm_df_relax_{tag}.png', x_ext=x_ext, y_ext=y_ext, save_dir=_THIS_DIR)

# Comparison: raw vs relax at representative heights
sel_heights = [5.0, 4.0, 3.5, 3.0, 2.5, 2.0]
sel_iz = [int(round((args.z_start - h) / args.dz)) for h in sel_heights]
sel_iz = [iz for iz in sel_iz if 0 <= iz < nz]
nsel = len(sel_iz)
kw_base = dict(origin='lower', cmap='bwr', aspect='auto',
               extent=[x_ext[0], x_ext[1], y_ext[0], y_ext[1]])
fig, axes = plt.subplots(2, nsel, figsize=(3*nsel, 6))
fig.suptitle(f"PTCDA raw vs relax [{tag}]", fontsize=11)
for col, iz in enumerate(sel_iz):
    h = probe_heights[iz]
    norm_raw   = safe_norm(Fz_raw[:,:,iz])
    norm_relax = safe_norm(Fz_relax[:,:,iz])
    im0 = axes[0,col].imshow(Fz_raw[:,:,iz].T,   norm=norm_raw,   **kw_base)
    im1 = axes[1,col].imshow(Fz_relax[:,:,iz].T, norm=norm_relax, **kw_base)
    axes[0,col].set_title(f"raw h={h:.1f}A +/-{norm_raw.vmax:.2g}", fontsize=7)
    axes[1,col].set_title(f"relax h={h:.1f}A +/-{norm_relax.vmax:.2g}", fontsize=7)
    plt.colorbar(im0, ax=axes[0,col], shrink=0.8); plt.colorbar(im1, ax=axes[1,col], shrink=0.8)
    for ax in [axes[0,col], axes[1,col]]: ax.tick_params(labelsize=5)
plt.tight_layout()
plt.savefig(os.path.join(_THIS_DIR, f'afm_comparison_{tag}.png'), dpi=110, bbox_inches='tight'); plt.close()
print(f"Saved afm_comparison_{tag}.png")

# Fz traces at atom positions
def nearest_scan_pixel(pos_xy, scan_p0, scan_da, scan_db, nx_s, ny_s):
    ix = int(round((pos_xy[0] - scan_p0[0]) / scan_da[0]))
    iy = int(round((pos_xy[1] - scan_p0[1]) / scan_db[1]))
    return max(0, min(nx_s-1, ix)), max(0, min(ny_s-1, iy))

atom_xy = afm.atoms_arr[:, :2]
centroid_xy = atom_xy.mean(axis=0)
dists = np.linalg.norm(atom_xy - centroid_xy, axis=1)
atom_pixels = [nearest_scan_pixel(atom_xy[ia], scan_p0, scan_da, scan_db, nx_s, ny_s)
               for ia in [int(np.argmin(dists)), int(np.argsort(dists)[len(dists)//2]), int(np.argmax(dists))]]
atom_labels = ['center-atom','mid-atom','edge-atom']

fig, axes = plt.subplots(2, 3, figsize=(14, 7))
fig.suptitle(f"PTCDA Fz traces [{tag}] (top=full, bottom=zoom +/-0.5)")
for col, ((ix,iy), lbl) in enumerate(zip(atom_pixels, atom_labels)):
    for row, (ax, ylim) in enumerate(zip(axes[:,col], [None, (-0.5, 0.5)])):
        ax.plot(probe_heights, Fz_raw[ix,iy,:],   color='steelblue', label='raw FF', lw=1.2)
        ax.plot(probe_heights, Fz_relax[ix,iy,:], color='tomato', lw=1.2, ls='--', label='PP relax')
        ax.axhline(0, color='k', lw=0.5); ax.axvline(3.2, color='gray', lw=0.8, ls=':', label='R0=3.2A')
        ax.set_xlabel('h (A)'); ax.set_ylabel('Fz (eV/A)'); ax.set_title(lbl if row==0 else f"{lbl} [zoom]", fontsize=8)
        ax.legend(fontsize=6); ax.invert_xaxis()
        if ylim: ax.set_ylim(ylim)
plt.tight_layout()
plt.savefig(os.path.join(_THIS_DIR, f'afm_traces_atoms_{tag}.png'), dpi=110, bbox_inches='tight'); plt.close()
print(f"Saved afm_traces_atoms_{tag}.png")

# LJ vs Morse comparison (if both npy files exist)
other_tag = tag.replace('morse','lj') if 'morse' in tag else tag.replace('lj','morse')
Fz_other_path = os.path.join(_THIS_DIR, f'Fz_relax_{other_tag}.npy')
hts_other_path = os.path.join(_THIS_DIR, f'probe_heights_{other_tag}.npy')
if os.path.exists(Fz_other_path) and os.path.exists(hts_other_path):
    Fz_other  = np.load(Fz_other_path)
    hts_other = np.load(hts_other_path)
    lbl_this  = 'Morse' if 'morse' in tag else 'LJ'
    lbl_other = 'LJ'    if 'morse' in tag else 'Morse'
    sel_iz_cmp = [int(round((args.z_start - h) / args.dz)) for h in [4.0, 3.5, 3.0, 2.5, 2.0, 1.6]]
    sel_iz_cmp = [iz for iz in sel_iz_cmp if 0 <= iz < nz]
    nsel_cmp = len(sel_iz_cmp)
    fig, axes = plt.subplots(2, nsel_cmp, figsize=(3*nsel_cmp, 6))
    fig.suptitle(f"PTCDA Fz: {lbl_this} vs {lbl_other}", fontsize=10)
    for col, iz in enumerate(sel_iz_cmp):
        h = probe_heights[iz]
        iz_other = min(int(np.argmin(np.abs(hts_other - h))), Fz_other.shape[2]-1)
        for row, (Fz_arr, lbl_row) in enumerate([(Fz_relax, lbl_this), (Fz_other, lbl_other)]):
            iz_use = iz if row==0 else iz_other
            norm = safe_norm(Fz_arr[:,:,iz_use])
            im = axes[row,col].imshow(Fz_arr[:,:,iz_use].T, norm=norm, **kw_base)
            axes[row,col].set_title(f"{lbl_row} h={h:.1f}A +/-{norm.vmax:.2g}", fontsize=7)
            plt.colorbar(im, ax=axes[row,col], shrink=0.8); axes[row,col].tick_params(labelsize=4)
    plt.tight_layout()
    plt.savefig(os.path.join(_THIS_DIR, f'afm_LJvsMorse_{tag}.png'), dpi=110, bbox_inches='tight'); plt.close()
    print(f"Saved afm_LJvsMorse_{tag}.png")

print("=== DONE ===")
