#!/usr/bin/env python3
"""
test_fdbm.py  –  FDBM AFM of PTCDA using Fireball density.
Thin top-level script calling reusable functions from pyBall.OCL.AFM.
Steps: 1)SCF 2)density projection 3)FFT Poisson 4)debug plots
       5)FDBM forces 6)PP relax 7)df 8)compare with Morse
"""
import sys, os, argparse
import numpy as np

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_ROOT     = os.path.realpath(os.path.join(_THIS_DIR, '..', '..'))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from pyBall.OCL.AFM import (
    run_fireball_scf, setup_density_grid, build_neutral_atom_rho,
    project_density_grids, fft_poisson, compute_fdbm_forcefield,
    make_fdbm_force_func, pp_relax_2d, compute_df,
    plot_slices, plot_grid_Fz, safe_norm,
)

XYZ_PATH     = os.path.join(_ROOT, 'cpp', 'common_resources', 'xyz', 'PTCDA.xyz')
_FDATA_HCNO  = os.path.join(_ROOT, 'tests', 'Fireball', 'Fdata_HCNO')
_FDATA_BASIS = os.path.join(_FDATA_HCNO, 'basis')

# ── CLI ───────────────────────────────────────────────────────────────────────
ap = argparse.ArgumentParser()
ap.add_argument('--noplot',  action='store_true')
ap.add_argument('--nxy',     type=int, nargs=2, default=[60, 60])
ap.add_argument('--nscf',    type=int, default=200)
ap.add_argument('--dstep',   type=float, default=0.15)
ap.add_argument('--z_start', type=float, default=6.0)
ap.add_argument('--z_end',   type=float, default=0.4)
ap.add_argument('--dz',      type=float, default=0.2)
ap.add_argument('--A_pauli', type=float, default=16.0,   help='Pauli repulsion amplitude [eV*A^3/e^2]')
ap.add_argument('--C6_CO',   type=float, default=30.0,   help='CO-tip London C6 per atom [eV*A^6]')
ap.add_argument('--q_CO',    type=float, default=-0.05,  help='CO-tip net charge [e]')
ap.add_argument('--sigma_tip', type=float, default=0.7,  help='Gaussian tip width sigma [A]')
ap.add_argument('--es_model', choices=['mulliken','poisson'], default='poisson')
ap.add_argument('--force_model', choices=['gradient','fdbm'], default='fdbm')
ap.add_argument('--proj_notiled', action='store_true')
ap.add_argument('--proj_debug_exit', action='store_true')
ap.add_argument('--proj_debug_clear', action='store_true')
ap.add_argument('--proj_debug_return0', action='store_true')
ap.add_argument('--proj_debug_read_task', action='store_true')
ap.add_argument('--proj_debug_read_grid', action='store_true')
ap.add_argument('--stop_after_proj', action='store_true')
ap.add_argument('--verbose', type=int, default=0)
args = ap.parse_args()

nz         = int(round((args.z_start - args.z_end) / args.dz)) + 1
nx_s, ny_s = args.nxy
_step      = args.dstep
print(f"=== test_fdbm nxy={args.nxy} nz={nz} z={args.z_start}..{args.z_end} ===")

# ── STEP 1: Fireball SCF ─────────────────────────────────────────────────────
print("\n=== STEP 1: Fireball SCF ===")
scf = run_fireball_scf(XYZ_PATH, _FDATA_HCNO, nscf=args.nscf, verbosity=args.verbose)
atomTypes   = scf['atomTypes']
atomPos     = scf['atomPos']
rho_sparse  = scf['rho_sparse']
neighs      = scf['neighs']
q_mulliken  = scf['q_mulliken']
natoms      = scf['natoms']

# ── STEP 2: Density projection ───────────────────────────────────────────────
print("\n=== STEP 2: Density projection ===")
grid_spec, origin, ngrid = setup_density_grid(atomPos, step=_step)
rho_na_sparse = build_neutral_atom_rho(atomTypes, neighs, natoms)
proj_kwargs = dict(
    use_tiled=(not args.proj_notiled),
    debug_early_exit=args.proj_debug_exit,
    debug_clear_only=args.proj_debug_clear,
    debug_return0=args.proj_debug_return0,
    debug_read_task=args.proj_debug_read_task,
    debug_read_grid=args.proj_debug_read_grid,
)
rho_grid, rho_na_grid, rho_diff, projector = project_density_grids(
    rho_sparse, rho_na_sparse, neighs, atomTypes, atomPos, grid_spec,
    _FDATA_BASIS, step=_step, proj_kwargs=proj_kwargs, verbosity=args.verbose)
np.save(os.path.join(_THIS_DIR, 'rho_grid_fdbm.npy'), rho_grid)
np.save(os.path.join(_THIS_DIR, 'rho_na_grid_fdbm.npy'), rho_na_grid)
np.save(os.path.join(_THIS_DIR, 'rho_diff_grid_fdbm.npy'), rho_diff)

if args.stop_after_proj:
    print("=== STOP_AFTER_PROJ ==="); sys.exit(0)

# ── STEP 3: FFT Poisson ─────────────────────────────────────────────────────
print("\n=== STEP 3: FFT Poisson solver ===")
from scipy.ndimage import gaussian_filter
sigma_elec = 0.8 / _step
rho_smooth = gaussian_filter(rho_grid.astype(np.float64), sigma=sigma_elec).astype(np.float32)
dV = _step**3
print(f"rho_smooth: max={rho_smooth.max():.5f}  integral={rho_smooth.sum()*dV:.2f} e")
V_ES = fft_poisson(-rho_diff.astype(np.float64), _step)
print(f"V_ES(delta): range=[{V_ES.min():.3f},{V_ES.max():.3f}] eV")

# ── STEP 4: Debug plots ─────────────────────────────────────────────────────
if not args.noplot:
    plot_slices(rho_grid,    'PTCDA rho(r) [e/A^3] (Fireball STO)', 'fdbm_rho_slices.png', save_dir=_THIS_DIR)
    plot_slices(rho_smooth,  'PTCDA smoothed rho [e/A^3]',           'fdbm_VH_slices.png', save_dir=_THIS_DIR)
    plot_slices(rho_na_grid, 'Neutral-atom rho_NA(r) [e/A^3]',      'fdbm_rhoNA_slices.png', save_dir=_THIS_DIR)
    plot_slices(rho_diff,    'Delta_rho = rho_SCF - rho_NA [e/A^3]', 'fdbm_rhoDiff_slices.png', sym=True, cmap='bwr', save_dir=_THIS_DIR)
    plot_slices(V_ES,        'ESP from delta_rho (Poisson) [eV]',    'fdbm_VES_slices.png', sym=True, cmap='bwr', save_dir=_THIS_DIR)

# ── STEP 5: FDBM force field ────────────────────────────────────────────────
print("\n=== STEP 5: FDBM force field ===")
mol_z = 0.0
scan_x0 = atomPos[:,0].min() + (atomPos[:,0].max()-atomPos[:,0].min())*0.05
scan_y0 = atomPos[:,1].min() + (atomPos[:,1].max()-atomPos[:,1].min())*0.05
scan_dx = (atomPos[:,0].max()-atomPos[:,0].min())*0.9 / max(nx_s-1, 1)
scan_dy = (atomPos[:,1].max()-atomPos[:,1].min())*0.9 / max(ny_s-1, 1)
scan_xs = scan_x0 + np.arange(nx_s)*scan_dx
scan_ys = scan_y0 + np.arange(ny_s)*scan_dy
probe_heights = args.z_start + np.arange(nz)*(-args.dz)
np.save(os.path.join(_THIS_DIR, 'probe_heights_fdbm.npy'), probe_heights)

ff = compute_fdbm_forcefield(
    rho_grid, V_ES, origin, _step, atomPos,
    scan_xs, scan_ys, probe_heights, mol_z=mol_z,
    A_pauli=args.A_pauli, C6_CO=args.C6_CO, q_CO=args.q_CO, sigma_tip=args.sigma_tip,
    force_model=args.force_model, es_model=args.es_model,
    q_mulliken=q_mulliken, rho_smooth=rho_smooth)
FEs_raw = ff['FEs_raw'];  Fz_raw = ff['Fz_raw']
np.save(os.path.join(_THIS_DIR, 'FEs_raw_fdbm.npy'), FEs_raw)
np.save(os.path.join(_THIS_DIR, 'Fz_raw_fdbm.npy'), Fz_raw)

# Energy field diagnostic plots
if not args.noplot and args.force_model == 'fdbm':
    plot_slices(ff['E_Pauli_field'], f'E_Pauli (A={args.A_pauli}) [eV]', 'fdbm_E_Pauli_slices.png', save_dir=_THIS_DIR)
    plot_slices(ff['E_ES_field'],    f'E_ES (q={args.q_CO}) [eV]', 'fdbm_E_ES_slices.png', sym=True, cmap='bwr', save_dir=_THIS_DIR)

# ── STEP 6: PP relaxation ───────────────────────────────────────────────────
print("\n=== STEP 6: PP relaxation ===")
grads_pauli = ff.get('grads_E_Pauli', None)
grads_es    = ff.get('grads_E_ES', None)
force_func = make_fdbm_force_func(
    origin, _step, grads_pauli, grads_es, atomPos,
    C6_CO=args.C6_CO, force_model=args.force_model,
    A_pauli=args.A_pauli, q_CO=args.q_CO,
    grads_rho=ff.get('grads_rho'), grads_VES=ff.get('grads_VES'),
    es_model=args.es_model, q_mulliken=q_mulliken)
FEs_relax = pp_relax_2d(force_func, scan_xs, scan_ys, probe_heights,
                         mol_z=mol_z, step=_step)
Fz_relax = FEs_relax[:,:,:,2]
np.save(os.path.join(_THIS_DIR, 'FEs_relax_fdbm.npy'), FEs_relax)
np.save(os.path.join(_THIS_DIR, 'Fz_relax_fdbm.npy'), Fz_relax)

# df = -dFz/dz
df_raw   = compute_df(Fz_raw,   args.dz)
df_relax = compute_df(Fz_relax, args.dz)
np.save(os.path.join(_THIS_DIR, 'df_raw_fdbm.npy'),   df_raw)
np.save(os.path.join(_THIS_DIR, 'df_relax_fdbm.npy'), df_relax)
print(f"df_relax: min={df_relax.min():.4f}  max={df_relax.max():.4f}")

# ── STEP 7: Plots ───────────────────────────────────────────────────────────
if args.noplot:
    print("=== DONE (--noplot) ==="); sys.exit(0)

print("\n=== STEP 7: Plots ===")
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

x_ext = [float(scan_xs[0]), float(scan_xs[-1])]
y_ext = [float(scan_ys[0]), float(scan_ys[-1])]

plot_grid_Fz(Fz_raw,   probe_heights, 'Fz BEFORE PP relax', 'fdbm_Fz_raw.png',   x_ext=x_ext, y_ext=y_ext, save_dir=_THIS_DIR)
plot_grid_Fz(Fz_relax, probe_heights, 'Fz AFTER PP relax',  'fdbm_Fz_relax.png',  x_ext=x_ext, y_ext=y_ext, save_dir=_THIS_DIR)
plot_grid_Fz(df_relax, probe_heights, 'df=-dFz/dz RELAX',   'fdbm_df_relax.png',  x_ext=x_ext, y_ext=y_ext, save_dir=_THIS_DIR)

# Component traces at representative atoms
Fz_pauli = ff['F_pauli'][:,:,:,2]; Fz_vdw = ff['F_vdw'][:,:,:,2]; Fz_es = ff['F_es'][:,:,:,2]
centroid = atomPos[:,:2].mean(axis=0)
dists = np.linalg.norm(atomPos[:,:2]-centroid, axis=1)
ia_c = int(np.argmin(dists)); ia_e = int(np.argmax(dists)); ia_m = int(np.argsort(dists)[natoms//2])
def to_pixel(axy):
    ix = int(round((axy[0]-scan_x0)/scan_dx)); iy = int(round((axy[1]-scan_y0)/scan_dy))
    return max(0,min(nx_s-1,ix)), max(0,min(ny_s-1,iy))
trace_pixels = [to_pixel(atomPos[ia,:2]) for ia in [ia_c,ia_m,ia_e]]
trace_labels = ['center-atom','mid-atom','edge-atom']
kw_b = dict(origin='lower', cmap='bwr', aspect='auto', extent=[x_ext[0],x_ext[1],y_ext[0],y_ext[1]])
fig, axes = plt.subplots(2, 3, figsize=(14, 7))
fig.suptitle(f'FDBM Fz components [{args.force_model}]', fontsize=11)
for col, ((ix,iy), lbl) in enumerate(zip(trace_pixels, trace_labels)):
    for row,(ax,ylim) in enumerate(zip(axes[:,col],[None,(-0.5,0.5)])):
        ax.plot(probe_heights, Fz_pauli[ix,iy,:], color='purple',    lw=1.0, label='Pauli')
        ax.plot(probe_heights, Fz_vdw[ix,iy,:],   color='royalblue', lw=1.0, label='vdW')
        ax.plot(probe_heights, Fz_es[ix,iy,:],     color='seagreen',  lw=1.0, label=f"ES[{args.es_model}]")
        ax.plot(probe_heights, Fz_raw[ix,iy,:],    color='k',         lw=1.0, ls='--', label='sum')
        ax.axhline(0,color='k',lw=0.5); ax.set_xlabel('h (A)'); ax.set_ylabel('Fz (eV/A)')
        ax.set_title(lbl if row==0 else f"{lbl} [zoom]",fontsize=8); ax.legend(fontsize=6); ax.invert_xaxis()
        if ylim: ax.set_ylim(ylim)
plt.tight_layout()
plt.savefig(os.path.join(_THIS_DIR,'fdbm_components_traces_raw.png'), dpi=120, bbox_inches='tight'); plt.close()
print('Saved fdbm_components_traces_raw.png')

# FDBM vs Morse comparison (if Morse npy files present)
morse_path = os.path.join(_THIS_DIR,'Fz_relax_morse_fixedQ.npy')
morse_ht_path = os.path.join(_THIS_DIR,'probe_heights_morse_fixedQ.npy')
if os.path.exists(morse_path) and os.path.exists(morse_ht_path):
    hts_morse = np.load(morse_ht_path)
    df_morse_path = os.path.join(_THIS_DIR,'df_relax_morse_fixedQ.npy')
    df_morse = np.load(df_morse_path) if os.path.exists(df_morse_path) else None
    fdbm_heights = [3.6, 3.4, 3.2, 3.0, 2.8]
    morse_heights = [3.8, 3.6, 3.4, 3.2, 3.0]
    ncols = min(len(fdbm_heights), len(morse_heights))
    fig, axes = plt.subplots(2, ncols, figsize=(3*ncols, 6))
    fig.suptitle("FDBM df vs Morse df (heights matched by ring-feature onset)", fontsize=8)
    for col in range(ncols):
        h_f = fdbm_heights[col];  iz_f = int(np.argmin(np.abs(probe_heights - h_f)))
        h_m = morse_heights[col]; iz_m = int(np.argmin(np.abs(hts_morse - h_m)))
        fdbm_data = df_relax if iz_f < df_relax.shape[2] else Fz_relax
        morse_data = df_morse if df_morse is not None else np.load(morse_path)
        for row, (arr, iz_use, lbl, h_use) in enumerate([
                (fdbm_data, iz_f, 'FDBM df', h_f),
                (morse_data, iz_m, 'Morse df', h_m)]):
            vabs = max(float(np.percentile(np.abs(arr[:,:,iz_use]), 99)), 1e-6)
            norm = TwoSlopeNorm(vmin=-vabs, vcenter=0, vmax=vabs)
            im = axes[row,col].imshow(arr[:,:,iz_use].T, norm=norm, **kw_b)
            axes[row,col].set_title(f"{lbl} h={h_use:.1f}A +/-{vabs:.3g}", fontsize=7)
            plt.colorbar(im, ax=axes[row,col], shrink=0.8); axes[row,col].tick_params(labelsize=4)
    plt.tight_layout()
    plt.savefig(os.path.join(_THIS_DIR,'fdbm_vs_morse.png'), dpi=110, bbox_inches='tight'); plt.close()
    print("Saved fdbm_vs_morse.png")

# Fz approach traces
fig,axes = plt.subplots(2, 3, figsize=(14, 7))
fig.suptitle("PTCDA FDBM Fz traces (top=full, bottom=zoom +/-0.5)")
for col,((ix,iy),lbl) in enumerate(zip(trace_pixels, trace_labels)):
    for row,(ax,ylim) in enumerate(zip(axes[:,col],[None,(-0.5,0.5)])):
        ax.plot(probe_heights, Fz_raw[ix,iy,:],   color='steelblue', label='raw FF', lw=1.2)
        ax.plot(probe_heights, Fz_relax[ix,iy,:], color='tomato', label='PP relax', lw=1.2, ls='--')
        ax.axhline(0,color='k',lw=0.5); ax.axvline(3.2,color='gray',lw=0.8,ls=':',label='R0=3.2A')
        ax.set_xlabel('h (A)'); ax.set_ylabel('Fz (eV/A)')
        ax.set_title(lbl if row==0 else f"{lbl} [zoom]", fontsize=8); ax.legend(fontsize=6); ax.invert_xaxis()
        if ylim: ax.set_ylim(ylim)
plt.tight_layout()
plt.savefig(os.path.join(_THIS_DIR,'fdbm_traces.png'), dpi=110, bbox_inches='tight'); plt.close()
print("Saved fdbm_traces.png")

print("\n=== DONE ===")
