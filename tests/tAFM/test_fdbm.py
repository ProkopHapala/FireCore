#!/usr/bin/env python3
"""
test_fdbm.py  –  FDBM AFM of PTCDA using Fireball density.
Steps: 1)SCF 2)density projection 3)FFT Poisson 4)debug plots
       5)FDBM forces 6)PP relax 7)df 8)compare with Morse
"""
import sys, os, argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from scipy.ndimage import map_coordinates

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_ROOT     = os.path.realpath(os.path.join(_THIS_DIR, '..', '..'))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

XYZ_PATH    = os.path.join(_ROOT, 'cpp', 'common_resources', 'xyz', 'PTCDA.xyz')
_FDATA_HCNO = os.path.join(_ROOT, 'tests', 'Fireball', 'Fdata_HCNO')
_FDATA_BASIS = os.path.join(_FDATA_HCNO, 'basis')
COULOMB_CONST = 14.3996448915

ap = argparse.ArgumentParser()
ap.add_argument('--noplot',  action='store_true')
ap.add_argument('--nxy',     type=int, nargs=2, default=[60, 60])
ap.add_argument('--nscf',    type=int, default=200)
ap.add_argument('--dstep',   type=float, default=0.15)
ap.add_argument('--z_start', type=float, default=6.0)
ap.add_argument('--z_end',   type=float, default=0.4)
ap.add_argument('--dz',      type=float, default=0.2)
ap.add_argument('--A_pauli', type=float, default=3.0,    help='Pauli repulsion amplitude [eV·Å³]')
ap.add_argument('--C6_CO',   type=float, default=10.0,   help='CO-tip London C6 per atom [eV·Å^6]')
ap.add_argument('--q_CO',    type=float, default=-0.05,  help='CO-tip net charge [e]')
args = ap.parse_args()

nz        = int(round((args.z_start - args.z_end) / args.dz)) + 1
nx_s, ny_s = args.nxy
print(f"=== test_fdbm nxy={args.nxy} nz={nz} z={args.z_start}..{args.z_end} ===")

# ── STEP 1: Fireball SCF ──────────────────────────────────────────────────────
print("\n=== STEP 1: Fireball SCF ===")
_ORIG_CWD = os.getcwd()
_FBALL_CWD = os.path.join(_ROOT, 'tests', 'pyFireball')
os.chdir(_FBALL_CWD)
_fdata_local = os.path.join(_FBALL_CWD, 'Fdata')
if not os.path.exists(os.path.join(_fdata_local, 'info.dat')):
    if os.path.lexists(_fdata_local):
        os.unlink(_fdata_local)
    os.symlink(_FDATA_HCNO, _fdata_local)
    print(f"Created Fdata symlink → {_FDATA_HCNO}")

import pyBall.FireCore as fc

ELEM_Z = {'H':1,'C':6,'N':7,'O':8}
with open(XYZ_PATH) as f:
    lines = f.readlines()
natoms    = int(lines[0])
atomTypes = []
atomPos   = []
for line in lines[2:2+natoms]:
    p = line.split()
    atomTypes.append(ELEM_Z[p[0]])
    atomPos.append([float(p[1]), float(p[2]), float(p[3])])
atomTypes = np.array(atomTypes, dtype=np.int32)
atomPos   = np.array(atomPos,   dtype=np.float64)
print(f"PTCDA: {natoms} atoms  z=[{atomPos[:,2].min():.2f},{atomPos[:,2].max():.2f}]")

fc.setVerbosity(0)
fc.preinit()
fc.init(atomTypes, atomPos)
print(f"SCF (max {args.nscf} iter)...")
fc.SCF(atomPos, nmax_scf=args.nscf)

dims   = fc.get_HS_dims()
neighs = fc.get_HS_neighs(dims)
neighs = fc.get_rho_sparse(dims, data=neighs)
rho_sparse = neighs.rho
print(f"rho_sparse shape={rho_sparse.shape}  |max|={np.abs(rho_sparse).max():.4f}")

# firecore_getCharges returns dimension(natoms) populations (Mulliken or Lowdin)
# pass (natoms,1) so the C pointer points to contiguous natoms doubles
charges2d = np.zeros((natoms, 1), dtype=np.float64)
fc.getCharges(charges2d)
q_pop = charges2d[:,0]              # VALENCE electron populations per atom
# Fireball DFTB: only valence electrons. Z_valence: H=1, C=4, N=5, O=6, S=6
_Z_to_Zval = {1:1, 6:4, 7:5, 8:6, 16:6}
z_val = np.array([_Z_to_Zval.get(int(z), int(z)) for z in atomTypes], dtype=float)
q_mulliken = q_pop - z_val     # net Mulliken charge = valence_pop - Z_valence (≈0 for neutral)
print(f"Mulliken charges: sum={q_mulliken.sum():.4f}  range=[{q_mulliken.min():.3f},{q_mulliken.max():.3f}]")
os.chdir(_ORIG_CWD)

# ── STEP 2: Project density ───────────────────────────────────────────────────
print("\n=== STEP 2: Density projection ===")
from pyBall.FireballOCL import Grid as ocl_grid

_step   = args.dstep
_margin = 4.0
BLOCK   = 8
pos_min = atomPos.min(axis=0) - _margin
pos_max = atomPos.max(axis=0) + np.array([_margin, _margin, _margin + 6.0])
span    = pos_max - pos_min
ngrid   = (np.ceil(np.ceil(span/_step)/BLOCK).astype(int)*BLOCK)
total_span = ngrid * _step
origin  = (0.5*(pos_min+pos_max) - 0.5*total_span).astype(np.float32)
print(f"Grid: origin={origin.round(2)} step={_step} ngrid={ngrid} span={total_span.round(2)}")

grid_spec = {
    'origin': origin,
    'dA': [_step,0.,0.], 'dB': [0.,_step,0.], 'dC': [0.,0.,_step],
    'ngrid': ngrid.astype(int)
}
RCUT = {1:3.5, 6:5.0, 7:5.0, 8:5.0}
atoms_dict = {
    'pos':  atomPos,
    'Rcut': np.array([RCUT.get(int(z),4.5) for z in atomTypes]),
    'type': atomTypes,
}
projector = ocl_grid.GridProjector(fdata_dir=_FDATA_BASIS)
try:
    projector.load_basis(sorted(set(atomTypes.tolist())))
except Exception as e:
    print(f"[WARN] basis load: {e}")

print("Projecting density to grid...")
rho_grid = projector.project(rho_sparse, neighs, atoms_dict, grid_spec, nMaxAtom=64)
nx_d, ny_d, nz_d = rho_grid.shape
dV = _step**3
print(f"rho_grid shape={rho_grid.shape} range=[{rho_grid.min():.5f},{rho_grid.max():.5f}]")
print(f"  Integrated electrons = {rho_grid.sum()*dV:.2f}")
np.save(os.path.join(_THIS_DIR,'rho_grid_fdbm.npy'), rho_grid)

# ── STEP 3: FFT Poisson ────────────────────────────────────────────────────────
print("\n=== STEP 3: FFT Poisson solver ===")

def fft_poisson(rho, step):
    nx, ny, nz = rho.shape
    rho_k = np.fft.fftn(rho)
    kx = 2*np.pi * np.fft.fftfreq(nx, d=step)
    ky = 2*np.pi * np.fft.fftfreq(ny, d=step)
    kz = 2*np.pi * np.fft.fftfreq(nz, d=step)
    KX,KY,KZ = np.meshgrid(kx,ky,kz,indexing='ij')
    k2 = KX**2+KY**2+KZ**2; k2[0,0,0]=1.0
    V_k = 4.0*np.pi*COULOMB_CONST*rho_k/k2; V_k[0,0,0]=0.0
    return np.real(np.fft.ifftn(V_k)).astype(np.float32)

# Smooth rho for Pauli gradient (no Poisson needed - ES handled via Mulliken charges below)
from scipy.ndimage import gaussian_filter
sigma_elec = 0.8 / _step    # mild smoothing for gradient stability
rho_smooth = gaussian_filter(rho_grid.astype(np.float64), sigma=sigma_elec).astype(np.float32)
print(f"rho_smooth: max={rho_smooth.max():.5f}  integral={rho_smooth.sum()*dV:.2f} e")
# No V_ES from Poisson: Mulliken charges used in atom loop for electrostatics
V_H = np.zeros(rho_grid.shape, dtype=np.float32)   # placeholder for plot compat

# ── STEP 4: Debug plots ─────────────────────────────────────────────────────
def safe_norm(d2d, pct=99):
    v = max(float(np.percentile(np.abs(d2d),pct)),1e-8)
    return TwoSlopeNorm(vmin=-v,vcenter=0,vmax=v)

def plot_slices(data, title, fname, sym=False, cmap='magma'):
    nx,ny,nz = data.shape
    cx,cy,cz = nx//2, ny//2, nz//2
    if sym: cmap='bwr'
    fig, axes = plt.subplots(2,3,figsize=(16,8)); fig.suptitle(title)
    norm = safe_norm(data) if sym else None
    kw = dict(origin='lower',cmap=cmap,aspect='auto',norm=norm)
    for ax,sl,tl in zip(axes[0],
        [data[cx,:,:].T, data[:,cy,:].T, data[:,:,cz].T],
        [f'ix={cx} (YZ)', f'iy={cy} (XZ)', f'iz={cz} (XY)']):
        im=ax.imshow(sl,**kw); ax.set_title(tl); plt.colorbar(im,ax=ax,shrink=0.8)
    axes[1,0].plot(data[cx,cy,:]); axes[1,0].set_xlabel('iz'); axes[1,0].set_title('z-profile center')
    axes[1,1].plot(data[:,cy,cz]); axes[1,1].set_xlabel('ix'); axes[1,1].set_title('x-profile center')
    axes[1,2].plot(data[cx,:,cz]); axes[1,2].set_xlabel('iy'); axes[1,2].set_title('y-profile center')
    for ax in axes[1]: ax.axhline(0,color='k',lw=0.5)
    plt.tight_layout()
    plt.savefig(os.path.join(_THIS_DIR,fname),dpi=90,bbox_inches='tight'); plt.close()
    print(f"Saved {fname}")

if not args.noplot:
    plot_slices(rho_grid, 'PTCDA density ρ(r) [e/Å³] (Fireball STO)', 'fdbm_rho_slices.png')
    plot_slices(rho_smooth, 'PTCDA smoothed density ρ_smooth [e/Å³]', 'fdbm_VH_slices.png')

# ── STEP 5: FDBM force field ──────────────────────────────────────────────────
print("\n=== STEP 5: FDBM force field ===")

mol_z = 0.0
scan_x0 = atomPos[:,0].min() + (atomPos[:,0].max()-atomPos[:,0].min())*0.05
scan_y0 = atomPos[:,1].min() + (atomPos[:,1].max()-atomPos[:,1].min())*0.05
scan_dx = (atomPos[:,0].max()-atomPos[:,0].min())*0.9/max(nx_s-1,1)
scan_dy = (atomPos[:,1].max()-atomPos[:,1].min())*0.9/max(ny_s-1,1)
scan_xs = scan_x0 + np.arange(nx_s)*scan_dx
scan_ys = scan_y0 + np.arange(ny_s)*scan_dy
probe_heights = args.z_start + np.arange(nz)*(-args.dz)
np.save(os.path.join(_THIS_DIR,'probe_heights_fdbm.npy'), probe_heights)

XX,YY,ZZ = np.meshgrid(scan_xs, scan_ys, probe_heights+mol_z, indexing='ij')
flat_pos = np.stack([XX.ravel(),YY.ravel(),ZZ.ravel()],axis=1)   # (N,3)
grid_c   = (flat_pos - origin) / _step

def interp3(g, c): return map_coordinates(g, c.T, order=1, mode='nearest').astype(np.float32)

print("  Computing gradients...")
for axis, label in [(0,'x'),(1,'y'),(2,'z')]:
    pass  # compute below per-axis

grads_rho = [np.gradient(rho_smooth, _step, axis=a).astype(np.float32) for a in range(3)]
grads_VH  = [np.zeros_like(grads_rho[a]) for a in range(3)]   # unused (no V_H Poisson)

print("  Interpolating gradients at scan positions...")
F_fdbm = np.zeros((flat_pos.shape[0], 4), dtype=np.float32)
# F_Pauli = -A * grad(rho_smooth)  [repulsive density overlap]
for a in range(3):
    F_fdbm[:,a] += -args.A_pauli * interp3(grads_rho[a], grid_c)

print("  Atom loop: London vdW + Mulliken electrostatics...")
RA2_VDW = 1.5**2   # vdW softening
RA2_ES  = 2.0**2   # ES softening (avoids divergence, Mulliken at short range)
for ia in range(natoms):
    dr  = flat_pos - atomPos[ia]
    r2  = np.sum(dr**2,axis=1)
    r8_vdw = (r2 + RA2_VDW)**4
    r3_es  = (r2 + RA2_ES)**1.5
    for a in range(3):
        F_fdbm[:,a] -= args.C6_CO * 6.0 * dr[:,a] / r8_vdw          # London attractive
        F_fdbm[:,a] += COULOMB_CONST * q_mulliken[ia] * args.q_CO * dr[:,a] / r3_es   # Mulliken ES

FEs_raw = F_fdbm.reshape(nx_s, ny_s, nz, 4)
Fz_raw  = FEs_raw[:,:,:,2]
print(f"Fz_raw: min={Fz_raw.min():.4f}  max={Fz_raw.max():.4f}  mean={Fz_raw.mean():.4f} eV/Å")
np.save(os.path.join(_THIS_DIR,'FEs_raw_fdbm.npy'),FEs_raw)
np.save(os.path.join(_THIS_DIR,'Fz_raw_fdbm.npy'), Fz_raw)

# ── STEP 6: PP relaxation ─────────────────────────────────────────────────────
print("\n=== STEP 6: PP relaxation ===")
# Correct approach: compute forces in physical coordinates from density grid.
# For each height iz, build a 2D force map on (scan_region + lateral margin) at fixed z,
# then relax probes laterally using that map.
K_LAT   = 0.03    # eV/Å²  lateral spring stiffness
N_RELAX = 30      # relaxation iterations per height
PP_MARGIN = 2.0   # Å extra around scan region for relaxation

XX2, YY2 = np.meshgrid(scan_xs, scan_ys, indexing='ij')   # (nx_s, ny_s)

def fdbm_forces_at(positions):
    """Compute (Fx, Fy, Fz) at arbitrary positions (N,3) in Å. Returns (N,3) float32."""
    gc = (positions - origin) / _step
    F = np.zeros((len(positions), 3), dtype=np.float64)
    for a in range(3):
        F[:,a] += -args.A_pauli * interp3(grads_rho[a], gc).astype(np.float64)
    for ia in range(natoms):
        dr  = positions - atomPos[ia]
        r2  = np.sum(dr**2, axis=1)
        r8_vdw = (r2 + RA2_VDW)**4
        r3_es  = (r2 + RA2_ES)**1.5
        for a in range(3):
            F[:,a] -= args.C6_CO * 6.0 * dr[:,a] / r8_vdw
            F[:,a] += COULOMB_CONST * q_mulliken[ia] * args.q_CO * dr[:,a] / r3_es
    return F.astype(np.float32)

# Pre-build 2D lateral force maps per height for fast PP relaxation
# Grid covers scan region + PP_MARGIN, same step as density grid
pp_x0 = scan_xs[0]  - PP_MARGIN;  pp_x1 = scan_xs[-1]  + PP_MARGIN
pp_y0 = scan_ys[0]  - PP_MARGIN;  pp_y1 = scan_ys[-1]  + PP_MARGIN
pp_nx = int(np.ceil((pp_x1-pp_x0)/_step)) + 1
pp_ny = int(np.ceil((pp_y1-pp_y0)/_step)) + 1
pp_xs = pp_x0 + np.arange(pp_nx)*_step
pp_ys = pp_y0 + np.arange(pp_ny)*_step
PP_X, PP_Y = np.meshgrid(pp_xs, pp_ys, indexing='ij')   # (pp_nx, pp_ny)

def interp_phys_2d(F2d, x0, step, px, py):
    """Interpolate (pp_nx,pp_ny) field at physical positions px,py using grid coords."""
    cx = np.clip((px - x0) / step, 0, F2d.shape[0]-1.001)
    cy = np.clip((py - pp_y0) / step, 0, F2d.shape[1]-1.001)
    return map_coordinates(F2d, [cx.ravel(), cy.ravel()],
                           order=1, mode='nearest').reshape(px.shape).astype(np.float32)

FEs_relax = np.zeros_like(FEs_raw)

for iz in range(nz):
    probe_z  = probe_heights[iz] + mol_z
    PP_Z     = np.full_like(PP_X, probe_z)
    pp_flat  = np.stack([PP_X.ravel(), PP_Y.ravel(), PP_Z.ravel()], axis=1)
    print(f"  iz={iz:2d} h={probe_heights[iz]:.1f} Å  probe_z={probe_z:.2f} Å  building {pp_nx}x{pp_ny} force map...")
    FF = fdbm_forces_at(pp_flat)
    FF_x = FF[:,0].reshape(pp_nx, pp_ny)
    FF_y = FF[:,1].reshape(pp_nx, pp_ny)
    FF_z = FF[:,2].reshape(pp_nx, pp_ny)
    # Probe starts at scan equilibrium positions
    probe_x = XX2.astype(np.float32).copy()
    probe_y = YY2.astype(np.float32).copy()
    vx = np.zeros_like(probe_x)
    vy = np.zeros_like(probe_y)
    for _ in range(N_RELAX):
        Fx_s = interp_phys_2d(FF_x, pp_x0, _step, probe_x, probe_y) - K_LAT*(probe_x - XX2)
        Fy_s = interp_phys_2d(FF_y, pp_x0, _step, probe_x, probe_y) - K_LAT*(probe_y - YY2)
        vx = 0.8*vx + 0.3*Fx_s;  probe_x += vx*0.3
        vy = 0.8*vy + 0.3*Fy_s;  probe_y += vy*0.3
    # Record Fz (and lateral forces) at relaxed positions
    FEs_relax[:,:,iz,0] = interp_phys_2d(FF_x, pp_x0, _step, probe_x, probe_y)
    FEs_relax[:,:,iz,1] = interp_phys_2d(FF_y, pp_x0, _step, probe_x, probe_y)
    FEs_relax[:,:,iz,2] = interp_phys_2d(FF_z, pp_x0, _step, probe_x, probe_y)

Fz_relax = FEs_relax[:,:,:,2]
print(f"Fz_relax: min={Fz_relax.min():.4f}  max={Fz_relax.max():.4f}  mean={Fz_relax.mean():.4f} eV/Å")
np.save(os.path.join(_THIS_DIR,'FEs_relax_fdbm.npy'),FEs_relax)
np.save(os.path.join(_THIS_DIR,'Fz_relax_fdbm.npy'), Fz_relax)

# df = -dFz/dz
df_raw   = -np.gradient(Fz_raw,   abs(args.dz), axis=2)
df_relax = -np.gradient(Fz_relax, abs(args.dz), axis=2)
np.save(os.path.join(_THIS_DIR,'df_raw_fdbm.npy'),   df_raw)
np.save(os.path.join(_THIS_DIR,'df_relax_fdbm.npy'), df_relax)
print(f"df_relax: min={df_relax.min():.4f}  max={df_relax.max():.4f}")

# ── STEP 7: Plots ─────────────────────────────────────────────────────────────
if args.noplot:
    print("=== DONE (--noplot) ==="); sys.exit(0)

print("\n=== STEP 7: Plots ===")
x_ext = [float(scan_xs[0]), float(scan_xs[-1])]
y_ext = [float(scan_ys[0]), float(scan_ys[-1])]

def plot_grid_Fz(Fz, heights, label, fname, ncols=7):
    nz_p = len(heights)
    nrows = int(np.ceil(nz_p/ncols))
    fig,axes = plt.subplots(nrows,ncols,figsize=(2.5*ncols,2.8*nrows))
    axes = np.array(axes).reshape(nrows,ncols)
    fig.suptitle(f"PTCDA FDBM {label} (eV/Å) [per-slice]",fontsize=10)
    kw = dict(origin='lower',cmap='bwr',aspect='auto',
              extent=[x_ext[0],x_ext[1],y_ext[0],y_ext[1]])
    for k in range(nz_p):
        r,c = divmod(k,ncols); ax=axes[r,c]
        vabs = max(float(np.percentile(np.abs(Fz[:,:,k]),99)),1e-6)
        norm = TwoSlopeNorm(vmin=-vabs,vcenter=0,vmax=vabs)
        im = ax.imshow(Fz[:,:,k].T, norm=norm, **kw)
        ax.set_title(f"h={heights[k]:.1f}Å ±{vabs:.2g}",fontsize=6); ax.tick_params(labelsize=4)
        plt.colorbar(im,ax=ax,shrink=0.8)
    for k in range(nz_p,nrows*ncols):
        r,c=divmod(k,ncols); axes[r,c].set_visible(False)
    plt.tight_layout()
    plt.savefig(os.path.join(_THIS_DIR,fname),dpi=90,bbox_inches='tight'); plt.close()
    print(f"Saved {fname}")

plot_grid_Fz(Fz_raw,   probe_heights,'Fz BEFORE PP relax','fdbm_Fz_raw.png')
plot_grid_Fz(Fz_relax, probe_heights,'Fz AFTER PP relax', 'fdbm_Fz_relax.png')
plot_grid_Fz(df_relax, probe_heights,'df=-dFz/dz RELAX',  'fdbm_df_relax.png')

# Comparison raw vs relax at representative heights
sel_heights = [5.0,4.0,3.5,3.0,2.5,2.0]
sel_iz = [int(round((args.z_start-h)/args.dz)) for h in sel_heights]
sel_iz = [iz for iz in sel_iz if 0<=iz<nz]
nsel = len(sel_iz)
fig,axes = plt.subplots(2,nsel,figsize=(3*nsel,6))
fig.suptitle("PTCDA FDBM: raw vs PP-relaxed Fz",fontsize=11)
kw_b = dict(origin='lower',cmap='bwr',aspect='auto',
            extent=[x_ext[0],x_ext[1],y_ext[0],y_ext[1]])
for col,iz in enumerate(sel_iz):
    h=probe_heights[iz]
    for row,Fz_arr,lbl in [(0,Fz_raw,'raw'),(1,Fz_relax,'relax')]:
        vabs=max(float(np.percentile(np.abs(Fz_arr[:,:,iz]),99)),1e-6)
        norm=TwoSlopeNorm(vmin=-vabs,vcenter=0,vmax=vabs)
        im=axes[row,col].imshow(Fz_arr[:,:,iz].T,norm=norm,**kw_b)
        axes[row,col].set_title(f"{lbl} h={h:.1f}Å ±{vabs:.2g}",fontsize=7)
        plt.colorbar(im,ax=axes[row,col],shrink=0.8)
        axes[row,col].tick_params(labelsize=5)
plt.tight_layout()
plt.savefig(os.path.join(_THIS_DIR,'fdbm_comparison.png'),dpi=110,bbox_inches='tight'); plt.close()
print("Saved fdbm_comparison.png")

# FDBM vs Morse comparison (if Morse file present)
# Shows FDBM at its best heights (~1.6-2.4 Å) vs Morse at its best heights (~3.0-3.8 Å)
# (FDBM height scale is shifted ~1.5 Å lower due to diffuse STO density)
morse_path = os.path.join(_THIS_DIR,'Fz_relax_morse_fixedQ.npy')
morse_ht_path = os.path.join(_THIS_DIR,'probe_heights_morse_fixedQ.npy')
if os.path.exists(morse_path) and os.path.exists(morse_ht_path):
    hts_morse = np.load(morse_ht_path)
    df_morse_path = os.path.join(_THIS_DIR,'df_relax_morse_fixedQ.npy')
    df_morse = np.load(df_morse_path) if os.path.exists(df_morse_path) else None
    # FDBM ring features at h=1.6-2.4 Å; Morse at h=3.0-3.8 Å (STO density shift ~1.5 Å)
    fdbm_heights = [2.4, 2.2, 2.0, 1.8, 1.6]
    morse_heights = [3.8, 3.6, 3.4, 3.2, 3.0]
    ncols = min(len(fdbm_heights), len(morse_heights))
    fig, axes = plt.subplots(2, ncols, figsize=(3*ncols, 6))
    fig.suptitle("FDBM df (h=1.6-2.4Å) vs Morse df (h=3.0-3.8Å)\n(heights matched by ring-feature onset; STO density shifts FDBM ~1.5Å lower)",fontsize=8)
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
            axes[row,col].set_title(f"{lbl} h={h_use:.1f}Å ±{vabs:.3g}", fontsize=7)
            plt.colorbar(im, ax=axes[row,col], shrink=0.8)
            axes[row,col].tick_params(labelsize=4)
    plt.tight_layout()
    plt.savefig(os.path.join(_THIS_DIR,'fdbm_vs_morse.png'),dpi=110,bbox_inches='tight'); plt.close()
    print("Saved fdbm_vs_morse.png")

# Fz approach traces at representative pixels
centroid = atomPos[:,:2].mean(axis=0)
dists = np.linalg.norm(atomPos[:,:2]-centroid,axis=1)
ia_c = int(np.argmin(dists)); ia_e = int(np.argmax(dists))
ia_m = int(np.argsort(dists)[natoms//2])
def to_pixel(axy):
    ix=int(round((axy[0]-scan_x0)/scan_dx)); iy=int(round((axy[1]-scan_y0)/scan_dy))
    return max(0,min(nx_s-1,ix)), max(0,min(ny_s-1,iy))
trace_pixels = [to_pixel(atomPos[ia,:2]) for ia in [ia_c,ia_m,ia_e]]
trace_labels = ['center-atom','mid-atom','edge-atom']
fig,axes = plt.subplots(2,3,figsize=(14,7))
fig.suptitle("PTCDA FDBM Fz traces (top=full, bottom=zoomed ±0.5)")
for col,((ix,iy),lbl) in enumerate(zip(trace_pixels,trace_labels)):
    for row,(ax,ylim) in enumerate(zip(axes[:,col],[None,(-0.5,0.5)])):
        ax.plot(probe_heights,Fz_raw[ix,iy,:],   color='steelblue',label='raw FF',lw=1.2)
        ax.plot(probe_heights,Fz_relax[ix,iy,:], color='tomato',label='PP relax',lw=1.2,ls='--')
        ax.axhline(0,color='k',lw=0.5); ax.axvline(3.2,color='gray',lw=0.8,ls=':',label='R0=3.2Å')
        ax.set_xlabel('Probe height (Å)',fontsize=8); ax.set_ylabel('Fz (eV/Å)',fontsize=8)
        ax.set_title(lbl if row==0 else f"{lbl} [zoom]",fontsize=8)
        ax.legend(fontsize=6); ax.invert_xaxis()
        if ylim: ax.set_ylim(ylim); ax.axhspan(ylim[0],0,alpha=0.04,color='blue')
plt.tight_layout()
plt.savefig(os.path.join(_THIS_DIR,'fdbm_traces.png'),dpi=110,bbox_inches='tight'); plt.close()
print("Saved fdbm_traces.png")

print("\n=== DONE ===")
