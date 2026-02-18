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
ap.add_argument('--es_model', choices=['mulliken','poisson'], default='mulliken', help='Electrostatics model: mulliken point charges or Poisson from delta density')
ap.add_argument('--proj_notiled', action='store_true', help='Use non-tiled density projection kernel (debug GPU crashes)')
ap.add_argument('--proj_debug_exit', action='store_true', help='Build projection kernels with DEBUG_EARLY_EXIT=1 (debug GPU crashes)')
ap.add_argument('--proj_debug_clear', action='store_true', help='Build projection kernels with DEBUG_CLEAR_ONLY=1 (debug GPU crashes)')
ap.add_argument('--proj_debug_return0', action='store_true', help='Build projection kernels with DEBUG_RETURN0=1 (debug GPU crashes)')
ap.add_argument('--proj_debug_read_task', action='store_true', help='Build projection kernels with DEBUG_READ_TASK=1 (debug GPU crashes)')
ap.add_argument('--proj_debug_read_grid', action='store_true', help='Build projection kernels with DEBUG_READ_GRID=1 (debug GPU crashes)')
ap.add_argument('--proj_ntasks', type=int, default=0, help='Limit number of projection tasks (blocks) for debugging; 0=all')
ap.add_argument('--stop_after_proj', action='store_true', help='Exit after density projection step (for crash bisection)')
ap.add_argument('--verbose', type=int, default=0, help='Verbosity level for GridProjector (0=quiet, 1=debug, 2=timing)')
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

# DEBUG: quick neighbor list sanity (is there a self-neighbor entry?)
try:
    nj0 = neighs.neigh_j.reshape(natoms, -1)
    print(f"neigh_j[0,:10]={nj0[0,:10]}")
except Exception as e:
    print(f"[WARN] cannot print neigh_j head: {e}")

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
# Wf cutoffs in Angstrom (.wf rcutoff is in Bohr; abohr=0.529177) + 0.2 Å margin
RCUT = {1:2.3, 6:2.6, 7:2.6, 8:2.5}
atoms_dict = {
    'pos':  atomPos,
    'Rcut': np.array([RCUT.get(int(z),4.5) for z in atomTypes]),
    'type': atomTypes,
}
projector = ocl_grid.GridProjector(
    fdata_dir=_FDATA_BASIS,
    debug_early_exit=args.proj_debug_exit,
    debug_clear_only=args.proj_debug_clear,
    debug_return0=args.proj_debug_return0,
    debug_read_task=args.proj_debug_read_task,
    debug_read_grid=args.proj_debug_read_grid,
    verbosity=args.verbose,
)
try:
    projector.load_basis(sorted(set(atomTypes.tolist())))
except Exception as e:
    print(f"[WARN] basis load: {e}")

print("Projecting density to grid...")
proj_tasks = None
if args.proj_ntasks and args.proj_ntasks > 0:
    print(f"[DEBUG] Building tasks on host and slicing to first {args.proj_ntasks} tasks")
    tasks_np, task_atoms_np = projector.build_tasks(atoms_dict, grid_spec, nMaxAtom=64)
    nt = min(int(args.proj_ntasks), len(tasks_np))
    proj_tasks = (tasks_np[:nt].copy(), task_atoms_np[:nt].copy())
    print(f"[DEBUG] proj_tasks: nt={nt} na_max={int(tasks_np['na'][:nt].max()) if nt>0 else 0}")
rho_grid = projector.project(rho_sparse, neighs, atoms_dict, grid_spec, tasks=proj_tasks, nMaxAtom=64, use_tiled=(not args.proj_notiled))
nx_d, ny_d, nz_d = rho_grid.shape
dV = _step**3
print(f"rho_grid shape={rho_grid.shape} range=[{rho_grid.min():.5f},{rho_grid.max():.5f}]")
print(f"  Integrated electrons = {rho_grid.sum()*dV:.2f}")
np.save(os.path.join(_THIS_DIR,'rho_grid_fdbm.npy'), rho_grid)

# Build a simple neutral-atom (promolecule) density on the same sparse layout:
# Only on-site (i,i) blocks are set, off-diagonal are zero.
print("Building neutral-atom density matrix (on-site only)...")
rho_na = np.zeros_like(rho_sparse, dtype=np.float32)
neigh_j = neighs.neigh_j.reshape(natoms, -1)
def _onsite_occ(Z):
    # Returns diag occupations for (s,py,pz,px) in the kernel convention used in Grid.cl.
    # This is a crude promolecule guess (valence only).
    if Z == 1:  # H: 1s1
        return np.array([1.0, 0.0, 0.0, 0.0], dtype=np.float32)
    if Z == 6:  # C: 2s2 2p2
        return np.array([2.0, 2.0/3.0, 2.0/3.0, 2.0/3.0], dtype=np.float32)
    if Z == 7:  # N: 2s2 2p3
        return np.array([2.0, 1.0, 1.0, 1.0], dtype=np.float32)
    if Z == 8:  # O: 2s2 2p4
        return np.array([2.0, 4.0/3.0, 4.0/3.0, 4.0/3.0], dtype=np.float32)
    # fallback: put all valence into s
    return np.array([float(_Z_to_Zval.get(int(Z), int(Z))), 0.0, 0.0, 0.0], dtype=np.float32)

for i in range(natoms):
    # find self-neighbor slot where neigh_j == i+1
    slots = np.where(neigh_j[i] == (i+1))[0]
    if len(slots) == 0:
        raise RuntimeError(f"No self-neighbor slot found for atom i={i}")
    iself = int(slots[0])
    occ = _onsite_occ(int(atomTypes[i]))
    rho_na[i, iself, :, :] = 0.0
    rho_na[i, iself, 0, 0] = occ[0]
    rho_na[i, iself, 1, 1] = occ[1]
    rho_na[i, iself, 2, 2] = occ[2]
    rho_na[i, iself, 3, 3] = occ[3]

rho_na_grid = projector.project(rho_na, neighs, atoms_dict, grid_spec, tasks=proj_tasks, nMaxAtom=64, use_tiled=(not args.proj_notiled))
print(f"rho_na_grid: range=[{rho_na_grid.min():.5f},{rho_na_grid.max():.5f}] integral={rho_na_grid.sum()*dV:.2f} e")
rho_diff = (rho_grid - rho_na_grid).astype(np.float32)
print(f"rho_diff: range=[{rho_diff.min():.5f},{rho_diff.max():.5f}] integral={rho_diff.sum()*dV:.2f} e")
np.save(os.path.join(_THIS_DIR,'rho_na_grid_fdbm.npy'), rho_na_grid)
np.save(os.path.join(_THIS_DIR,'rho_diff_grid_fdbm.npy'), rho_diff)

if args.stop_after_proj:
    print("=== STOP_AFTER_PROJ ===")
    sys.exit(0)

def _clip3(i, n):
    return max(0, min(n-1, int(i)))

def plot_local_zooms(rho, title, fname, atom_indices, half_box_A=3.0):
    """Show XY zooms around selected atoms at atom z-plane + +/- a few slices."""
    half = int(round(half_box_A / _step))
    fig, axes = plt.subplots(len(atom_indices), 3, figsize=(12, 3.2*len(atom_indices)))
    if len(atom_indices) == 1:
        axes = np.array([axes])
    fig.suptitle(title, fontsize=11)
    for row, ia in enumerate(atom_indices):
        gc = (atomPos[ia] - origin) / _step
        ix, iy, iz = (_clip3(round(gc[0]), nx_d), _clip3(round(gc[1]), ny_d), _clip3(round(gc[2]), nz_d))
        slx = slice(max(0, ix-half), min(nx_d, ix+half+1))
        sly = slice(max(0, iy-half), min(ny_d, iy+half+1))
        # show iz-2, iz, iz+2 (clipped)
        izs = [_clip3(iz-2, nz_d), iz, _clip3(iz+2, nz_d)]
        for col, iz0 in enumerate(izs):
            data = rho[slx, sly, iz0].T
            vmax = float(np.percentile(data, 99.5))
            im = axes[row, col].imshow(data, origin='lower', cmap='magma', vmin=0.0, vmax=max(vmax, 1e-8), aspect='equal')
            axes[row, col].set_title(f"ia={ia} Z={atomTypes[ia]} iz={iz0}  max~{data.max():.3g}", fontsize=8)
            plt.colorbar(im, ax=axes[row, col], shrink=0.8)
            axes[row, col].tick_params(labelsize=6)
    plt.tight_layout()
    plt.savefig(os.path.join(_THIS_DIR, fname), dpi=120, bbox_inches='tight'); plt.close()
    print(f"Saved {fname}")

def radial_profile(rho, rmax_A=8.0, nbins=200, ia=0):
    """Radial average rho(r) around atom ia from local subgrid (for smoothness/decay diagnostics)."""
    gc = (atomPos[ia] - origin) / _step
    ix, iy, iz = (int(round(gc[0])), int(round(gc[1])), int(round(gc[2])))
    R = int(np.ceil(rmax_A / _step))
    x0, x1 = max(0, ix-R), min(nx_d, ix+R+1)
    y0, y1 = max(0, iy-R), min(ny_d, iy+R+1)
    z0, z1 = max(0, iz-R), min(nz_d, iz+R+1)
    xs = (np.arange(x0, x1) * _step + origin[0]).astype(np.float64)
    ys = (np.arange(y0, y1) * _step + origin[1]).astype(np.float64)
    zs = (np.arange(z0, z1) * _step + origin[2]).astype(np.float64)
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing='ij')
    d = np.sqrt((X-atomPos[ia,0])**2 + (Y-atomPos[ia,1])**2 + (Z-atomPos[ia,2])**2)
    v = rho[x0:x1, y0:y1, z0:z1].astype(np.float64)
    d = d.ravel(); v = v.ravel()
    m = d <= rmax_A
    d = d[m]; v = v[m]
    bins = np.linspace(0.0, rmax_A, nbins+1)
    ib = np.clip(np.searchsorted(bins, d, side='right')-1, 0, nbins-1)
    sum_v = np.bincount(ib, weights=v, minlength=nbins)
    cnt_v = np.bincount(ib, minlength=nbins)
    avg = sum_v / np.maximum(cnt_v, 1)
    rc = 0.5*(bins[:-1] + bins[1:])
    return rc, avg, cnt_v

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
# Poisson ES from delta-density (SCF - neutral atoms)
V_ES = fft_poisson(-rho_diff.astype(np.float64), _step).astype(np.float32)   # charge density = -rho_e
print(f"V_ES(delta): range=[{V_ES.min():.3f},{V_ES.max():.3f}] eV")
V_H = V_ES

# Density projection diagnostics (rings/steps around H etc.)
if not args.noplot:
    print("  Density projection diagnostics...")
    # pick one H (if present) and one C near center
    iHs = np.where(atomTypes == 1)[0]
    iCs = np.where(atomTypes == 6)[0]
    iaH = int(iHs[0]) if len(iHs) > 0 else 0
    centroid = atomPos[:,:2].mean(axis=0)
    if len(iCs) > 0:
        dCs = np.linalg.norm(atomPos[iCs,:2] - centroid[None,:], axis=1)
        iaC = int(iCs[int(np.argmin(dCs))])
    else:
        iaC = int(np.argmin(np.linalg.norm(atomPos[:,:2] - centroid[None,:], axis=1)))
    # Use a larger box to catch any hard Rcut artifacts (rims/shells)
    plot_local_zooms(rho_grid,   'ρ projection local zooms (raw)',    'fdbm_rho_zooms_raw.png',    [iaH, iaC], half_box_A=6.0)
    plot_local_zooms(rho_smooth, 'ρ projection local zooms (smooth)', 'fdbm_rho_zooms_smooth.png', [iaH, iaC], half_box_A=6.0)
    for ia, tag in [(iaH, 'H'), (iaC, 'C')]:
        rc, avg, cnt = radial_profile(rho_grid, rmax_A=8.0, nbins=220, ia=ia)
        fig, axes = plt.subplots(1, 2, figsize=(10, 3.2))
        fig.suptitle(f"Radial profile around atom ia={ia} Z={atomTypes[ia]} ({tag})", fontsize=10)
        axes[0].plot(rc, avg, lw=1.2)
        axes[0].set_xlabel('r (Å)'); axes[0].set_ylabel('⟨ρ⟩ (e/Å³)')
        axes[0].set_title('linear')
        axes[1].semilogy(rc, np.maximum(avg, 1e-12), lw=1.2)
        axes[1].set_xlabel('r (Å)'); axes[1].set_ylabel('⟨ρ⟩ (e/Å³)')
        axes[1].set_title('log')
        for ax in axes:
            ax.grid(True, alpha=0.2)
        plt.tight_layout()
        fn = f"fdbm_rho_radial_{tag}.png"
        plt.savefig(os.path.join(_THIS_DIR, fn), dpi=130, bbox_inches='tight'); plt.close()
        print(f"Saved {fn}")

    # 1D z-lineouts above the chosen atoms (this avoids contamination by other atoms)
    def lineout_z(rho, ia, zmax_A=10.0):
        gc = (atomPos[ia] - origin) / _step
        ix = _clip3(round(gc[0]), nx_d)
        iy = _clip3(round(gc[1]), ny_d)
        iz0 = _clip3(round(gc[2]), nz_d)
        npts = int(round(zmax_A / _step))
        izs = np.arange(iz0, min(nz_d, iz0+npts))
        zsA = (izs * _step + origin[2]) - atomPos[ia,2]
        return zsA.astype(np.float64), rho[ix, iy, izs].astype(np.float64)

    fig, axes = plt.subplots(1, 2, figsize=(11, 3.2))
    fig.suptitle('1D density lineout along +z above atoms (checks smooth decay / no steps)', fontsize=10)
    for ax, (ia, tag) in zip(axes, [(iaH,'H'), (iaC,'C')]):
        zA, v = lineout_z(rho_grid, ia, zmax_A=10.0)
        ax.plot(zA, v, lw=1.2)
        ax.set_xlabel('z above atom (Å)'); ax.set_ylabel('ρ (e/Å³)')
        ax.set_title(f"ia={ia} Z={atomTypes[ia]} ({tag})")
        ax.grid(True, alpha=0.2)
    plt.tight_layout()
    plt.savefig(os.path.join(_THIS_DIR,'fdbm_rho_lineout_z.png'), dpi=130, bbox_inches='tight'); plt.close()
    print('Saved fdbm_rho_lineout_z.png')

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
    plot_slices(rho_na_grid, 'Neutral-atom promolecule ρ_NA(r) [e/Å³]', 'fdbm_rhoNA_slices.png')
    plot_slices(rho_diff, 'Difference density Δρ=ρ_SCF-ρ_NA [e/Å³]', 'fdbm_rhoDiff_slices.png', sym=True, cmap='bwr')
    plot_slices(V_ES, 'Electrostatic potential from Δρ (Poisson) [eV]', 'fdbm_VES_slices.png', sym=True, cmap='bwr')

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
grads_VES = [np.gradient(V_ES, _step, axis=a).astype(np.float32) for a in range(3)]

print("  Interpolating gradients at scan positions...")
F_pauli = np.zeros((flat_pos.shape[0], 3), dtype=np.float32)
F_vdw   = np.zeros((flat_pos.shape[0], 3), dtype=np.float32)
F_es    = np.zeros((flat_pos.shape[0], 3), dtype=np.float32)
F_es_poiss = np.zeros((flat_pos.shape[0], 3), dtype=np.float32)
# F_Pauli = -A * grad(rho_smooth)  [repulsive density overlap]
for a in range(3):
    F_pauli[:,a] += -args.A_pauli * interp3(grads_rho[a], grid_c)

print("  Atom loop: London vdW + Mulliken electrostatics...")
RA2_VDW = 1.5**2   # vdW softening
RA2_ES  = 2.0**2   # ES softening (avoids divergence, Mulliken at short range)
for ia in range(natoms):
    dr  = flat_pos - atomPos[ia]
    r2  = np.sum(dr**2,axis=1)
    r8_vdw = (r2 + RA2_VDW)**4
    r3_es  = (r2 + RA2_ES)**1.5
    for a in range(3):
        F_vdw[:,a] -= args.C6_CO * 6.0 * dr[:,a] / r8_vdw          # London attractive
        F_es[:,a]  += COULOMB_CONST * q_mulliken[ia] * args.q_CO * dr[:,a] / r3_es   # Mulliken ES

# Poisson ES force (from delta density): F = q_tip * (-∇V)
if args.es_model == 'poisson':
    for a in range(3):
        F_es_poiss[:,a] += -args.q_CO * interp3(grads_VES[a], grid_c)

F_fdbm = np.zeros((flat_pos.shape[0], 4), dtype=np.float32)
F_fdbm[:,:3] = F_pauli + F_vdw + (F_es_poiss if args.es_model=='poisson' else F_es)

FEs_raw = F_fdbm.reshape(nx_s, ny_s, nz, 4)
Fz_raw  = FEs_raw[:,:,:,2]
print(f"Fz_raw: min={Fz_raw.min():.4f}  max={Fz_raw.max():.4f}  mean={Fz_raw.mean():.4f} eV/Å")
np.save(os.path.join(_THIS_DIR,'FEs_raw_fdbm.npy'),FEs_raw)
np.save(os.path.join(_THIS_DIR,'Fz_raw_fdbm.npy'), Fz_raw)

# Component diagnostics: 1D profiles of Pauli/vdW/ES together
if not args.noplot:
    Fz_pauli = F_pauli[:,2].reshape(nx_s, ny_s, nz)
    Fz_vdw   = F_vdw[:,2].reshape(nx_s, ny_s, nz)
    Fz_es    = F_es[:,2].reshape(nx_s, ny_s, nz)
    Fz_esp    = F_es_poiss[:,2].reshape(nx_s, ny_s, nz)
    # pick three representative atom-centered pixels (center/mid/edge)
    centroid = atomPos[:,:2].mean(axis=0)
    dists = np.linalg.norm(atomPos[:,:2]-centroid,axis=1)
    ia_c = int(np.argmin(dists)); ia_e = int(np.argmax(dists))
    ia_m = int(np.argsort(dists)[natoms//2])
    def to_pixel(axy):
        ix=int(round((axy[0]-scan_x0)/scan_dx)); iy=int(round((axy[1]-scan_y0)/scan_dy))
        return max(0,min(nx_s-1,ix)), max(0,min(ny_s-1,iy))
    trace_pixels = [to_pixel(atomPos[ia,:2]) for ia in [ia_c,ia_m,ia_e]]
    trace_labels = ['center-atom','mid-atom','edge-atom']
    fig, axes = plt.subplots(2, 3, figsize=(14, 7))
    fig.suptitle('FDBM Fz components vs height (top=full, bottom=zoom ±0.5)', fontsize=11)
    for col, ((ix,iy), lbl) in enumerate(zip(trace_pixels, trace_labels)):
        for row,(ax,ylim) in enumerate(zip(axes[:,col],[None,(-0.5,0.5)])):
            ax.plot(probe_heights, Fz_pauli[ix,iy,:], color='purple',    lw=1.0, label='Pauli')
            ax.plot(probe_heights, Fz_vdw[ix,iy,:],   color='royalblue', lw=1.0, label='vdW')
            ax.plot(probe_heights, (Fz_esp if args.es_model=='poisson' else Fz_es)[ix,iy,:], color='seagreen',  lw=1.0, label=f"ES[{args.es_model}]")
            ax.plot(probe_heights, Fz_raw[ix,iy,:],   color='k',         lw=1.0, ls='--', label='sum')
            ax.axhline(0,color='k',lw=0.5)
            ax.set_xlabel('Probe height (Å)',fontsize=8); ax.set_ylabel('Fz (eV/Å)',fontsize=8)
            ax.set_title(lbl if row==0 else f"{lbl} [zoom]",fontsize=8)
            ax.legend(fontsize=6); ax.invert_xaxis()
            if ylim:
                ax.set_ylim(ylim); ax.axhspan(ylim[0],0,alpha=0.04,color='blue')
            ax.grid(True, alpha=0.15)
    plt.tight_layout()
    plt.savefig(os.path.join(_THIS_DIR,'fdbm_components_traces_raw.png'), dpi=120, bbox_inches='tight'); plt.close()
    print('Saved fdbm_components_traces_raw.png')

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
