#!/usr/bin/env python3
"""
run_pyocl_fdbm.py - Self-contained PyOpenCL FDBM AFM implementation

Step-by-step implementation with rich debugging and validation at each step.
Uses Fireball only for SCF density matrix; all subsequent operations in PyOpenCL.

Steps:
  1. Density projection to grid (rho_SCF, rho_NA, delta_rho)
  2. Electrostatics via Poisson equation (V from delta_rho)
  3. Pauli repulsion via density convolution
  4. Electrostatics tip-sample convolution
  5. Dispersion (Lennard-Jones C6/r^6)
  6. Final composed calculation with PP relaxation

Input: pentacene.xyz (sample molecule at z=0)
Debug outputs: saved to debug/step{1-6}/ directories
"""
import sys, os, argparse, json, time
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.realpath(os.path.join(_THIS_DIR, '..', '..', '..'))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from pyBall import FireCore as fc
from pyBall.FireballOCL import Grid as ocl_grid
from pyBall.OCL.AFM import (
    fft_poisson, build_gaussian_tip, compute_df,
    plot_slices, plot_grid_Fz, safe_norm
)

# Directories
_FDATA_HCNO = os.path.join(_ROOT, 'tests', 'Fireball', 'Fdata_HCNO')
_FDATA_BASIS = os.path.join(_FDATA_HCNO, 'basis')
_XYZ_PATH = os.path.join(_THIS_DIR, 'pentacene.xyz')
_DEBUG_DIR = os.path.join(_THIS_DIR, 'debug')

# Atomic constants for neutral atom occupations
Z_TO_ZVAL = {1: 1, 6: 4, 7: 5, 8: 6, 16: 6}
RCUT_DEFAULT = {1: 2.3, 6: 2.6, 7: 2.6, 8: 2.5}

# COULOMB_CONST from OCL_AFM.py
COULOMB_CONST = 14.3996448915

def setup_debug_dirs():
    for step in ['step1_density', 'step2_electrostatics', 'step3_pauli',
                 'step4_electrostatics_conv', 'step5_dispersion', 'step5b_lj_comparison', 'step6_composed']:
        os.makedirs(os.path.join(_DEBUG_DIR, step), exist_ok=True)

def load_xyz(fname):
    with open(fname, 'r') as f:
        lines = f.readlines()
    natoms = int(lines[0])
    atomTypes = []; atomPos = []
    for line in lines[2:2+natoms]:
        p = line.split()
        sym = p[0]
        z = 6 if sym == 'C' else 1 if sym == 'H' else 8 if sym == 'O' else 7
        atomTypes.append(z)
        atomPos.append([float(p[1]), float(p[2]), float(p[3])])
    return np.array(atomTypes, dtype=np.int32), np.array(atomPos, dtype=np.float64)

def log_step(step_dir, msg, also_print=True):
    log_path = os.path.join(step_dir, 'log.txt')
    with open(log_path, 'a') as f:
        f.write(msg + '\n')
    if also_print:
        print(msg)

def save_npy(step_dir, name, arr):
    path = os.path.join(step_dir, name)
    np.save(path, arr)
    print(f"  Saved {name} shape={arr.shape} dtype={arr.dtype}")

def plot_field_slices(data, title, fname, step_dir, z_slice=2.0, origin=None, step=None,
                      sym=False, cmap='magma', per_slice_range=True, zvals=None):
    """Plot XY slice at z, XZ slice through center, and line profile.
    
    per_slice_range=True: each image uses its own vmin/vmax based on actual slice data.
    This is essential for XY views where signal is much weaker than global max.
    """
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle(title)
    nx, ny, nz = data.shape
    
    # Compute physical extent if origin/step provided
    extent_xy = None
    extent_xz = None
    if origin is not None and step is not None:
        extent_xy = [float(origin[0]), float(origin[0]) + (nx-1)*step,
                     float(origin[1]), float(origin[1]) + (ny-1)*step]
        extent_xz = [float(origin[0]), float(origin[0]) + (nx-1)*step,
                     float(origin[2]), float(origin[2]) + (nz-1)*step]
    
    # XY slice at z=2.0A (or nearest)
    if origin is not None and step is not None:
        iz = int(np.clip(np.round((z_slice - origin[2]) / step), 0, nz-1))
    else:
        iz = nz // 2
    
    # XY slice data
    slice_xy = data[:, :, iz]
    if per_slice_range:
        vmin_xy, vmax_xy = slice_xy.min(), slice_xy.max()
    else:
        vmin_xy, vmax_xy = data.min(), data.max()
    if sym:
        vabs = max(abs(vmin_xy), abs(vmax_xy), 1e-12)
        vmin_xy, vmax_xy = -vabs, vabs
    im0 = axes[0].imshow(slice_xy.T, origin='lower', cmap=cmap, vmin=vmin_xy, vmax=vmax_xy,
                         extent=extent_xy, aspect='equal')
    axes[0].set_title(f"XY z={z_slice:.1f}Å (iz={iz})  range=[{vmin_xy:.2e},{vmax_xy:.2e}]")
    axes[0].set_xlabel("x [Å]")
    axes[0].set_ylabel("y [Å]")
    fig.colorbar(im0, ax=axes[0])
    
    # XZ slice through center (y=ny//2)
    iy = ny // 2
    slice_xz = data[:, iy, :]
    if per_slice_range:
        vmin_xz, vmax_xz = slice_xz.min(), slice_xz.max()
    else:
        vmin_xz, vmax_xz = data.min(), data.max()
    if sym:
        vabs = max(abs(vmin_xz), abs(vmax_xz), 1e-12)
        vmin_xz, vmax_xz = -vabs, vabs
    im1 = axes[1].imshow(slice_xz.T, origin='lower', cmap=cmap, vmin=vmin_xz, vmax=vmax_xz,
                         extent=extent_xz, aspect='equal')
    axes[1].set_title(f"XZ y-center  range=[{vmin_xz:.2e},{vmax_xz:.2e}]")
    axes[1].set_xlabel("x [Å]")
    axes[1].set_ylabel("z [Å]")
    fig.colorbar(im1, ax=axes[1])
    
    # Line profile above center atom
    ix, iy_c = nx // 2, ny // 2
    zlabel = np.arange(nz) if zvals is None else zvals
    axes[2].plot(zlabel, data[ix, iy_c, :], 'k-', lw=1.0)
    axes[2].set_title(f"Line profile at center (ix={ix},iy={iy_c})")
    if zvals is not None:
        axes[2].set_xlabel("z [Å]")
    else:
        axes[2].set_xlabel("z index")
    axes[2].set_ylabel("Value")
    axes[2].grid(True, alpha=0.3)
    
    plt.tight_layout()
    out = os.path.join(step_dir, fname)
    plt.savefig(out, dpi=120, bbox_inches='tight'); plt.close()
    print(f"  Saved {fname}")

def check_field_sanity(field, name, min_abs=1e-6, max_abs=1e+6):
    """Raise ValueError if field range is suspiciously small or large."""
    fmin, fmax = float(field.min()), float(field.max())
    fmaxabs = max(abs(fmin), abs(fmax))
    if fmaxabs < min_abs:
        raise ValueError(f"{name}: range [{fmin:.3e},{fmax:.3e}] near zero (<{min_abs})")
    if fmaxabs > max_abs:
        raise ValueError(f"{name}: range [{fmin:.3e},{fmax:.3e}] exceeds threshold (>{max_abs})")
    return True

def plot_line_profile(data, title, fname, step_dir, ix=None, iy=None):
    """Plot line profile along z at (ix, iy)."""
    nx, ny, nz = data.shape
    if ix is None: ix = nx // 2
    if iy is None: iy = ny // 2
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(data[ix, iy, :], 'k-', lw=1.5)
    ax.set_title(f"{title} (ix={ix}, iy={iy})")
    ax.set_xlabel("z index")
    ax.set_ylabel("Value")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    out = os.path.join(step_dir, fname)
    plt.savefig(out, dpi=120, bbox_inches='tight'); plt.close()
    print(f"  Saved {fname}")

# ═══════════════════════════════════════════════════════════════════════════════
# Step 1: Density Projection
# ═══════════════════════════════════════════════════════════════════════════════

def step1_density_projection(atomTypes, atomPos, step=0.15, margin=4.0, z_extra=6.0,
                            block=8, nscf=200, fdata_dir=_FDATA_HCNO,
                            fdata_basis=_FDATA_BASIS, verbosity=0, stop_after=False):
    step_dir = os.path.join(_DEBUG_DIR, 'step1_density')
    log = lambda msg: log_step(step_dir, msg)
    
    log("=" * 60)
    log("STEP 1: Density Projection to Grid")
    log("=" * 60)
    
    # 1a. Run Fireball SCF
    log("\n1a. Running Fireball SCF...")
    _fball_cwd = os.path.join(_ROOT, 'tests', 'pyFireball')
    orig_cwd = os.getcwd()
    os.chdir(_fball_cwd)
    
    fdata_local = os.path.join(_fball_cwd, 'Fdata')
    if not os.path.exists(os.path.join(fdata_local, 'info.dat')):
        if os.path.lexists(fdata_local):
            os.unlink(fdata_local)
        os.symlink(fdata_dir, fdata_local)
        log(f"Created Fdata symlink -> {fdata_dir}")
    
    fc.setVerbosity(verbosity)
    fc.preinit()
    fc.init(atomTypes, atomPos)
    fc.SCF(atomPos, nmax_scf=nscf)
    
    dims = fc.get_HS_dims()
    neighs = fc.get_HS_neighs(dims)
    neighs = fc.get_rho_sparse(dims, data=neighs)
    rho_sparse = neighs.rho
    log(f"rho_sparse shape={rho_sparse.shape} max_abs={np.abs(rho_sparse).max():.4f}")
    
    # Mulliken charges
    charges2d = np.zeros((len(atomTypes), 1), dtype=np.float64)
    fc.getCharges(charges2d)
    q_pop = charges2d[:, 0]
    z_val = np.array([Z_TO_ZVAL.get(int(z), int(z)) for z in atomTypes], dtype=float)
    q_mulliken = (q_pop - z_val).astype(np.float32)
    log(f"Mulliken charges: sum={q_mulliken.sum():.4f} range=[{q_mulliken.min():.3f},{q_mulliken.max():.3f}]")
    
    os.chdir(orig_cwd)
    
    # 1b. Setup grid
    log("\n1b. Setting up density grid...")
    pos_min = atomPos.min(axis=0) - margin
    pos_max = atomPos.max(axis=0) + np.array([margin, margin, margin + z_extra])
    span = pos_max - pos_min
    ngrid = (np.ceil(np.ceil(span / step) / block).astype(int) * block)
    total_span = ngrid * step
    origin = (0.5*(pos_min + pos_max) - 0.5*total_span).astype(np.float32)
    grid_spec = {
        'origin': origin,
        'dA': [step, 0., 0.], 'dB': [0., step, 0.], 'dC': [0., 0., step],
        'ngrid': ngrid.astype(int),
    }
    log(f"Grid: origin={origin.round(2)} ngrid={ngrid} span={total_span.round(2)}")
    
    # 1c. Build neutral atom density matrix
    log("\n1c. Building neutral-atom density matrix...")
    rho_na = np.zeros_like(rho_sparse, dtype=np.float32)
    neigh_j = neighs.neigh_j.reshape(len(atomTypes), -1)
    for i in range(len(atomTypes)):
        slots = np.where(neigh_j[i] == (i+1))[0]
        if len(slots) == 0:
            raise RuntimeError(f"No self-neighbor for atom i={i}")
        iself = int(slots[0])
        occ = _onsite_occ(int(atomTypes[i]))
        rho_na[i, iself, :, :] = 0.0
        for k in range(4):
            rho_na[i, iself, k, k] = occ[k]
    log("Neutral-atom density matrix built.")
    
    # 1d. Project to grid
    log("\n1d. Projecting densities to grid...")
    projector = ocl_grid.GridProjector(fdata_dir=fdata_basis, verbosity=verbosity)
    projector.load_basis(sorted(set(atomTypes.tolist())))
    atoms_dict = {
        'pos': atomPos,
        'Rcut': np.array([RCUT_DEFAULT.get(int(z), 4.5) for z in atomTypes]),
        'type': atomTypes,
    }
    
    log("  Projecting SCF density...")
    rho_grid = projector.project(rho_sparse, neighs, atoms_dict, grid_spec, nMaxAtom=64, use_tiled=True)
    dV = step**3
    q_rho = rho_grid.sum() * dV
    log(f"  rho_grid: shape={rho_grid.shape} range=[{rho_grid.min():.5f},{rho_grid.max():.5f}]")
    log(f"  Integrated electrons (rho_SCF) = {q_rho:.2f} (expected: {z_val.sum():.1f})")
    
    log("  Projecting neutral-atom density...")
    rho_na_grid = projector.project(rho_na, neighs, atoms_dict, grid_spec, nMaxAtom=64, use_tiled=True)
    q_na = rho_na_grid.sum() * dV
    log(f"  rho_NA: shape={rho_na_grid.shape} range=[{rho_na_grid.min():.5f},{rho_na_grid.max():.5f}]")
    log(f"  Integrated electrons (rho_NA) = {q_na:.2f} (expected: {z_val.sum():.1f})")
    
    rho_diff = (rho_grid - rho_na_grid).astype(np.float32)
    q_diff = rho_diff.sum() * dV
    log(f"  delta_rho: range=[{rho_diff.min():.5f},{rho_diff.max():.5f}]")
    log(f"  Integrated delta charge = {q_diff:.4f} (should be ~0)")
    
    # Save outputs
    save_npy(step_dir, 'rho_grid.npy', rho_grid)
    save_npy(step_dir, 'rho_na_grid.npy', rho_na_grid)
    save_npy(step_dir, 'rho_diff.npy', rho_diff)
    np.save(os.path.join(step_dir, 'origin.npy'), origin)
    np.save(os.path.join(step_dir, 'ngrid.npy'), ngrid)
    
    # Plot diagnostics
    plot_field_slices(rho_grid, "SCF Density [e/Å³]", "step1_rho_slices.png", step_dir,
                      z_slice=2.0, origin=origin, step=step)
    plot_field_slices(rho_na_grid, "Neutral Atom Density [e/Å³]", "step1_rhoNA_slices.png", step_dir,
                      z_slice=2.0, origin=origin, step=step)
    plot_field_slices(rho_diff, "Delta Density = rho_SCF - rho_NA [e/Å³]", "step1_rhoDiff_slices.png",
                      step_dir, z_slice=2.0, origin=origin, step=step, sym=True, cmap='bwr')
    
    # Max projections
    plot_slices(rho_grid, "SCF Density Max Projections", "step1_rho_maxproj.png", save_dir=step_dir)
    
    # Log validation
    log("\n--- Validation ---")
    log(f"Charge conservation check: {abs(q_diff):.4f} (target: < 0.1)")
    if abs(q_diff) > 0.5:
        log("WARNING: Charge neutrality violated! Check projection settings.")
    
    # Sanity checks
    check_field_sanity(rho_grid, "rho_grid", min_abs=1e-12, max_abs=1e6)
    check_field_sanity(rho_diff, "rho_diff", min_abs=1e-12, max_abs=1e6)
    
    log("\nStep 1 COMPLETE.")
    
    if stop_after:
        log("STOP_AFTER_STEP1 set. Exiting.")
        sys.exit(0)
    
    return rho_grid, rho_na_grid, rho_diff, origin, ngrid, grid_spec, q_mulliken, atomPos, atomTypes

def _onsite_occ(Z):
    if Z == 1:  return np.array([1.0, 0.0, 0.0, 0.0], dtype=np.float32)
    if Z == 6:  return np.array([2.0, 2.0/3, 2.0/3, 2.0/3], dtype=np.float32)
    if Z == 7:  return np.array([2.0, 1.0, 1.0, 1.0], dtype=np.float32)
    if Z == 8:  return np.array([2.0, 4.0/3, 4.0/3, 4.0/3], dtype=np.float32)
    return np.array([float(Z_TO_ZVAL.get(int(Z), int(Z))), 0.0, 0.0, 0.0], dtype=np.float32)

def get_co_tip_density(grid_spec, fdata_dir=_FDATA_HCNO, fdata_basis=_FDATA_BASIS,
                       step=0.15, nscf=50, verbosity=0):
    """
    Run Fireball SCF on CO molecule and project density to the given grid_spec.
    Uses subprocess because Fireball cannot reallocate for different molecules.
    Returns (rho_total, rho_delta, projector) for use as tip density.
    
    O is placed at grid origin (probe position). C is at z=1.13A.
    """
    co_dir = os.path.join(_DEBUG_DIR, 'co_tip')
    os.makedirs(co_dir, exist_ok=True)
    
    # Check cache
    cache_total = os.path.join(co_dir, 'co_rho_total.npy')
    cache_delta = os.path.join(co_dir, 'co_rho_delta.npy')
    cache_meta = os.path.join(co_dir, 'co_meta.json')
    
    # Build cache key from grid_spec + nscf + step
    grid_key = {
        'ngrid': grid_spec['ngrid'].tolist(),
        'origin': grid_spec['origin'].tolist(),
        'step': step,
        'nscf': nscf,
    }
    
    use_cache = False
    if os.path.exists(cache_total) and os.path.exists(cache_delta) and os.path.exists(cache_meta):
        with open(cache_meta, 'r') as f:
            cached_key = json.load(f)
        if cached_key == grid_key:
            use_cache = True
            print(f"\n[CO Tip] Using cached CO densities from {co_dir}")
    
    if not use_cache:
        print(f"\n[CO Tip] Computing CO density via subprocess (nscf={nscf})...")
        # Serialize grid_spec for subprocess
        def _tolist(x):
            return x.tolist() if hasattr(x, 'tolist') else list(x)
        grid_spec_json = json.dumps({
            'origin': _tolist(grid_spec['origin']),
            'dA': _tolist(grid_spec['dA']),
            'dB': _tolist(grid_spec['dB']),
            'dC': _tolist(grid_spec['dC']),
            'ngrid': _tolist(grid_spec['ngrid']),
        })
        
        script = os.path.join(_THIS_DIR, 'compute_co_tip.py')
        cmd = [
            sys.executable, script,
            co_dir, grid_spec_json, str(step), str(nscf),
            fdata_dir, fdata_basis
        ]
        env = os.environ.copy()
        env['PYTHONPATH'] = _ROOT + ':' + env.get('PYTHONPATH', '')
        
        import subprocess
        result = subprocess.run(cmd, capture_output=True, text=True, env=env, cwd=_THIS_DIR)
        if result.returncode != 0:
            print("STDOUT:", result.stdout)
            print("STDERR:", result.stderr)
            raise RuntimeError(f"CO tip subprocess failed with code {result.returncode}")
        print(result.stdout)
        
        # Save cache key
        with open(cache_meta, 'w') as f:
            json.dump(grid_key, f)
    
    rho_total_co = np.load(cache_total)
    rho_delta_co = np.load(cache_delta)
    print(f"  CO rho_total: shape={rho_total_co.shape} range=[{rho_total_co.min():.4f},{rho_total_co.max():.4f}]")
    print(f"  CO rho_delta: shape={rho_delta_co.shape} range=[{rho_delta_co.min():.4f},{rho_delta_co.max():.4f}]")
    
    return rho_total_co, rho_delta_co, None

# ═══════════════════════════════════════════════════════════════════════════════
# Step 2: Electrostatics via Poisson
# ═══════════════════════════════════════════════════════════════════════════════

def step2_electrostatics(rho_diff, step, origin, ngrid, stop_after=False):
    step_dir = os.path.join(_DEBUG_DIR, 'step2_electrostatics')
    log = lambda msg: log_step(step_dir, msg)
    
    log("=" * 60)
    log("STEP 2: Electrostatics via Poisson Equation")
    log("=" * 60)
    
    nx, ny, nz = rho_diff.shape
    log(f"Input delta_rho: shape={rho_diff.shape}")
    
    # FFT Poisson solver
    log("\n2a. Solving Poisson equation via FFT...")
    rho_k = np.fft.fftn(rho_diff)
    kx = 2*np.pi * np.fft.fftfreq(nx, d=step)
    ky = 2*np.pi * np.fft.fftfreq(ny, d=step)
    kz = 2*np.pi * np.fft.fftfreq(nz, d=step)
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij')
    k2 = KX**2 + KY**2 + KZ**2
    k2[0, 0, 0] = 1.0  # avoid division by zero
    V_k = 4.0 * np.pi * COULOMB_CONST * rho_k / k2
    V_k[0, 0, 0] = 0.0
    V_ES = np.real(np.fft.ifftn(V_k)).astype(np.float32)
    
    log(f"V_ES: shape={V_ES.shape} range=[{V_ES.min():.3f},{V_ES.max():.3f}] eV")
    
    # Validation: verify Poisson equation by finite difference Laplacian
    log("\n2b. Validating Poisson solution...")
    from scipy.ndimage import laplace
    lapV = laplace(V_ES.astype(np.float64)) / (step**2)
    rhs = -4.0 * np.pi * COULOMB_CONST * rho_diff.astype(np.float64)
    residual = lapV - rhs
    log(f"  Laplacian(V) residual: max={np.abs(residual).max():.3f} mean={np.abs(residual).mean():.6f}")
    
    # Save outputs
    save_npy(step_dir, 'V_ES.npy', V_ES)
    
    # Plot diagnostics
    plot_field_slices(V_ES, "Electrostatic Potential [eV]", "step2_VES_slices.png",
                      step_dir, z_slice=2.0, origin=origin, step=step, sym=True, cmap='bwr')
    plot_line_profile(V_ES, "V_ES line profile", "step2_VES_lineprofile.png", step_dir)
    
    # Poisson check plot - symmetric ranges per subplot
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle("Poisson Equation Validation")
    
    def _symshow(ax, data_slice, title):
        vabs = max(abs(float(data_slice.min())), abs(float(data_slice.max())), 1e-12)
        im = ax.imshow(data_slice.T, origin='lower', cmap='bwr', vmin=-vabs, vmax=vabs)
        ax.set_title(f"{title}\nrange=[{-vabs:.2e},{vabs:.2e}]")
        fig.colorbar(im, ax=ax)
        return im
    
    # Find z-slice where |rho| is maximal (molecule plane), not nz//2
    iz_rho = int(np.argmax(np.sum(np.abs(rho_diff), axis=(0,1))))
    z_rho = origin[2] + iz_rho * step
    log(f"  Poisson check at z={z_rho:.2f}Å (iz={iz_rho}) where |rho| is maximal")
    
    lapV_slice = lapV[:, :, iz_rho]
    rhs_slice  = rhs[:, :, iz_rho]
    # Assert that ranges are physically consistent (within 20% of each other)
    vabs_lapV = max(abs(float(lapV_slice.min())), abs(float(lapV_slice.max())))
    vabs_rhs  = max(abs(float(rhs_slice.min())),  abs(float(rhs_slice.max())))
    ratio = max(vabs_lapV, vabs_rhs) / min(vabs_lapV, vabs_rhs) if min(vabs_lapV, vabs_rhs) > 1e-12 else 1.0
    assert ratio < 1.2, f"Poisson validation FAILED: |lapV|={vabs_lapV:.3e} vs |rhs|={vabs_rhs:.3e} (ratio={ratio:.2f}). Check units/step scaling."
    log(f"  Poisson range check: |lapV|={vabs_lapV:.3e} |rhs|={vabs_rhs:.3e} ratio={ratio:.3f}  [OK]")
    
    _symshow(axes[0], lapV_slice, "Laplacian(V) [eV/Å²]")
    _symshow(axes[1], rhs_slice, f"-4πρ [eV/Å²] at z={z_rho:.1f}Å")
    _symshow(axes[2], residual[:, :, iz_rho], f"Residual (max_err={np.abs(residual).max():.2e})")
    
    plt.tight_layout()
    plt.savefig(os.path.join(step_dir, "step2_VES_poisson_check.png"), dpi=120, bbox_inches='tight')
    plt.close()
    print("  Saved step2_VES_poisson_check.png")
    
    # Sanity check
    check_field_sanity(V_ES, "V_ES", min_abs=1e-6, max_abs=1e6)
    
    log("\nStep 2 COMPLETE.")
    if stop_after:
        log("STOP_AFTER_STEP2 set. Exiting.")
        sys.exit(0)
    return V_ES

# ═══════════════════════════════════════════════════════════════════════════════
# Step 3: Pauli Repulsion via Density Convolution
# ═══════════════════════════════════════════════════════════════════════════════

def step3_pauli(rho_grid, origin, step, scan_xs, scan_ys, probe_heights,
                A_pauli=3.0, sigma_tip=0.7, tip_density=None, stop_after=False,
                beta_pauli=1.0):
    """
    Compute Pauli repulsion energy field via FFT convolution.
    
    tip_density: if provided, use as tip density (e.g., true CO total density).
                 If None, use Gaussian tip with sigma_tip.
    beta_pauli: exponent applied to rho_grid before convolution (default 1.0).
                Values > 1.0 sharpen the repulsive core; < 1.0 soften it.
    """
    step_dir = os.path.join(_DEBUG_DIR, 'step3_pauli')
    log = lambda msg: log_step(step_dir, msg)
    
    log("=" * 60)
    log("STEP 3: Pauli Repulsion via Density Convolution")
    log("=" * 60)
    
    nx, ny, nz = rho_grid.shape
    dV = step**3
    
    # Build tip density
    if tip_density is not None:
        log("\n3a. Using TRUE CO density as tip (from Fireball SCF)...")
        # FFT convolution requires tip at array index 0. CO was projected to grid center.
        # Circularly shift from grid center to index 0.
        nx_t, ny_t, nz_t = tip_density.shape
        cx_t, cy_t, cz_t = nx_t // 2, ny_t // 2, nz_t // 2
        tip_kernel = np.roll(np.roll(np.roll(tip_density.astype(np.float64), -cx_t, axis=0), -cy_t, axis=1), -cz_t, axis=2)
        log(f"  Tip density: shape={tip_kernel.shape} range=[{tip_kernel.min():.4f},{tip_kernel.max():.4f}]")
        log(f"  Shifted from grid center ({cx_t},{cy_t},{cz_t}) to index 0 for FFT")
    else:
        log("\n3a. Building Gaussian tip density kernel...")
        tip_kernel = build_gaussian_tip((nx, ny, nz), step, sigma_tip)
    
    # FFT convolution: E_pauli = A * convolution(rho_sample^beta, rho_tip^beta)
    log("\n3b. Computing Pauli energy via FFT convolution...")
    log(f"  A_pauli={A_pauli}  beta_pauli={beta_pauli}")
    if abs(beta_pauli - 1.0) > 1e-6:
        rho_eff = np.power(rho_grid.astype(np.float64), beta_pauli)
        tip_eff = np.power(np.maximum(tip_kernel, 0.0), beta_pauli)  # tip density is non-negative
        log(f"  rho^beta:  range=[{rho_eff.min():.6f},{rho_eff.max():.6f}]")
        log(f"  tip^beta:  range=[{tip_eff.min():.6f},{tip_eff.max():.6f}]")
    else:
        rho_eff = rho_grid.astype(np.float64)
        tip_eff = tip_kernel
    rho_k = np.fft.fftn(rho_eff)
    tip_k = np.fft.fftn(tip_eff)
    E_pauli_field = A_pauli * dV * np.real(np.fft.ifftn(rho_k * tip_k)).astype(np.float32)
    
    log(f"E_Pauli: shape={E_pauli_field.shape} range=[{E_pauli_field.min():.3f},{E_pauli_field.max():.3f}] eV")
    
    # Compute gradient (forces) via FFT
    log("\n3c. Computing Pauli forces via FFT gradient...")
    grads = []
    for a, freq in enumerate([np.fft.fftfreq(nx, d=step), np.fft.fftfreq(ny, d=step), np.fft.fftfreq(nz, d=step)]):
        k = 2*np.pi * freq
        shape = [1, 1, 1]
        shape[a] = len(k)
        k_arr = k.reshape(shape)
        grad_k = 1j * k_arr * np.fft.fftn(E_pauli_field.astype(np.float64))
        grad = np.real(np.fft.ifftn(grad_k)).astype(np.float32)
        grads.append(grad)
    grads_E_Pauli = np.stack(grads, axis=-1)  # (nx,ny,nz,3)
    
    log(f"F_Pauli: range=[{grads_E_Pauli.min():.3f},{grads_E_Pauli.max():.3f}] eV/Å")
    
    # Save outputs
    save_npy(step_dir, 'E_Pauli_field.npy', E_pauli_field)
    save_npy(step_dir, 'grads_E_Pauli.npy', grads_E_Pauli)
    
    # Plot diagnostics
    plot_field_slices(E_pauli_field, f"Pauli Energy [eV]  (A={A_pauli}, beta={beta_pauli})", "step3_Epauli_slices.png",
                      step_dir, z_slice=2.0, origin=origin, step=step)
    
    # Fz component
    Fz_pauli = -grads_E_Pauli[:, :, :, 2]
    plot_field_slices(Fz_pauli, f"Pauli Fz [eV/Å]  (A={A_pauli}, beta={beta_pauli})", "step3_Fz_pauli_slices.png",
                      step_dir, z_slice=2.0, origin=origin, step=step, sym=True, cmap='bwr')
    
    # Sanity checks
    check_field_sanity(E_pauli_field, "E_pauli_field", min_abs=1e-6, max_abs=1e6)
    check_field_sanity(grads_E_Pauli, "grads_E_Pauli", min_abs=1e-6, max_abs=1e6)
    
    log("\nStep 3 COMPLETE.")
    if stop_after:
        log("STOP_AFTER_STEP3 set. Exiting.")
        sys.exit(0)
    return E_pauli_field, grads_E_Pauli

# ═══════════════════════════════════════════════════════════════════════════════
# Step 4: Electrostatics Tip-Sample Convolution
# ═══════════════════════════════════════════════════════════════════════════════

def step4_electrostatics_conv(V_ES, origin, step, scan_xs, scan_ys, probe_heights,
                             q_CO=-0.05, sigma_tip=0.7, tip_charge_density=None, stop_after=False):
    """
    Compute electrostatic tip-sample interaction via FFT convolution.
    
    tip_charge_density: if provided, use as tip charge distribution (e.g., CO delta density).
                        If None, use Gaussian scaled by q_CO.
    """
    step_dir = os.path.join(_DEBUG_DIR, 'step4_electrostatics_conv')
    log = lambda msg: log_step(step_dir, msg)
    
    log("=" * 60)
    log("STEP 4: Electrostatics Tip-Sample Convolution")
    log("=" * 60)
    
    nx, ny, nz = V_ES.shape
    dV = step**3
    
    # Build tip charge density
    if tip_charge_density is not None:
        log("\n4a. Using TRUE CO delta-density as tip charge...")
        # FFT convolution requires tip at array index 0. CO was projected to grid center.
        # Circularly shift from grid center to index 0.
        nx_t, ny_t, nz_t = tip_charge_density.shape
        cx_t, cy_t, cz_t = nx_t // 2, ny_t // 2, nz_t // 2
        tip_charge = np.roll(np.roll(np.roll(tip_charge_density.astype(np.float64), -cx_t, axis=0), -cy_t, axis=1), -cz_t, axis=2)
        log(f"  Tip charge: shape={tip_charge.shape} range=[{tip_charge.min():.4f},{tip_charge.max():.4f}]")
        log(f"  Integral (should be ~0 for neutral CO): {tip_charge.sum() * dV:.4f} e")
        log(f"  Shifted from grid center ({cx_t},{cy_t},{cz_t}) to index 0 for FFT")
    else:
        log("\n4a. Building Gaussian tip charge density...")
        tip_kernel = build_gaussian_tip((nx, ny, nz), step, sigma_tip)
        tip_charge = q_CO * tip_kernel  # charged tip density
    
    # Convolution: E_es = convolution(V_ES, rho_tip)
    log("\n4b. Computing electrostatic energy via convolution...")
    V_k = np.fft.fftn(V_ES.astype(np.float64))
    tip_k = np.fft.fftn(tip_charge)
    E_ES_field = dV * np.real(np.fft.ifftn(V_k * tip_k)).astype(np.float32)
    
    log(f"E_ES: shape={E_ES_field.shape} range=[{E_ES_field.min():.3f},{E_ES_field.max():.3f}] eV")
    
    # Compute gradient (forces)
    log("\n4c. Computing electrostatic forces...")
    grads = []
    for a, freq in enumerate([np.fft.fftfreq(nx, d=step), np.fft.fftfreq(ny, d=step), np.fft.fftfreq(nz, d=step)]):
        k = 2*np.pi * freq
        shape = [1, 1, 1]
        shape[a] = len(k)
        k_arr = k.reshape(shape)
        grad_k = 1j * k_arr * np.fft.fftn(E_ES_field.astype(np.float64))
        grad = np.real(np.fft.ifftn(grad_k)).astype(np.float32)
        grads.append(grad)
    grads_E_ES = np.stack(grads, axis=-1)
    
    log(f"F_ES: range=[{grads_E_ES.min():.3f},{grads_E_ES.max():.3f}] eV/Å")
    
    # Save outputs
    save_npy(step_dir, 'E_ES_field.npy', E_ES_field)
    save_npy(step_dir, 'grads_E_ES.npy', grads_E_ES)
    
    # Plot diagnostics
    plot_field_slices(E_ES_field, "Electrostatic Energy [eV]", "step4_EES_slices.png",
                      step_dir, z_slice=2.0, origin=origin, step=step, sym=True, cmap='bwr')
    Fz_es = -grads_E_ES[:, :, :, 2]
    plot_field_slices(Fz_es, "Electrostatic Fz [eV/Å]", "step4_Fz_ES_slices.png",
                      step_dir, z_slice=2.0, origin=origin, step=step, sym=True, cmap='bwr')
    
    # Sanity checks
    check_field_sanity(E_ES_field, "E_ES_field", min_abs=1e-6, max_abs=1e6)
    check_field_sanity(grads_E_ES, "grads_E_ES", min_abs=1e-6, max_abs=1e6)
    
    log("\nStep 4 COMPLETE.")
    if stop_after:
        log("STOP_AFTER_STEP4 set. Exiting.")
        sys.exit(0)
    return E_ES_field, grads_E_ES

# ═══════════════════════════════════════════════════════════════════════════════
# Step 5: Dispersion (Lennard-Jones C6/r^6)
# ═══════════════════════════════════════════════════════════════════════════════

def step5_dispersion(atomPos, atomTypes, origin, step, ngrid,
                     C6_CO=30.0, mol_z=0.0, stop_after=False):
    """Compute vdW on the density grid (same grid as Pauli/ES), then interpolate to scan in step6."""
    step_dir = os.path.join(_DEBUG_DIR, 'step5_dispersion')
    log = lambda msg: log_step(step_dir, msg)
    
    log("=" * 60)
    log("STEP 5: Dispersion (C6/r^6) on density grid")
    log("=" * 60)
    
    nx, ny, nz = int(ngrid[0]), int(ngrid[1]), int(ngrid[2])
    natoms = len(atomPos)
    
    # Build C6 per atom (simplified: C6 for C=15, H=1, tip_C6=30)
    C6_atom = np.array([15.0 if z == 6 else 1.0 if z == 1 else 10.0 for z in atomTypes])
    
    log(f"\n5a. Computing pairwise C6/r^6 dispersion on {nx}x{ny}x{nz} density grid...")
    log(f"  Atoms: {natoms}, C6 range=[{C6_atom.min():.1f},{C6_atom.max():.1f}]")
    
    # Density grid coordinates
    xs = origin[0] + np.arange(nx) * step
    ys = origin[1] + np.arange(ny) * step
    zs = origin[2] + np.arange(nz) * step
    XX, YY, ZZ = np.meshgrid(xs, ys, zs, indexing='ij')
    
    E_vdw = np.zeros((nx, ny, nz), dtype=np.float32)
    
    # Use regularized C6/r^6 with RA2 softening
    RA2_VDW = 1.5**2
    log(f"  Using vdW regularization: RA2={RA2_VDW:.2f} Å²")
    
    # Vectorized computation per atom
    for ia in range(natoms):
        ap = atomPos[ia]
        C6_eff = np.sqrt(C6_atom[ia] * C6_CO)
        dx = XX - ap[0]; dy = YY - ap[1]; dz = ZZ - ap[2]
        r2 = dx**2 + dy**2 + dz**2
        E_vdw -= C6_eff / (r2 + RA2_VDW)**3  # -C6/(r^2+RA2)^3
    
    log(f"E_vdw: shape={E_vdw.shape} range=[{E_vdw.min():.3f},{E_vdw.max():.3f}] eV")
    
    # Compute energy gradients on density grid (same step in all directions)
    log("\n5b. Computing vdW energy gradients on density grid...")
    grads_vdw = [np.gradient(E_vdw, step, axis=0),
                 np.gradient(E_vdw, step, axis=1),
                 np.gradient(E_vdw, step, axis=2)]
    grads_E_vdw = np.stack(grads_vdw, axis=-1).astype(np.float32)  # ∇E
    
    log(f"grads_E_vdw: range=[{grads_E_vdw.min():.3f},{grads_E_vdw.max():.3f}] eV/Å")
    
    # Save outputs
    save_npy(step_dir, 'E_vdw_field.npy', E_vdw)
    save_npy(step_dir, 'grads_E_vdw.npy', grads_E_vdw)
    
    # Plot diagnostics - slice at molecule z-plane (where atoms are)
    mol_z_mid = float(np.median(atomPos[:, 2]))
    plot_field_slices(E_vdw, "Dispersion Energy [eV]", "step5_Evdw_slices.png",
                      step_dir, z_slice=mol_z_mid, origin=origin, step=step, cmap='viridis')
    
    # Sanity checks
    check_field_sanity(E_vdw, "E_vdw", min_abs=1e-6, max_abs=1e6)
    check_field_sanity(grads_E_vdw, "grads_E_vdw", min_abs=1e-6, max_abs=1e6)
    
    log("\nStep 5 COMPLETE.")
    if stop_after:
        log("STOP_AFTER_STEP5 set. Exiting.")
        sys.exit(0)
    return E_vdw, grads_E_vdw

# ═══════════════════════════════════════════════════════════════════════════════
# Step 5b: LJ/Morse Comparison (for validation against FDBM)
# ═══════════════════════════════════════════════════════════════════════════════

def compute_lj_forces(scan_xs, scan_ys, probe_heights, atomPos, atomTypes, q_mulliken,
                      q_CO=-0.05, sigma_tip=3.0, eps_tip=0.003, r_cut=12.0):
    """
    Compute Lennard-Jones + point-charge forces on scan grid for side-by-side
    comparison with FDBM decomposition.
    
    Returns Fz on scan grid for: repulsive (Pauli-like), attractive (London-like),
    electrostatic, and total.
    """
    step_dir = os.path.join(_DEBUG_DIR, 'step5b_lj_comparison')
    log = lambda msg: log_step(step_dir, msg)
    log("=" * 60)
    log("STEP 5b: LJ + Point-Charge Comparison")
    log("=" * 60)
    
    nx_s, ny_s, nz_s = len(scan_xs), len(scan_ys), len(probe_heights)
    
    # LJ parameters (UFF-like, in eV and Å)
    # Tip apex modeled as oxygen-like atom
    lj_params = {
        1:  (2.89, 0.0025),   # H: sigma, epsilon
        6:  (3.40, 0.00377),  # C
        7:  (3.26, 0.00376),  # N
        8:  (3.00, 0.00296),  # O
    }
    
    sigma_tip_A = float(sigma_tip)
    eps_tip_A = float(eps_tip)
    
    # Pre-allocate
    Fz_rep = np.zeros((nx_s, ny_s, nz_s), dtype=np.float32)
    Fz_att = np.zeros((nx_s, ny_s, nz_s), dtype=np.float32)
    Fz_es_lj = np.zeros((nx_s, ny_s, nz_s), dtype=np.float32)
    
    log(f"Scan grid: {nx_s}x{ny_s}x{nz_s}")
    log(f"Atoms: {len(atomPos)}  Tip: sigma={sigma_tip_A:.2f}Å eps={eps_tip_A:.4f}eV q={q_CO:.3f}e")
    
    XX, YY, ZZ = np.meshgrid(scan_xs, scan_ys, probe_heights, indexing='ij')
    
    for ia, (Z, pos, q) in enumerate(zip(atomTypes, atomPos, q_mulliken)):
        sigma_s, eps_s = lj_params.get(int(Z), (3.40, 0.00377))
        # Lorentz-Berthelot mixing
        sigma = 0.5 * (sigma_s + sigma_tip_A)
        eps = np.sqrt(eps_s * eps_tip_A)
        
        dx = XX - float(pos[0])
        dy = YY - float(pos[1])
        dz = ZZ - float(pos[2])
        r2 = dx*dx + dy*dy + dz*dz
        r2 = np.maximum(r2, 0.25)  # regularization: min r=0.5Å
        r = np.sqrt(r2)
        
        sr6 = (sigma / r)**6
        sr12 = sr6 * sr6
        
        # LJ force z-component: Fz = (48*eps*sr12 - 24*eps*sr6) * dz / r^2
        fac = (48.0 * eps * sr12 - 24.0 * eps * sr6) / r2
        Fz_total_lj = fac * dz
        
        # Decompose
        Fz_rep += (48.0 * eps * sr12 / r2 * dz).astype(np.float32)
        Fz_att += (-24.0 * eps * sr6 / r2 * dz).astype(np.float32)
        
        # Electrostatic: Fz = q_atom * q_tip * 14.4 * dz / r^3
        if abs(q_CO) > 1e-12 and abs(q) > 1e-12:
            Fz_es_lj += (q * q_CO * COULOMB_CONST * dz / (r2 * r)).astype(np.float32)
    
    Fz_total = Fz_rep + Fz_att + Fz_es_lj
    
    log(f"Fz_rep:   range=[{Fz_rep.min():.4f},{Fz_rep.max():.4f}] eV/Å")
    log(f"Fz_att:   range=[{Fz_att.min():.4f},{Fz_att.max():.4f}] eV/Å")
    log(f"Fz_es:    range=[{Fz_es_lj.min():.4f},{Fz_es_lj.max():.4f}] eV/Å")
    log(f"Fz_total: range=[{Fz_total.min():.4f},{Fz_total.max():.4f}] eV/Å")
    
    # Save
    save_npy(step_dir, 'Fz_lj_repulsive.npy', Fz_rep)
    save_npy(step_dir, 'Fz_lj_attractive.npy', Fz_att)
    save_npy(step_dir, 'Fz_lj_es.npy', Fz_es_lj)
    save_npy(step_dir, 'Fz_lj_total.npy', Fz_total)
    
    log("Step 5b COMPLETE.")
    return Fz_rep, Fz_att, Fz_es_lj, Fz_total

# ═══════════════════════════════════════════════════════════════════════════════
# Step 6: Final Composed Calculation
# ═══════════════════════════════════════════════════════════════════════════════

def step6_composed(E_pauli, grads_pauli, E_es, grads_es, E_vdw, grads_vdw,
                   origin, step, scan_xs, scan_ys, probe_heights, atomPos,
                   A_pauli=3.0, beta_pauli=1.0, q_CO=-0.05, stop_after=False,
                   lj_forces=None):
    """
    Compose all force components on the SCAN grid.
    
    Steps 3-5 compute Pauli/ES/vdW on density grid (nx_d,ny_d,nz_d).
    We interpolate all three to scan positions before composition.
    """
    from scipy.ndimage import map_coordinates
    step_dir = os.path.join(_DEBUG_DIR, 'step6_composed')
    log = lambda msg: log_step(step_dir, msg)
    
    log("=" * 60)
    log("STEP 6: Final Composed Calculation")
    log("=" * 60)
    
    nx_s, ny_s, nz_s = len(scan_xs), len(scan_ys), len(probe_heights)
    log(f"Scan grid: {nx_s}x{ny_s}x{nz_s}")
    log(f"Density grid: {E_pauli.shape}")
    
    # Build scan positions
    XX, YY, ZZ = np.meshgrid(scan_xs, scan_ys, probe_heights, indexing='ij')
    flat_pos = np.stack([XX.ravel(), YY.ravel(), ZZ.ravel()], axis=1)
    grid_c = (flat_pos - origin) / step  # fractional coordinates on density grid
    
    # Interpolate Pauli, ES, and vdW forces from density grid to scan positions
    log("\n6a. Interpolating Pauli/ES/vdW forces to scan positions...")
    F_pauli_scan = np.zeros((flat_pos.shape[0], 3), dtype=np.float32)
    F_es_scan = np.zeros((flat_pos.shape[0], 3), dtype=np.float32)
    F_vdw_scan = np.zeros((flat_pos.shape[0], 3), dtype=np.float32)
    
    for a in range(3):
        F_pauli_scan[:, a] = -map_coordinates(grads_pauli[:, :, :, a], grid_c.T, order=1, mode='nearest')
        F_es_scan[:, a]    = -map_coordinates(grads_es[:, :, :, a],    grid_c.T, order=1, mode='nearest')
        F_vdw_scan[:, a]   = -map_coordinates(grads_vdw[:, :, :, a],   grid_c.T, order=1, mode='nearest')
    
    # Total force
    F_total_scan = F_pauli_scan + F_vdw_scan + F_es_scan
    
    # Reshape to scan grid
    F_pauli = F_pauli_scan.reshape(nx_s, ny_s, nz_s, 3)
    F_es = F_es_scan.reshape(nx_s, ny_s, nz_s, 3)
    F_vdw = F_vdw_scan.reshape(nx_s, ny_s, nz_s, 3)
    F_total = F_total_scan.reshape(nx_s, ny_s, nz_s, 3)
    
    # Extract Fz
    Fz = F_total[:, :, :, 2]
    Fz_pauli = F_pauli[:, :, :, 2]
    Fz_es = F_es[:, :, :, 2]
    Fz_vdw = F_vdw[:, :, :, 2]
    
    log(f"Fz_pauli: range=[{Fz_pauli.min():.4f},{Fz_pauli.max():.4f}] eV/Å")
    log(f"Fz_es:    range=[{Fz_es.min():.4f},{Fz_es.max():.4f}] eV/Å")
    log(f"Fz_vdw:   range=[{Fz_vdw.min():.4f},{Fz_vdw.max():.4f}] eV/Å")
    log(f"Fz_total: range=[{Fz.min():.4f},{Fz.max():.4f}] eV/Å")
    
    # ═══════════════════════════════════════════════════════════════════════════
    # POTENTIAL ENERGY MAPS (V = E_pauli + E_es + E_vdw on scan grid)
    # ═══════════════════════════════════════════════════════════════════════════
    log("\n6a1. Computing potential energy maps on scan grid...")
    
    # Interpolate Pauli, ES, and vdW energy fields from density grid to scan grid
    E_pauli_scan = map_coordinates(E_pauli, grid_c.T, order=1, mode='nearest').reshape(nx_s, ny_s, nz_s)
    E_es_scan = map_coordinates(E_es, grid_c.T, order=1, mode='nearest').reshape(nx_s, ny_s, nz_s)
    E_vdw_scan = map_coordinates(E_vdw, grid_c.T, order=1, mode='nearest').reshape(nx_s, ny_s, nz_s)
    
    V_total = E_pauli_scan + E_es_scan + E_vdw_scan
    
    log(f"V_pauli:  range=[{E_pauli_scan.min():.4f},{E_pauli_scan.max():.4f}] eV")
    log(f"V_es:     range=[{E_es_scan.min():.4f},{E_es_scan.max():.4f}] eV")
    log(f"V_vdw:    range=[{E_vdw_scan.min():.4f},{E_vdw_scan.max():.4f}] eV")
    log(f"V_total:  range=[{V_total.min():.4f},{V_total.max():.4f}] eV")
    
    # Save potential maps
    save_npy(step_dir, 'V_pauli.npy', E_pauli_scan)
    save_npy(step_dir, 'V_es.npy', E_es_scan)
    save_npy(step_dir, 'V_vdw.npy', E_vdw_scan)
    save_npy(step_dir, 'V_total.npy', V_total)
    
    # Plot potential decomposition at z=5.0Å (XY view)
    log("\n6a2. Plotting potential energy decomposition at z=5.0Å...")
    iz_5 = np.argmin(np.abs(probe_heights - 5.0))
    if iz_5 >= nz_s:
        iz_5 = nz_s // 2
    h5 = probe_heights[iz_5]
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    fig.suptitle(f"Potential Energy Decomposition at z={h5:.1f}Å  (A_pauli={A_pauli}, beta={beta_pauli})")
    
    titles = ['Pauli (repulsive)', 'Electrostatic', 'vdW (attractive)', 'TOTAL']
    fields = [E_pauli_scan, E_es_scan, E_vdw_scan, V_total]
    cmaps = ['magma', 'bwr', 'viridis', 'magma']
    
    # Scan extent for all XY plots in step6
    scan_extent = [float(scan_xs[0]), float(scan_xs[-1]), float(scan_ys[0]), float(scan_ys[-1])]
    
    for ax, field, title, cmap in zip(axes.flat, fields, titles, cmaps):
        vmin, vmax = field[:, :, iz_5].min(), field[:, :, iz_5].max()
        if title in ['Electrostatic', 'TOTAL']:
            vabs = max(abs(vmin), abs(vmax), 1e-12)
            vmin, vmax = -vabs, vabs
            cmap = 'bwr'
        im = ax.imshow(field[:, :, iz_5].T, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax,
                       extent=scan_extent, aspect='equal')
        ax.set_title(f"{title}\nrange=[{vmin:.3e},{vmax:.3e}] eV")
        ax.set_xlabel("x [Å]")
        ax.set_ylabel("y [Å]")
        plt.colorbar(im, ax=ax, shrink=0.8)
    
    plt.tight_layout()
    plt.savefig(os.path.join(step_dir, "step6_V_decomposition_xy_z5.png"), dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved step6_V_decomposition_xy_z5.png")
    
    # Plot potential maps at several heights - per-image adaptive ranges
    log("\n6a3. Plotting V_total at several heights...")
    heights_to_plot = [probe_heights[0], probe_heights[nz_s//4], probe_heights[nz_s//2], probe_heights[3*nz_s//4], probe_heights[-1]]
    n_plots = len(heights_to_plot)
    fig, axes = plt.subplots(1, n_plots, figsize=(4*n_plots, 4))
    fig.suptitle("Total Potential V_total [eV] - per-image adaptive range")
    for ax, h in zip(axes.flat if n_plots > 1 else [axes], heights_to_plot):
        iz = np.argmin(np.abs(probe_heights - h))
        vmin, vmax = V_total[:, :, iz].min(), V_total[:, :, iz].max()
        # Use magma for all-positive or bwr for bipolar
        if vmin < -1e-6:
            vabs = max(abs(vmin), abs(vmax), 1e-12)
            im = ax.imshow(V_total[:, :, iz].T, origin='lower', cmap='bwr', vmin=-vabs, vmax=vabs,
                           extent=scan_extent, aspect='equal')
        else:
            im = ax.imshow(V_total[:, :, iz].T, origin='lower', cmap='magma', vmin=vmin, vmax=vmax,
                           extent=scan_extent, aspect='equal')
        ax.set_title(f"z={h:.1f}Å\n[{vmin:.2e},{vmax:.2e}]")
        ax.set_xlabel("x [Å]")
        ax.set_ylabel("y [Å]")
        plt.colorbar(im, ax=ax, shrink=0.8)
    plt.tight_layout()
    plt.savefig(os.path.join(step_dir, "step6_V_total_heights.png"), dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved step6_V_total_heights.png")
    
    # PP relaxation using FIRE-like algorithm from AFM.py
    log("\n6b. Probe particle relaxation...")
    try:
        from pyBall.OCL.AFM import pp_relax_2d
        
        # Build force_func that interpolates from scan grid
        from scipy.ndimage import map_coordinates
        def force_func(positions):
            """Interpolate forces at arbitrary positions from scan-grid force field."""
            positions = np.asarray(positions, dtype=np.float32)
            # Map positions to scan grid fractional coordinates
            ix = (positions[:, 0] - scan_xs[0]) / (scan_xs[-1] - scan_xs[0]) * (nx_s - 1)
            iy = (positions[:, 1] - scan_ys[0]) / (scan_ys[-1] - scan_ys[0]) * (ny_s - 1)
            # For z, find nearest height slice
            iz = np.argmin(np.abs(positions[:, 2:3] - probe_heights), axis=1)
            iz = np.clip(iz, 0, nz_s - 1)
            
            # Use map_coordinates for x,y interpolation at each z slice
            forces = np.zeros((len(positions), 3), dtype=np.float32)
            for a in range(3):
                # Interpolate each force component across x,y for all z slices
                for i, (xi, yi, zi) in enumerate(zip(ix, iy, iz)):
                    # Simple bilinear interpolation at nearest z
                    forces[i, a] = map_coordinates(F_total[:, :, :, a], 
                                                   [[xi], [yi], [zi]], 
                                                   order=1, mode='nearest')[0]
            return forces
        
        FEs_relax = pp_relax_2d(force_func, scan_xs, scan_ys, probe_heights, 
                                mol_z=0.0, K_LAT=0.03, N_RELAX=30, step=abs(scan_xs[1]-scan_xs[0]))
        Fz_relax = FEs_relax[:, :, :, 2]
        log(f"Fz_relax: range=[{Fz_relax.min():.4f},{Fz_relax.max():.4f}] eV/Å")
    except Exception as e:
        log(f"WARNING: PP relaxation failed ({e}), using raw forces")
        Fz_relax = Fz.copy()
    
    # Frequency shift df = -dFz/dz
    if nz_s > 1:
        dz = abs(probe_heights[1] - probe_heights[0])
        df = compute_df(Fz_relax, dz)
        log(f"df: range=[{df.min():.4f},{df.max():.4f}]")
    else:
        df = np.zeros_like(Fz_relax)
        log("WARNING: Only 1 z-layer, df=0")
    
    # Save outputs
    save_npy(step_dir, 'Fz_raw.npy', Fz)
    save_npy(step_dir, 'Fz_relax.npy', Fz_relax)
    save_npy(step_dir, 'df.npy', df)
    save_npy(step_dir, 'F_pauli.npy', F_pauli)
    save_npy(step_dir, 'F_es.npy', F_es)
    save_npy(step_dir, 'F_vdw.npy', F_vdw)
    
    # Plot component traces at center
    log("\n6c. Plotting component traces...")
    ix, iy = nx_s // 2, ny_s // 2
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle("Force Components at Center")
    
    ax = axes[0]
    ax.plot(probe_heights, Fz_pauli[ix, iy, :], 'purple', label='Pauli', lw=1.5)
    ax.plot(probe_heights, Fz_es[ix, iy, :], 'seagreen', label='ES', lw=1.5)
    ax.plot(probe_heights, Fz_vdw[ix, iy, :], 'royalblue', label='vdW', lw=1.5)
    ax.plot(probe_heights, Fz[ix, iy, :], 'k--', label='Total', lw=2.0)
    ax.axhline(0, color='k', lw=0.5)
    ax.set_xlabel('Height [Å]'); ax.set_ylabel('Fz [eV/Å]')
    ax.set_title('Full range'); ax.legend(); ax.invert_xaxis()
    
    ax = axes[1]
    ax.plot(probe_heights, Fz_pauli[ix, iy, :], 'purple', label='Pauli', lw=1.5)
    ax.plot(probe_heights, Fz_es[ix, iy, :], 'seagreen', label='ES', lw=1.5)
    ax.plot(probe_heights, Fz_vdw[ix, iy, :], 'royalblue', label='vdW', lw=1.5)
    ax.plot(probe_heights, Fz[ix, iy, :], 'k--', label='Total', lw=2.0)
    ax.axhline(0, color='k', lw=0.5)
    ax.set_xlabel('Height [Å]'); ax.set_ylabel('Fz [eV/Å]')
    ax.set_title('Zoom ±0.5'); ax.set_ylim(-0.5, 0.5); ax.legend(); ax.invert_xaxis()
    
    plt.tight_layout()
    plt.savefig(os.path.join(step_dir, "step6_component_traces.png"), dpi=120, bbox_inches='tight')
    plt.close()
    print("  Saved step6_component_traces.png")
    
    # Side-by-side FDBM vs LJ comparison plot at multiple heights
    if lj_forces is not None:
        log("\n6d. Plotting FDBM vs LJ side-by-side comparison...")
        Fz_lj_rep, Fz_lj_att, Fz_lj_es, Fz_lj_total = lj_forces
        
        fdmb_fields = [Fz_pauli, Fz_es, Fz_vdw, Fz]
        lj_fields = [Fz_lj_rep, Fz_lj_es, Fz_lj_att, Fz_lj_total]
        titles = ['Pauli / Repulsive', 'Electrostatic', 'London / Attractive', 'TOTAL']
        cmaps_fdmb = ['magma', 'bwr', 'viridis', 'magma']
        cmaps_lj = ['magma', 'bwr', 'viridis', 'magma']
        
        # Plot at z=5.0Å (original) and additional heights 3.0, 3.5, 4.0Å
        comp_heights = [5.0, 4.0, 3.5, 3.0]
        for h_comp in comp_heights:
            iz_comp = np.argmin(np.abs(probe_heights - h_comp))
            if iz_comp >= nz_s:
                iz_comp = nz_s // 2
            h_actual = probe_heights[iz_comp]
            
            fig, axes = plt.subplots(2, 4, figsize=(18, 9))
            fig.suptitle(f"FDBM (A={A_pauli},beta={beta_pauli}) vs LJ  z={h_actual:.1f}Å")
            
            # Per-image color scale (not common per column)
            for col, (f_f, f_l, title, cmap) in enumerate(zip(fdmb_fields, lj_fields, titles, cmaps_fdmb)):
                # FDBM row: auto-scale to its own range
                vmin_f = f_f[:, :, iz_comp].min(); vmax_f = f_f[:, :, iz_comp].max()
                if title == 'Electrostatic':
                    vabs_f = max(abs(vmin_f), abs(vmax_f), 1e-12)
                    vmin_f, vmax_f = -vabs_f, vabs_f
                im = axes[0, col].imshow(f_f[:, :, iz_comp].T, origin='lower', cmap='bwr' if title == 'Electrostatic' else cmap,
                                         vmin=vmin_f, vmax=vmax_f, extent=scan_extent, aspect='equal')
                axes[0, col].set_title(f"FDBM {title}\n[{vmin_f:.2e},{vmax_f:.2e}]")
                axes[0, col].set_xlabel("x [Å]"); axes[0, col].set_ylabel("y [Å]")
                plt.colorbar(im, ax=axes[0, col], shrink=0.7)
                
                # LJ row: auto-scale to its own range
                vmin_l = f_l[:, :, iz_comp].min(); vmax_l = f_l[:, :, iz_comp].max()
                if title == 'Electrostatic':
                    vabs_l = max(abs(vmin_l), abs(vmax_l), 1e-12)
                    vmin_l, vmax_l = -vabs_l, vabs_l
                im = axes[1, col].imshow(f_l[:, :, iz_comp].T, origin='lower', cmap='bwr' if title == 'Electrostatic' else cmap,
                                         vmin=vmin_l, vmax=vmax_l, extent=scan_extent, aspect='equal')
                axes[1, col].set_title(f"LJ {title}\n[{vmin_l:.2e},{vmax_l:.2e}]")
                axes[1, col].set_xlabel("x [Å]"); axes[1, col].set_ylabel("y [Å]")
                plt.colorbar(im, ax=axes[1, col], shrink=0.7)
            
            plt.tight_layout()
            fname = f"step6_fdbm_vs_lj_decomposition_z{h_actual:.1f}.png"
            plt.savefig(os.path.join(step_dir, fname), dpi=150, bbox_inches='tight')
            plt.close()
            print(f"  Saved {fname}")
        
        # Side-by-side traces comparison
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        fig.suptitle("FDBM vs LJ Force Traces at Center")
        
        ax = axes[0]
        ax.plot(probe_heights, Fz_pauli[ix, iy, :], 'purple', label='FDBM Pauli', lw=1.5)
        ax.plot(probe_heights, Fz_es[ix, iy, :], 'seagreen', label='FDBM ES', lw=1.5)
        ax.plot(probe_heights, Fz_vdw[ix, iy, :], 'royalblue', label='FDBM London', lw=1.5)
        ax.plot(probe_heights, Fz[ix, iy, :], 'k--', label='FDBM Total', lw=2.0)
        ax.axhline(0, color='k', lw=0.5)
        ax.set_xlabel('Height [Å]'); ax.set_ylabel('Fz [eV/Å]')
        ax.set_title('FDBM'); ax.legend(); ax.invert_xaxis()
        
        ax = axes[1]
        ax.plot(probe_heights, Fz_lj_rep[ix, iy, :], 'purple', label='LJ Repulsive', lw=1.5)
        ax.plot(probe_heights, Fz_lj_es[ix, iy, :], 'seagreen', label='LJ ES', lw=1.5)
        ax.plot(probe_heights, Fz_lj_att[ix, iy, :], 'royalblue', label='LJ Attractive', lw=1.5)
        ax.plot(probe_heights, Fz_lj_total[ix, iy, :], 'k--', label='LJ Total', lw=2.0)
        ax.axhline(0, color='k', lw=0.5)
        ax.set_xlabel('Height [Å]'); ax.set_ylabel('Fz [eV/Å]')
        ax.set_title('LJ'); ax.legend(); ax.invert_xaxis()
        
        plt.tight_layout()
        plt.savefig(os.path.join(step_dir, "step6_fdbm_vs_lj_traces.png"), dpi=120, bbox_inches='tight')
        plt.close()
        print("  Saved step6_fdbm_vs_lj_traces.png")
    
    # df maps at different heights
    if nz_s > 1:
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        fig.suptitle("Frequency Shift df = -dFz/dz")
        for i, iz in enumerate([0, nz_s//4, nz_s//2, 3*nz_s//4, nz_s-2, nz_s-1]):
            ax = axes.flat[i]
            vabs = max(abs(float(df[:, :, iz].min())), abs(float(df[:, :, iz].max())), 1e-6)
            norm = TwoSlopeNorm(vmin=-vabs, vcenter=0, vmax=vabs)
            im = ax.imshow(df[:, :, iz].T, origin='lower', cmap='bwr', norm=norm)
            ax.set_title(f"h={probe_heights[iz]:.2f}Å")
            plt.colorbar(im, ax=ax, shrink=0.8)
        plt.tight_layout()
        plt.savefig(os.path.join(step_dir, "step6_df_maps.png"), dpi=120, bbox_inches='tight')
        plt.close()
        print("  Saved step6_df_maps.png")
    
    log("\nStep 6 COMPLETE.")
    log("=" * 60)
    # Sanity checks
    check_field_sanity(Fz, "Fz", min_abs=1e-8, max_abs=1e6)
    check_field_sanity(Fz_relax, "Fz_relax", min_abs=1e-8, max_abs=1e6)
    
    log("ALL STEPS COMPLETE.")
    log("=" * 60)
    
    return Fz, Fz_relax, df

# ═══════════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    ap = argparse.ArgumentParser(description="PyOpenCL FDBM AFM step-by-step")
    ap.add_argument('--step', type=float, default=0.10, help='Density grid step [Å] (default 0.1)')
    ap.add_argument('--nscf', type=int, default=200, help='SCF iterations')
    ap.add_argument('--A_pauli', type=float, default=3.0, help='Pauli amplitude')
    ap.add_argument('--beta_pauli', type=float, default=1.0, help='Exponent applied to sample density before Pauli convolution (default 1.0, try 1.2)')
    ap.add_argument('--C6_CO', type=float, default=30.0, help='CO-tip C6')
    ap.add_argument('--q_CO', type=float, default=-0.05, help='CO-tip charge')
    ap.add_argument('--sigma_tip', type=float, default=0.7, help='Tip Gaussian width')
    ap.add_argument('--tip_model', type=str, default='gaussian', choices=['gaussian', 'co'],
                   help='Tip model: gaussian or true CO density from Fireball')
    ap.add_argument('--co_nscf', type=int, default=50, help='SCF iterations for CO tip density')
    ap.add_argument('--z_start', type=float, default=6.0, help='Scan start height')
    ap.add_argument('--z_end', type=float, default=0.4, help='Scan end height')
    ap.add_argument('--dz', type=float, default=0.2, help='Scan height step')
    ap.add_argument('--nxy', type=int, nargs=2, default=None, help='Scan grid (nx,ny). If not set, computed from --scan_step')
    ap.add_argument('--scan_step', type=float, default=0.1, help='Target scan pixel size in Angstrom (default 0.1A)')
    ap.add_argument('--stop_after', type=int, default=0, help='Stop after step N (1-5)')
    ap.add_argument('--verbosity', type=int, default=1)
    args = ap.parse_args()
    
    setup_debug_dirs()
    
    # Load molecule
    atomTypes, atomPos = load_xyz(_XYZ_PATH)
    print(f"Loaded pentacene: {len(atomTypes)} atoms")
    print(f"  z-range: [{atomPos[:,2].min():.2f}, {atomPos[:,2].max():.2f}]")
    
    # Scan grid
    nz = int(round((args.z_start - args.z_end) / args.dz)) + 1
    probe_heights = args.z_start + np.arange(nz) * (-args.dz)
    # Scan bounding box: molecule extent + 3Å margin on each side
    x_span = atomPos[:,0].max() - atomPos[:,0].min() + 6.0  # +3Å each side
    y_span = atomPos[:,1].max() - atomPos[:,1].min() + 6.0
    if args.nxy is not None:
        nx_scan, ny_scan = args.nxy[0], args.nxy[1]
    else:
        nx_scan = max(20, int(round(x_span / args.scan_step)))
        ny_scan = max(20, int(round(y_span / args.scan_step)))
    scan_xs = np.linspace(atomPos[:,0].min() - 3.0, atomPos[:,0].max() + 3.0, nx_scan)
    scan_ys = np.linspace(atomPos[:,1].min() - 3.0, atomPos[:,1].max() + 3.0, ny_scan)
    dx_scan = scan_xs[1] - scan_xs[0] if len(scan_xs) > 1 else args.scan_step
    dy_scan = scan_ys[1] - scan_ys[0] if len(scan_ys) > 1 else args.scan_step
    print(f"Scan grid: {nx_scan}x{ny_scan}x{nz}  dx={dx_scan:.3f}A dy={dy_scan:.3f}A  h={args.z_start}..{args.z_end}")
    
    # Step 1: Density projection
    rho_grid, rho_na_grid, rho_diff, origin, ngrid, grid_spec, q_mulliken, atomPos, atomTypes = \
        step1_density_projection(atomTypes, atomPos, step=args.step, nscf=args.nscf,
                                verbosity=args.verbosity, stop_after=(args.stop_after==1))
    
    # Step 2: Electrostatics
    V_ES = step2_electrostatics(rho_diff, args.step, origin, ngrid, stop_after=(args.stop_after==2))
    
    # Compute tip density if using true CO model
    tip_density_total = None
    tip_density_delta = None
    if args.tip_model == 'co':
        print(f"\n{'='*60}")
        print("COMPUTING TRUE CO TIP DENSITY")
        print(f"{'='*60}")
        tip_density_total, tip_density_delta, _ = get_co_tip_density(
            grid_spec, fdata_dir=_FDATA_HCNO, fdata_basis=_FDATA_BASIS,
            step=args.step, nscf=args.co_nscf, verbosity=args.verbosity)
        # Save CO diagnostic plots
        co_dir = os.path.join(_DEBUG_DIR, 'co_tip')
        os.makedirs(co_dir, exist_ok=True)
        
        # CO projected at grid center; find center z for diagnostic slice
        co_z_center = float(origin[2]) + (int(ngrid[2]) // 2) * args.step
        # CO total density slices
        plot_field_slices(tip_density_total, "CO Total Density (rho_SCF) [e/Å³]",
                          "co_rho_total_slices.png", co_dir, z_slice=co_z_center, origin=origin, step=args.step)
        # CO delta density slices (symmetric)
        plot_field_slices(tip_density_delta, "CO Delta Density (rho_SCF - rho_NA) [e/Å³]",
                          "co_rho_delta_slices.png", co_dir, z_slice=co_z_center, origin=origin, step=args.step,
                          sym=True, cmap='bwr')
        # CO delta density xz cut through axis (y=0, i.e. center)
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        fig.suptitle("CO Delta Density xz View (through molecular axis)")
        ny_d = tip_density_delta.shape[1]
        iy_axis = ny_d // 2
        vabs = max(abs(float(tip_density_delta.min())), abs(float(tip_density_delta.max())), 1e-12)
        nx_d, nz_d = tip_density_delta.shape[0], tip_density_delta.shape[2]
        # Extent for xz: x-axis and z-axis in Å
        ext_xz = [float(origin[0]), float(origin[0]) + (nx_d-1)*args.step,
                    float(origin[2]), float(origin[2]) + (nz_d-1)*args.step]
        im0 = axes[0].imshow(tip_density_delta[:, iy_axis, :].T, origin='lower', cmap='bwr',
                             vmin=-vabs, vmax=vabs, extent=ext_xz, aspect='equal')
        axes[0].set_title(f"CO delta-rho xz (iy={iy_axis})\nrange=[{-vabs:.3e},{vabs:.3e}]")
        axes[0].set_xlabel("x [Å]")
        axes[0].set_ylabel("z [Å]")
        plt.colorbar(im0, ax=axes[0])
        # CO total density xz
        im1 = axes[1].imshow(tip_density_total[:, iy_axis, :].T, origin='lower', cmap='magma',
                             extent=ext_xz, aspect='equal')
        axes[1].set_title(f"CO total rho xz (iy={iy_axis})")
        axes[1].set_xlabel("x [Å]")
        axes[1].set_ylabel("z [Å]")
        plt.colorbar(im1, ax=axes[1])
        plt.tight_layout()
        plt.savefig(os.path.join(co_dir, "co_rho_xz_axis.png"), dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  Saved co_rho_xz_axis.png")
    
    # Step 3: Pauli
    E_pauli, grads_pauli = step3_pauli(rho_grid, origin, args.step, scan_xs, scan_ys, probe_heights,
                                       A_pauli=args.A_pauli, sigma_tip=args.sigma_tip,
                                       tip_density=tip_density_total,
                                       beta_pauli=args.beta_pauli,
                                       stop_after=(args.stop_after==3))
    
    # Step 4: Electrostatics convolution
    E_es, grads_es = step4_electrostatics_conv(V_ES, origin, args.step, scan_xs, scan_ys, probe_heights,
                                               q_CO=args.q_CO, sigma_tip=args.sigma_tip,
                                               tip_charge_density=tip_density_delta,
                                               stop_after=(args.stop_after==4))
    
    # Step 5: Dispersion (on density grid, like Pauli/ES)
    E_vdw, grads_vdw = step5_dispersion(atomPos, atomTypes, origin, args.step, ngrid,
                                       C6_CO=args.C6_CO, stop_after=(args.stop_after==5))
    
    # Step 5b: LJ comparison (for side-by-side validation)
    lj_forces = compute_lj_forces(scan_xs, scan_ys, probe_heights, atomPos, atomTypes, q_mulliken,
                                  q_CO=args.q_CO, sigma_tip=3.0, eps_tip=0.003)
    
    # Step 6: Composed
    Fz, Fz_relax, df = step6_composed(E_pauli, grads_pauli, E_es, grads_es, E_vdw, grads_vdw,
                                      origin, args.step, scan_xs, scan_ys, probe_heights, atomPos,
                                      A_pauli=args.A_pauli, beta_pauli=args.beta_pauli, q_CO=args.q_CO, lj_forces=lj_forces)
    
    print("\n=== ALL DONE ===")
    print(f"Debug outputs saved to: {_DEBUG_DIR}")

if __name__ == '__main__':
    main()
