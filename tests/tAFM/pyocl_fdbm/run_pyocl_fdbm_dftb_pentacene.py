#!/usr/bin/env python3
"""
run_pyocl_fdbm_dftb_pentacene.py - Step-by-step FDBM AFM with DFTB+ backend

Refined version with proper normalization and process isolation for DFTBcore.
"""

import sys, os, argparse, json, time, shutil, subprocess
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from pathlib import Path
from multiprocessing import Process, Queue

# Add project root to path
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.realpath(os.path.join(_THIS_DIR, '..', '..', '..'))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from pyBall.DFTB.DFTBcore import DFTBcore
from pyBall.DFTB.DFTBplusParser import (
    parse_basis_hsd_ang, parse_wfc_hsd, convert_wfc_to_species_list_ang,
    parse_detailed_xml_custom, parse_eigenvec_bin_custom, evec_to_kernel_coeffs,
    precompute_coeff_gather
)
from pyBall.DFTB.Grid_dftb import setup_gridprojector_from_dftb
from pyBall.OCL import AFM as afm
from pyBall.OCL import clUtils as clu
from pyBall import atomicUtils as au
from pyBall import dftb_utils as du

# ==============================================================================
# Constants & Config
# ==============================================================================

BOHR2ANG = 0.5291772109
COULOMB_CONST = 14.3996448915
ELEM_Z = {'H':1, 'C':6, 'N':7, 'O':8, 'P':15, 'S':16}

# Normalization factor for density in Angstrom grid
# If orbitals are normalized in Bohr, density integral in Angstrom is B^3.
# To get 1.0 in Angstrom, multiply density by B^-3.
B3_FACTOR = 1.0 / (BOHR2ANG**3)

def setup_debug_dirs(debug_dir):
    if os.path.exists(debug_dir):
        shutil.rmtree(debug_dir)
    for step in ['step1_density', 'step2_electrostatics', 'step3_pauli',
                 'step4_electrostatics_conv', 'step5_dispersion', 'step6_composed', 'co_tip']:
        os.makedirs(os.path.join(debug_dir, step), exist_ok=True)

def load_xyz_as_types_pos(fname):
    """Load XYZ and return (atomic_numbers, positions)."""
    pos, _, names, _, _ = au.load_xyz(fname)
    types = np.array([ELEM_Z.get(e, 6) for e in names], dtype=np.int32)
    pos = np.array(pos, dtype=np.float64)
    return types, pos

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
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle(title)
    nx, ny, nz = data.shape
    extent_xy = None; extent_xz = None
    if origin is not None and step is not None:
        extent_xy = [float(origin[0]), float(origin[0]) + (nx-1)*step,
                     float(origin[1]), float(origin[1]) + (ny-1)*step]
        extent_xz = [float(origin[0]), float(origin[0]) + (nx-1)*step,
                     float(origin[2]), float(origin[2]) + (nz-1)*step]
    iz = int(np.clip(np.round((z_slice - origin[2]) / step), 0, nz-1)) if origin is not None else nz // 2
    
    slice_xy = data[:, :, iz]
    vabs = max(abs(float(slice_xy.min())), abs(float(slice_xy.max())), 1e-12)
    norm = TwoSlopeNorm(vmin=-vabs, vcenter=0, vmax=vabs) if sym else None
    im0 = axes[0].imshow(slice_xy.T, origin='lower', cmap=cmap, norm=norm, extent=extent_xy, aspect='equal')
    axes[0].set_title(f"XY z={z_slice:.1f}Å (iz={iz})")
    fig.colorbar(im0, ax=axes[0])
    
    iy = ny // 2
    slice_xz = data[:, iy, :]
    vabs = max(abs(float(slice_xz.min())), abs(float(slice_xz.max())), 1e-12)
    norm = TwoSlopeNorm(vmin=-vabs, vcenter=0, vmax=vabs) if sym else None
    im1 = axes[1].imshow(slice_xz.T, origin='lower', cmap=cmap, norm=norm, extent=extent_xz, aspect='equal')
    axes[1].set_title(f"XZ y-center")
    fig.colorbar(im1, ax=axes[1])
    
    ix, iy_c = nx // 2, ny // 2
    zlabel = np.arange(nz) if zvals is None else zvals
    axes[2].plot(zlabel, data[ix, iy_c, :], 'k-', lw=1.0)
    axes[2].set_title(f"Line profile at center")
    axes[2].set_xlabel("z [Å]" if zvals is not None else "z index")
    axes[2].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(step_dir, fname), dpi=120, bbox_inches='tight'); plt.close()

# ==============================================================================
# DFTB Helpers (Isolated in Process)
# ==============================================================================

def scf_worker(atomTypes, atomPos, work_dir, queue, slako_prefix):
    """Function to be run in a separate process."""
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    inv_z = {v: k for k, v in ELEM_Z.items()}
    enames = [inv_z[z] for z in atomTypes]
    au.save_xyz(str(work_dir / "geom.xyz"), enames, atomPos)

    species_used = sorted(set(atomTypes.tolist()))
    max_ang = "\n    ".join([f'{inv_z[z]} = "p"' if z > 1 else f'{inv_z[z]} = "s"' for z in species_used])

    hsd_content = f"""
Geometry = xyzFormat {{
  <<< "geom.xyz"
}}
Hamiltonian = DFTB {{
  SCC = Yes
  SlaterKosterFiles = Type2FileNames {{
    Prefix = "{slako_prefix}"
    Separator = "-"
    Suffix = ".skf"
  }}
  MaxAngularMomentum {{
    {max_ang}
  }}
}}
Analysis {{
  WriteEigenvectors = Yes
}}
Options {{
  WriteDetailedXml = Yes
}}
"""
    with open(work_dir / "dftb_in.hsd", "w") as f:
        f.write(hsd_content)
    
    orig_cwd = os.getcwd()
    os.chdir(work_dir)
    try:
        dftb = DFTBcore()
        dftb.init("dftb_in.hsd")
        dftb.run_scf()
        geo = parse_detailed_xml_custom("detailed.xml")
        evecs = parse_eigenvec_bin_custom("eigenvec.bin", geo['nstates'], geo['norb'])
        # Can't easily pass the whole geo/evecs back via queue (too big/complex)
        # Instead, we save them to disk or just use them here.
        # But we need them for projection in the main process.
        # Actually, let's just return the necessary pieces.
        queue.put((geo, evecs))
    except Exception as e:
        queue.put(e)
    finally:
        os.chdir(orig_cwd)

def run_isolated_scf(atomTypes, atomPos, work_dir, slako_prefix):
    q = Queue()
    p = Process(target=scf_worker, args=(atomTypes, atomPos, work_dir, q, slako_prefix))
    p.start()
    res = q.get()
    p.join()
    if isinstance(res, Exception):
        raise res
    return res

def project_dftb_density(geo, evecs, projector, atoms_dict, grid_spec, basis):
    natoms    = geo['natoms']
    occs      = geo['occupations'][:, 0, 0]
    occ_idx   = np.where(occs > 1e-6)[0]
    print(f"[project_dftb_density] {len(occ_idx)} occupied states")

    t0 = time.time()
    src_idx, dst_idx = precompute_coeff_gather(natoms, geo['species_per_atom'], geo['species_names'], basis)
    proj_ctx = projector.prepare_orbital_projection(atoms_dict, grid_spec)
    print(f"[project_dftb_density] Setup: {time.time()-t0:.3f}s")

    nx, ny, nz = proj_ctx['nx'], proj_ctx['ny'], proj_ctx['nz']
    rho_grid   = np.zeros((nx, ny, nz), dtype=np.float32)
    coeffs_flat = np.zeros(natoms * 4, dtype=np.float32)

    t1 = time.time()
    for idx, i in enumerate(occ_idx):
        coeffs_flat[:] = 0.0
        coeffs_flat[dst_idx] = evecs[i][src_idx]
        psi = projector.project_orbital_prepped(coeffs_flat, proj_ctx)
        rho_grid += occs[i] * (psi ** 2)
    n = len(occ_idx)
    print(f"[project_dftb_density] Loop: {time.time()-t1:.3f}s  ({(time.time()-t1)/n*1e3:.1f} ms/orbital)")
    return rho_grid * B3_FACTOR

def project_neutral_density(geo, projector, atoms_dict, grid_spec, basis):
    natoms          = geo['natoms']
    species_per_atom = geo['species_per_atom']
    species_names   = geo['species_names']
    sp_by_name      = {sp['name']: sp for sp in basis}

    src_idx, dst_idx = precompute_coeff_gather(natoms, species_per_atom, species_names, basis)
    proj_ctx = projector.prepare_orbital_projection(atoms_dict, grid_spec)

    nx, ny, nz    = proj_ctx['nx'], proj_ctx['ny'], proj_ctx['nz']
    rho_na_grid   = np.zeros((nx, ny, nz), dtype=np.float32)
    coeffs_flat   = np.zeros(natoms * 4, dtype=np.float32)

    # Neutral atomic orbital occupations:  {Z: {l: f_per_orbital}}
    OCC_NA = {1: {0: 1.0}, 6: {0: 2.0, 1: 2/3}, 7: {0: 2.0, 1: 1.0}, 8: {0: 2.0, 1: 4/3}}

    t0 = time.time()
    # dst_idx[k] = flat GPU index for eigenvector component k
    # For neutral density: orbital k has evec = e_k (identity row), so:
    #   coeffs_flat[dst_idx[k]] = 1.0, rest zero
    orb_offset = 0
    for ia in range(natoms):
        Z  = int(atoms_dict['type'][ia])
        sp = sp_by_name[species_names[species_per_atom[ia]]]
        occ_by_l = OCC_NA.get(Z, {0: 2.0, 1: 2/3})
        for orb in sp['orbitals']:
            l  = orb['l']
            nm = 2 * l + 1
            f_na = occ_by_l.get(l, 0.0)
            if f_na > 0:
                for m in range(nm):
                    coeffs_flat[:] = 0.0
                    coeffs_flat[dst_idx[orb_offset + m]] = 1.0
                    psi = projector.project_orbital_prepped(coeffs_flat, proj_ctx)
                    rho_na_grid += f_na * (psi ** 2)
            orb_offset += nm
    print(f"[project_neutral_density] Loop: {time.time()-t0:.3f}s")
    return rho_na_grid * B3_FACTOR

# ==============================================================================
# Step 1: Density Projection
# ==============================================================================

def step1_density_projection(atomTypes, atomPos, basis, slako_prefix, debug_dir, step=0.15, margin=4.0, z_extra=6.0, stop_after=False):
    step_dir = os.path.join(debug_dir, 'step1_density')
    log = lambda msg: log_step(step_dir, msg)
    log("="*60); log("STEP 1: Density Projection to Grid (DFTB+)"); log("="*60)

    log("\n1a. Running DFTB+ SCF (isolated process)...")
    t0 = time.time()
    geo, evecs = run_isolated_scf(atomTypes, atomPos, os.path.join(step_dir, "mol_work"), slako_prefix)
    print(f"[step1] DFTB SCF completed in {time.time()-t0:.2f}s")
    
    log("\n1b. Setting up density grid...")
    t1 = time.time()
    pos_ang = atomPos
    grid_spec, origin, ngrid = afm.setup_density_grid(pos_ang, step=step, margin=margin, z_extra=z_extra)
    log(f"Grid: origin={origin.round(2)} ngrid={ngrid}")
    print(f"[step1] Grid setup completed in {time.time()-t1:.2f}s")
    
    log("\n1c. Projecting densities...")
    t2 = time.time()
    dftb_data = {'coords_bohr': geo['coords_bohr'], 'species_per_atom': geo['species_per_atom'], 'species_names': geo['species_names']}
    print("[step1] Setting up projector...")
    t2a = time.time()
    projector, atoms_dict = setup_gridprojector_from_dftb(dftb_data, basis, verbosity=0)
    print(f"[step1] Projector setup completed in {time.time()-t2a:.2f}s")
    
    print("[step1] Projecting total density...")
    t2b = time.time()
    rho_grid = project_dftb_density(geo, evecs, projector, atoms_dict, grid_spec, basis)
    print(f"[step1] Total density projection completed in {time.time()-t2b:.2f}s")
    
    print("[step1] Projecting neutral density...")
    t2c = time.time()
    rho_na_grid = project_neutral_density(geo, projector, atoms_dict, grid_spec, basis)
    print(f"[step1] Neutral density projection completed in {time.time()-t2c:.2f}s")
    
    rho_diff = (rho_grid - rho_na_grid).astype(np.float32)
    print(f"[step1] Step 1 total time: {time.time()-t0:.2f}s")
    
    log(f"  rho_grid: integral={rho_grid.sum()*step**3:.4f} e")
    log(f"  rho_NA:   integral={rho_na_grid.sum()*step**3:.4f} e")
    log(f"  delta_rho: integral={rho_diff.sum()*step**3:.4f} e")
    
    save_npy(step_dir, 'rho_grid.npy', rho_grid)
    save_npy(step_dir, 'rho_na_grid.npy', rho_na_grid)
    save_npy(step_dir, 'rho_diff.npy', rho_diff)
    
    plot_field_slices(rho_grid, "SCF Density [e/Å³]", "step1_rho_slices.png", step_dir, origin=origin, step=step)
    plot_field_slices(rho_na_grid, "Neutral Atom Density [e/Å³]", "step1_rhoNA_slices.png", step_dir, origin=origin, step=step)
    plot_field_slices(rho_diff, "Delta Density [e/Å³]", "step1_rhoDiff_slices.png", step_dir, origin=origin, step=step, sym=True, cmap='bwr')
    
    if stop_after: sys.exit(0)
    return rho_grid, rho_na_grid, rho_diff, origin, ngrid, grid_spec, atomPos, atomTypes

# ==============================================================================
# Pipeline Steps (mostly identical to run_pyocl_fdbm.py)
# ==============================================================================

def step2_electrostatics(rho_diff, step, origin, ngrid, debug_dir):
    step_dir = os.path.join(debug_dir, 'step2_electrostatics')
    log = lambda msg: log_step(step_dir, msg)
    log("="*60); log("STEP 2: Electrostatics via Poisson"); log("="*60)
    V_ES = afm.fft_poisson(rho_diff, step)
    log(f"V_ES range: [{V_ES.min():.3f}, {V_ES.max():.3f}] eV")
    save_npy(step_dir, 'V_ES.npy', V_ES)
    plot_field_slices(V_ES, "Electrostatic Potential [eV]", "step2_VES_slices.png", step_dir, origin=origin, step=step, sym=True, cmap='bwr')
    return V_ES

def step3_pauli(rho_grid, origin, step, rho_tip_total, debug_dir, A_pauli=16.0, beta_pauli=1.0):
    step_dir = os.path.join(debug_dir, 'step3_pauli')
    log = lambda msg: log_step(step_dir, msg)
    log("="*60); log("STEP 3: Pauli Repulsion Convolution"); log("="*60)
    E_pauli_field, grads_E_Pauli = afm.compute_pauli_field(rho_grid, rho_tip_total, step, A_pauli=A_pauli, beta_pauli=beta_pauli)
    save_npy(step_dir, 'E_Pauli_field.npy', E_pauli_field)
    # Also save raw overlap for debugging
    dV = step**3
    nx_t, ny_t, nz_t = rho_tip_total.shape
    tip_kernel = np.roll(np.roll(np.roll(rho_tip_total[::-1,::-1,::-1], -(nx_t//2), axis=0), -(ny_t//2), axis=1), -(nz_t//2), axis=2)
    overlap_raw = dV * np.real(np.fft.ifftn(np.fft.fftn(rho_grid) * np.fft.fftn(tip_kernel))).astype(np.float32)
    save_npy(step_dir, 'overlap_raw.npy', overlap_raw)
    plot_field_slices(E_pauli_field, f"Pauli Energy [eV] (A={A_pauli:.1f}, b={beta_pauli:.3f})", "step3_Epauli_slices.png", step_dir, origin=origin, step=step)
    log(f"Pauli field range: [{E_pauli_field.min():.4f}, {E_pauli_field.max():.4f}] eV")
    return E_pauli_field, grads_E_Pauli

def step4_electrostatics_conv(V_ES, origin, step, rho_tip_delta, debug_dir):
    step_dir = os.path.join(debug_dir, 'step4_electrostatics_conv')
    log = lambda msg: log_step(step_dir, msg)
    log("="*60); log("STEP 4: Electrostatics Convolution"); log("="*60)
    E_es_field, grads_E_ES = afm.compute_es_conv_field(V_ES, rho_tip_delta, step)
    save_npy(step_dir, 'E_ES_field.npy', E_es_field)
    plot_field_slices(E_es_field, "ES Energy [eV]", "step4_EES_slices.png", step_dir, origin=origin, step=step, sym=True, cmap='bwr')
    return E_es_field, grads_E_ES

def step5_dispersion(atomPos, atomTypes, origin, step, ngrid, debug_dir, C6_CO=30.0):
    step_dir = os.path.join(debug_dir, 'step5_dispersion')
    log = lambda msg: log_step(step_dir, msg)
    log("="*60); log("STEP 5: Dispersion (C6/r^6)"); log("="*60)
    nx, ny, nz = [int(i) for i in ngrid[:3]]
    xs = origin[0] + np.arange(nx)*step; ys = origin[1] + np.arange(ny)*step; zs = origin[2] + np.arange(nz)*step
    XX, YY, ZZ = np.meshgrid(xs, ys, zs, indexing='ij')
    E_vdw = np.zeros((nx, ny, nz), dtype=np.float32)
    C6_atom = np.array([15.0 if z == 6 else 1.0 for z in atomTypes])
    RA2 = 1.5**2
    for ia in range(len(atomPos)):
        r2 = (XX-atomPos[ia,0])**2 + (YY-atomPos[ia,1])**2 + (ZZ-atomPos[ia,2])**2
        E_vdw -= np.sqrt(C6_atom[ia]*C6_CO) / (r2 + RA2)**3
    grads_E_vdw = np.stack([np.gradient(E_vdw, step, axis=i) for i in range(3)], axis=-1)
    save_npy(step_dir, 'E_vdw_field.npy', E_vdw)
    plot_field_slices(E_vdw, "Dispersion Energy [eV]", "step5_Evdw_slices.png", step_dir, origin=origin, step=step)
    return E_vdw, grads_E_vdw

def step6_composed(grads_pauli, grads_es, grads_vdw, origin, step, atomPos, scan_xs, scan_ys, heights, debug_dir, K_LAT=0.5):
    step_dir = os.path.join(debug_dir, 'step6_composed')
    log = lambda msg: log_step(step_dir, msg)
    log("="*60); log("STEP 6: Final Composed & Relaxation"); log("="*60)

    # F = -grad E
    F_total = -(grads_pauli + grads_es + grads_vdw)

    nx_s, ny_s = len(scan_xs), len(scan_ys)
    nz_s = len(heights)

    from scipy.ndimage import map_coordinates
    def force_func(positions):
        """Interpolate forces at arbitrary positions from scan-grid force field."""
        ix = (positions[:, 0] - origin[0]) / step
        iy = (positions[:, 1] - origin[1]) / step
        iz = (positions[:, 2] - origin[2]) / step
        coords = np.vstack([ix, iy, iz])
        fx = map_coordinates(F_total[..., 0], coords, order=1)
        fy = map_coordinates(F_total[..., 1], coords, order=1)
        fz = map_coordinates(F_total[..., 2], coords, order=1)
        return np.stack([fx, fy, fz], axis=-1)

    log("\n6b. Probe particle relaxation...")
    FEs_relax = afm.pp_relax_2d(force_func, scan_xs, scan_ys, heights, mol_z=atomPos[:,2].max(), K_LAT=K_LAT, N_RELAX=50, step=step)
    Fz_relax = FEs_relax[:,:,:,2]

    df = afm.compute_df(Fz_relax, heights[1]-heights[0])
    save_npy(step_dir, 'df.npy', df)
    afm.save_afm_images(df, scan_xs, scan_ys, heights, step_dir, prefix='df')

    log("Step 6 COMPLETE.")
    return df

# ==============================================================================
# Main
# ==============================================================================

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--step', type=float, default=0.15)
    parser.add_argument('--xyz', type=str, default='pentacene.xyz', help='Molecule XYZ file')
    parser.add_argument('--output_dir', type=str, default=None, help='Custom output directory name')
    parser.add_argument('--basis', type=str, default='mio-1-1', help='Basis set: mio-1-1 or 3ob-3-1')
    parser.add_argument('--slako_prefix', type=str, default=None, help='SLAKO files directory prefix')
    parser.add_argument('--basis_hsd', type=str, default=None, help='Basis HSD file path')
    parser.add_argument('--A_pauli', type=float, default=None, help='Pauli prefactor A (default: fitted per basis)')
    parser.add_argument('--beta_pauli', type=float, default=None, help='Pauli exponent beta (default: fitted per basis)')
    args = parser.parse_args()

    debug_dir = os.path.join(_THIS_DIR, args.output_dir) if args.output_dir else os.path.join(_THIS_DIR, f'debug_dftb_{args.basis.replace("-","_")}')

    # Resolve slako and basis_hsd defaults
    slako_prefix = args.slako_prefix if args.slako_prefix else du.SK_PATHS.get(args.basis, du.SK_PATHS.get('mio-1-1', ''))
    basis_hsd = args.basis_hsd if args.basis_hsd else du.WFC_HSD_PATHS.get(args.basis, du.WFC_HSD_PATHS.get('mio-1-1', ''))
    if not slako_prefix or not basis_hsd:
        raise ValueError(f"Cannot resolve paths for basis={args.basis}. Provide --slako_prefix and --basis_hsd.")

    A_pauli = args.A_pauli if args.A_pauli is not None else afm.PAULI_FITTED_DEFAULTS[args.basis]['A']
    beta_pauli = args.beta_pauli if args.beta_pauli is not None else afm.PAULI_FITTED_DEFAULTS[args.basis]['beta']

    print(f"Using basis: {args.basis}")
    print(f"SLAKO_PREFIX: {slako_prefix}")
    print(f"BASIS_HSD: {basis_hsd}")
    print(f"Pauli params: A={A_pauli:.2f}, beta={beta_pauli:.4f}")

    setup_debug_dirs(debug_dir)
    atomTypes, atomPos = load_xyz_as_types_pos(os.path.join(_THIS_DIR, args.xyz))

    # Load Basis once
    if args.basis == '3ob-3-1':
        basis_data = parse_wfc_hsd(basis_hsd)
        basis = convert_wfc_to_species_list_ang(basis_data, resolution_bohr=0.04)
    else:
        basis = parse_basis_hsd_ang(basis_hsd)
    print(f"Loaded basis for species: {[sp['name'] for sp in basis]}")

    # Step 1-2
    rho_grid, rho_na_grid, rho_diff, origin, ngrid, grid_spec, atomPos, atomTypes = \
        step1_density_projection(atomTypes, atomPos, basis, slako_prefix, debug_dir, step=args.step)
    V_ES = step2_electrostatics(rho_diff, args.step, origin, ngrid, debug_dir)

    # CO Tip Calculation (isolated process)
    log_step(os.path.join(debug_dir, 'co_tip'), "Computing CO Tip density...")
    geo_co, evec_co = run_isolated_scf(np.array([8, 6]), np.array([[0,0,0], [0,0,1.13]]), os.path.join(debug_dir, 'co_tip', 'co_work'), slako_prefix)
    dftb_co = {'coords_bohr': geo_co['coords_bohr'], 'species_per_atom': geo_co['species_per_atom'], 'species_names': geo_co['species_names']}
    projector_co, atoms_dict_co = setup_gridprojector_from_dftb(dftb_co, basis, verbosity=0)

    grid_center = origin + 0.5 * (ngrid - 1) * args.step
    geo_co['coords_bohr'] = (geo_co['coords_bohr'] * BOHR2ANG + grid_center) / BOHR2ANG
    atoms_dict_co['pos'] = geo_co['coords_bohr'] * BOHR2ANG

    rho_tip_total = project_dftb_density(geo_co, evec_co, projector_co, atoms_dict_co, grid_spec, basis)
    rho_tip_na = project_neutral_density(geo_co, projector_co, atoms_dict_co, grid_spec, basis)
    rho_tip_delta = rho_tip_total - rho_tip_na

    save_npy(os.path.join(debug_dir, 'co_tip'), 'co_rho_total.npy', rho_tip_total)
    save_npy(os.path.join(debug_dir, 'co_tip'), 'co_rho_delta.npy', rho_tip_delta)

    # Steps 3-5
    _, grads_pauli = step3_pauli(rho_grid, origin, args.step, rho_tip_total, debug_dir, A_pauli=A_pauli, beta_pauli=beta_pauli)
    _, grads_es    = step4_electrostatics_conv(V_ES, origin, args.step, rho_tip_delta, debug_dir)
    _, grads_vdw   = step5_dispersion(atomPos, atomTypes, origin, args.step, ngrid, debug_dir)

    # Step 6
    scan_xs = np.linspace(atomPos[:,0].min()-3, atomPos[:,0].max()+3, 100)
    scan_ys = np.linspace(atomPos[:,1].min()-3, atomPos[:,1].max()+3, 100)
    heights = np.arange(3.0, 5.5, 0.1)
    step6_composed(grads_pauli, grads_es, grads_vdw, origin, args.step, atomPos, scan_xs, scan_ys, heights, debug_dir)

    print(f"\nALL DONE. Results in {debug_dir}/")

if __name__ == "__main__":
    main()
