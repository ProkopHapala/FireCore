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

# ==============================================================================
# Constants & Config
# ==============================================================================

BOHR2ANG = 0.5291772109
COULOMB_CONST = 14.3996448915
ELEM_Z = {'H':1, 'C':6, 'N':7, 'O':8, 'P':15, 'S':16}

_XYZ_PATH = os.path.join(_THIS_DIR, 'pentacene.xyz')
_DEBUG_DIR = os.path.join(_THIS_DIR, 'debug_dftb_pentacene')  # Default, can be overridden
SLAKO_PREFIX = "/home/prokop/SIMULATIONS/dftbplus/slakos/mio-1-1/"
BASIS_HSD = "/home/prokop/git/FireCore/pyBall/DFTB/data/wfc.mio-1-1.hsd"

# Normalization factor for density in Angstrom grid
# If orbitals are normalized in Bohr, density integral in Angstrom is B^3.
# To get 1.0 in Angstrom, multiply density by B^-3.
B3_FACTOR = 1.0 / (BOHR2ANG**3) 

def setup_debug_dirs():
    if os.path.exists(_DEBUG_DIR):
        shutil.rmtree(_DEBUG_DIR)
    for step in ['step1_density', 'step2_electrostatics', 'step3_pauli',
                 'step4_electrostatics_conv', 'step5_dispersion', 'step6_composed', 'co_tip']:
        os.makedirs(os.path.join(_DEBUG_DIR, step), exist_ok=True)

def load_xyz(fname):
    with open(fname, 'r') as f:
        lines = f.readlines()
    natoms = int(lines[0])
    atomTypes = []; atomPos = []
    for line in lines[2:2+natoms]:
        p = line.split()
        sym = p[0]
        z = ELEM_Z[sym]
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

def scf_worker(atomTypes, atomPos, work_dir, queue):
    """Function to be run in a separate process."""
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    with open(work_dir / "geom.xyz", "w") as f:
        f.write(f"{len(atomTypes)}\nGenerated\n")
        inv_z = {v: k for k, v in ELEM_Z.items()}
        for i in range(len(atomTypes)):
            f.write(f"{inv_z[atomTypes[i]]} {atomPos[i,0]} {atomPos[i,1]} {atomPos[i,2]}\n")
    
    species_used = sorted(set(atomTypes.tolist()))
    inv_z = {v: k for k, v in ELEM_Z.items()}
    max_ang = "\n    ".join([f'{inv_z[z]} = "p"' if z > 1 else f'{inv_z[z]} = "s"' for z in species_used])

    hsd_content = f"""
Geometry = xyzFormat {{
  <<< "geom.xyz"
}}
Hamiltonian = DFTB {{
  SCC = Yes
  SlaterKosterFiles = Type2FileNames {{
    Prefix = "{SLAKO_PREFIX}"
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

def run_isolated_scf(atomTypes, atomPos, work_dir):
    q = Queue()
    p = Process(target=scf_worker, args=(atomTypes, atomPos, work_dir, q))
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

def step1_density_projection(atomTypes, atomPos, basis, step=0.15, margin=4.0, z_extra=6.0, stop_after=False):
    step_dir = os.path.join(_DEBUG_DIR, 'step1_density')
    log = lambda msg: log_step(step_dir, msg)
    log("="*60); log("STEP 1: Density Projection to Grid (DFTB+)"); log("="*60)
    
    log("\n1a. Running DFTB+ SCF (isolated process)...")
    t0 = time.time()
    geo, evecs = run_isolated_scf(atomTypes, atomPos, os.path.join(step_dir, "pentacene_work"))
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

def step2_electrostatics(rho_diff, step, origin, ngrid):
    step_dir = os.path.join(_DEBUG_DIR, 'step2_electrostatics')
    log = lambda msg: log_step(step_dir, msg)
    log("="*60); log("STEP 2: Electrostatics via Poisson"); log("="*60)
    V_ES = afm.fft_poisson(rho_diff, step)
    log(f"V_ES range: [{V_ES.min():.3f}, {V_ES.max():.3f}] eV")
    save_npy(step_dir, 'V_ES.npy', V_ES)
    plot_field_slices(V_ES, "Electrostatic Potential [eV]", "step2_VES_slices.png", step_dir, origin=origin, step=step, sym=True, cmap='bwr')
    return V_ES

def step3_pauli(rho_grid, origin, step, rho_tip_total, A_pauli=16.0):
    step_dir = os.path.join(_DEBUG_DIR, 'step3_pauli')
    log = lambda msg: log_step(step_dir, msg)
    log("="*60); log("STEP 3: Pauli Repulsion Convolution"); log("="*60)
    dV = step**3
    nx_t, ny_t, nz_t = rho_tip_total.shape
    tip_kernel = np.roll(np.roll(np.roll(rho_tip_total[::-1,::-1,::-1], -(nx_t//2), axis=0), -(ny_t//2), axis=1), -(nz_t//2), axis=2)
    E_pauli_field = A_pauli * dV * np.real(np.fft.ifftn(np.fft.fftn(rho_grid) * np.fft.fftn(tip_kernel))).astype(np.float32)
    grads_E_Pauli = np.stack([np.gradient(E_pauli_field, step, axis=i) for i in range(3)], axis=-1)
    save_npy(step_dir, 'E_Pauli_field.npy', E_pauli_field)
    plot_field_slices(E_pauli_field, "Pauli Energy [eV]", "step3_Epauli_slices.png", step_dir, origin=origin, step=step)
    return E_pauli_field, grads_E_Pauli

def step4_electrostatics_conv(V_ES, origin, step, rho_tip_delta):
    step_dir = os.path.join(_DEBUG_DIR, 'step4_electrostatics_conv')
    log = lambda msg: log_step(step_dir, msg)
    log("="*60); log("STEP 4: Electrostatics Convolution"); log("="*60)
    dV = step**3
    nx_t, ny_t, nz_t = rho_tip_delta.shape
    tip_kernel = np.roll(np.roll(np.roll(rho_tip_delta[::-1,::-1,::-1], -(nx_t//2), axis=0), -(ny_t//2), axis=1), -(nz_t//2), axis=2)
    E_es_field = dV * np.real(np.fft.ifftn(np.fft.fftn(V_ES) * np.fft.fftn(tip_kernel))).astype(np.float32)
    grads_E_ES = np.stack([np.gradient(E_es_field, step, axis=i) for i in range(3)], axis=-1)
    save_npy(step_dir, 'E_ES_field.npy', E_es_field)
    plot_field_slices(E_es_field, "ES Energy [eV]", "step4_EES_slices.png", step_dir, origin=origin, step=step, sym=True, cmap='bwr')
    return E_es_field, grads_E_ES

def step5_dispersion(atomPos, atomTypes, origin, step, ngrid, C6_CO=30.0):
    step_dir = os.path.join(_DEBUG_DIR, 'step5_dispersion')
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

def step6_composed(grads_pauli, grads_es, grads_vdw, origin, step, atomPos, scan_xs, scan_ys, heights, K_LAT=0.5):
    step_dir = os.path.join(_DEBUG_DIR, 'step6_composed')
    log = lambda msg: log_step(step_dir, msg)
    log("="*60); log("STEP 6: Final Composed & Relaxation"); log("="*60)
    
    # F = -grad E
    F_total = -(grads_pauli + grads_es + grads_vdw)
    
    nx_s, ny_s = len(scan_xs), len(scan_ys)
    nz_s = len(heights)
    
    from scipy.ndimage import map_coordinates
    def force_func(positions):
        """Interpolate forces at arbitrary positions from scan-grid force field."""
        # Map positions to grid fractional coordinates
        ix = (positions[:, 0] - origin[0]) / step
        iy = (positions[:, 1] - origin[1]) / step
        iz = (positions[:, 2] - origin[2]) / step
        
        coords = np.vstack([ix, iy, iz])
        fx = map_coordinates(F_total[..., 0], coords, order=1)
        fy = map_coordinates(F_total[..., 1], coords, order=1)
        fz = map_coordinates(F_total[..., 2], coords, order=1)
        return np.stack([fx, fy, fz], axis=-1)
    
    log("\n6b. Probe particle relaxation...")
    # pp_relax_2d(force_func, scan_xs, scan_ys, heights, mol_z=0.0, K_LAT=0.5, N_RELAX=100, step=0.1)
    FEs_relax = afm.pp_relax_2d(force_func, scan_xs, scan_ys, heights, mol_z=atomPos[:,2].max(), K_LAT=K_LAT, N_RELAX=50, step=step)
    Fz_relax = FEs_relax[:,:,:,2]
    
    df = afm.compute_df(Fz_relax, heights[1]-heights[0])
    save_npy(step_dir, 'df.npy', df)
    
    for i in [0, len(heights)//2, -1]:
        h = heights[i]
        plt.figure(figsize=(6,5))
        plt.imshow(df[:,:,i].T, origin='lower', extent=[scan_xs[0], scan_xs[-1], scan_ys[0], scan_ys[-1]], cmap='afmhot')
        plt.title(f"df at h={h:.1f} A")
        plt.colorbar(label="df [Hz]")
        plt.savefig(os.path.join(step_dir, f"df_h{h:.1f}.png"), dpi=150, bbox_inches='tight')
        plt.close()
        
    log("Step 6 COMPLETE.")
    return df

# ==============================================================================
# Main
# ==============================================================================

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--step', type=float, default=0.15)
    parser.add_argument('--output_dir', type=str, default=None, help='Custom output directory name (e.g., debug_dftb_pentacene_step0.1)')
    parser.add_argument('--basis', type=str, default='mio-1-1', help='Basis set: mio-1-1 or 3ob-3-1')
    args = parser.parse_args()
    
    # Override debug directory if custom output_dir specified
    global _DEBUG_DIR, SLAKO_PREFIX, BASIS_HSD
    if args.output_dir:
        _DEBUG_DIR = os.path.join(_THIS_DIR, args.output_dir)
    
    # Set basis-specific paths
    if args.basis == '3ob-3-1':
        SLAKO_PREFIX = "/home/prokop/SIMULATIONS/dftbplus/slakos/3ob-3-1/"
        BASIS_HSD = "/home/prokop/git/dftbplus/tests/grid/dftb_ptcda/wfc.3ob-3-1.hsd"
    elif args.basis == 'mio-1-1':
        SLAKO_PREFIX = "/home/prokop/SIMULATIONS/dftbplus/slakos/mio-1-1/"
        BASIS_HSD = "/home/prokop/git/dftbplus/tests/grid/dftb_ptcda/waveplot_in.hsd"
    else:
        raise ValueError(f"Unknown basis: {args.basis}. Use 'mio-1-1' or '3ob-3-1'")
    
    print(f"Using basis: {args.basis}")
    print(f"SLAKO_PREFIX: {SLAKO_PREFIX}")
    print(f"BASIS_HSD: {BASIS_HSD}")
    
    setup_debug_dirs()
    atomTypes, atomPos = load_xyz(_XYZ_PATH)
    
    # Load Basis once
    if args.basis == '3ob-3-1':
        # Use parse_wfc_hsd for 3ob-3-1 (wfc.*.hsd format) and convert
        basis_data = parse_wfc_hsd(BASIS_HSD)
        basis = convert_wfc_to_species_list_ang(basis_data, resolution_bohr=0.04)
    else:
        # Use parse_basis_hsd_ang for mio-1-1 (waveplot_in.hsd format)
        basis = parse_basis_hsd_ang(BASIS_HSD)
    print(f"Loaded basis for species: {[sp['name'] for sp in basis]}")
    
    # Step 1-2: Pentacene
    rho_grid, rho_na_grid, rho_diff, origin, ngrid, grid_spec, atomPos, atomTypes = \
        step1_density_projection(atomTypes, atomPos, basis, step=args.step)
    V_ES = step2_electrostatics(rho_diff, args.step, origin, ngrid)
    
    # CO Tip Calculation (isolated process)
    log_step(os.path.join(_DEBUG_DIR, 'co_tip'), "Computing CO Tip density...")
    geo_co, evec_co = run_isolated_scf(np.array([8, 6]), np.array([[0,0,0], [0,0,1.13]]), os.path.join(_DEBUG_DIR, 'co_tip', 'co_work'))
    dftb_co = {'coords_bohr': geo_co['coords_bohr'], 'species_per_atom': geo_co['species_per_atom'], 'species_names': geo_co['species_names']}
    projector_co, atoms_dict_co = setup_gridprojector_from_dftb(dftb_co, basis, verbosity=0)
    
    # We must ensure O is at grid center for convolution to work properly!
    # Tip atoms are projected at their positions within grid_spec.
    # In run_pyocl_fdbm.py, get_co_tip_density centers CO at grid center.
    grid_center = origin + 0.5 * (ngrid - 1) * args.step
    geo_co['coords_bohr'] = (geo_co['coords_bohr'] * BOHR2ANG + grid_center) / BOHR2ANG
    atoms_dict_co['pos'] = geo_co['coords_bohr'] * BOHR2ANG
    
    rho_tip_total = project_dftb_density(geo_co, evec_co, projector_co, atoms_dict_co, grid_spec, basis)
    rho_tip_na = project_neutral_density(geo_co, projector_co, atoms_dict_co, grid_spec, basis)
    rho_tip_delta = rho_tip_total - rho_tip_na
    
    save_npy(os.path.join(_DEBUG_DIR, 'co_tip'), 'co_rho_total.npy', rho_tip_total)
    save_npy(os.path.join(_DEBUG_DIR, 'co_tip'), 'co_rho_delta.npy', rho_tip_delta)
    
    # Steps 3-5
    _, grads_pauli = step3_pauli(rho_grid, origin, args.step, rho_tip_total)
    _, grads_es    = step4_electrostatics_conv(V_ES, origin, args.step, rho_tip_delta)
    _, grads_vdw   = step5_dispersion(atomPos, atomTypes, origin, args.step, ngrid)
    
    # Step 6
    scan_xs = np.linspace(atomPos[:,0].min()-3, atomPos[:,0].max()+3, 100)
    scan_ys = np.linspace(atomPos[:,1].min()-3, atomPos[:,1].max()+3, 100)
    heights = np.arange(3.0, 5.5, 0.1)
    step6_composed(grads_pauli, grads_es, grads_vdw, origin, args.step, atomPos, scan_xs, scan_ys, heights)
    
    print("\nALL DONE. Results in debug_dftb_pentacene/")

if __name__ == "__main__":
    main()
