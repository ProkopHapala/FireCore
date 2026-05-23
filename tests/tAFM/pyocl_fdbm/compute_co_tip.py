#!/usr/bin/env python3
"""
Standalone script to compute CO tip density via SCF + PyOpenCL projection.
Supports both FireCore and DFTB backends.
Must run in separate process because Fireball cannot reallocate for different molecules.

Usage:
    python compute_co_tip.py <output_dir> <grid_spec_json> <step> <nscf> <fdata_dir> <fdata_basis> [backend]

Outputs:
    <output_dir>/co_rho_total.npy
    <output_dir>/co_rho_delta.npy
"""
import sys, os, json, numpy as np

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.realpath(os.path.join(_THIS_DIR, '..', '..', '..'))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

# Parse backend argument (default to 'firecore' for backward compatibility)
backend = 'firecore'
if len(sys.argv) >= 8:
    backend = sys.argv[7].lower()

# Conditional imports based on backend
if backend == 'firecore':
    try:
        from pyBall import FireCore as fc
        from pyBall.FireballOCL import Grid as ocl_grid
    except (OSError, ImportError) as e:
        raise RuntimeError(f"[ERROR] FireCore backend requested but libFireCore.so not found: {e}")
        # print(f"[INFO] Falling back to DFTB backend...")
        # backend = 'dftb'
        # from pyBall.DFTB.DFTBcore import DFTBcore
        # from pyBall.DFTB import Grid_dftb as dg
elif backend == 'dftb':
    from pyBall.DFTB.DFTBcore import DFTBcore
    from pyBall.DFTB import Grid_dftb as dg
else:
    raise ValueError(f"Unknown backend: {backend}. Use 'firecore' or 'dftb'.")

Z_TO_ZVAL = {1: 1, 6: 4, 7: 5, 8: 6, 16: 6}
RCUT_DEFAULT = {1: 2.3, 6: 2.6, 7: 2.6, 8: 2.5}

def _onsite_occ(Z):
    if Z == 1:  return np.array([1.0, 0.0, 0.0, 0.0], dtype=np.float32)
    if Z == 6:  return np.array([2.0, 2.0/3, 2.0/3, 2.0/3], dtype=np.float32)
    if Z == 7:  return np.array([2.0, 1.0, 1.0, 1.0], dtype=np.float32)
    if Z == 8:  return np.array([2.0, 4.0/3, 4.0/3, 4.0/3], dtype=np.float32)
    return np.array([float(Z_TO_ZVAL.get(int(Z), int(Z))), 0.0, 0.0, 0.0], dtype=np.float32)

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

def main():
    if len(sys.argv) < 7:
        print("Usage: python compute_co_tip.py <out_dir> <grid_spec_json> <step> <nscf> <fdata_dir> <fdata_basis> [backend]")
        sys.exit(1)
    
    out_dir = sys.argv[1]
    grid_spec = json.loads(sys.argv[2])
    # Convert lists back to numpy arrays
    for key in ['origin', 'dA', 'dB', 'dC', 'ngrid']:
        if key in grid_spec:
            grid_spec[key] = np.array(grid_spec[key], dtype=np.float32 if key.startswith('d') else np.int32 if key=='ngrid' else np.float32)
    
    step = float(sys.argv[3])
    nscf = int(sys.argv[4])
    fdata_dir = sys.argv[5]
    fdata_basis = sys.argv[6]
    
    os.makedirs(out_dir, exist_ok=True)
    
    co_xyz = os.path.join(_ROOT, 'cpp', 'common_resources', 'xyz', 'CO.xyz')
    co_types, co_pos = load_xyz(co_xyz)
    
    # Center CO at grid center (so O is at center for FFT convolution)
    # FFT convolution places tip at array center, so O must be at grid center
    grid_center = np.array(grid_spec['origin']) + 0.5 * (np.array(grid_spec['ngrid']) - 1) * step
    # CO has O at origin in xyz, so shift all CO positions to grid_center
    co_pos = co_pos + grid_center
    print(f"  Grid center: {grid_center}")
    print(f"  CO positions (centered): O={co_pos[0]}, C={co_pos[1]}")
    print(f"  Using backend: {backend}")
    
    if backend == 'firecore':
        # Setup Fdata symlink
        _fball_cwd = os.path.join(_ROOT, 'tests', 'pyFireball')
        orig_cwd = os.getcwd()
        os.chdir(_fball_cwd)
        
        fdata_local = os.path.join(_fball_cwd, 'Fdata')
        if not os.path.exists(os.path.join(fdata_local, 'info.dat')):
            if os.path.lexists(fdata_local):
                os.unlink(fdata_local)
            os.symlink(fdata_dir, fdata_local)
        
        # Run Fireball SCF on CO
        print(f"[CO Tip Subprocess] Running Fireball SCF on CO (nscf={nscf})...")
        fc.setVerbosity(0)
        fc.preinit()
        fc.init(co_types, co_pos)
        fc.SCF(co_pos, nmax_scf=nscf)
        
        dims = fc.get_HS_dims()
        neighs = fc.get_HS_neighs(dims)
        neighs = fc.get_rho_sparse(dims, data=neighs)
        rho_sparse_co = neighs.rho
        
        # Build neutral atom density matrix
        rho_na_co = np.zeros_like(rho_sparse_co, dtype=np.float32)
        neigh_j = neighs.neigh_j.reshape(len(co_types), -1)
        for i in range(len(co_types)):
            slots = np.where(neigh_j[i] == (i+1))[0]
            if len(slots) == 0:
                raise RuntimeError(f"No self-neighbor for CO atom i={i}")
            iself = int(slots[0])
            occ = _onsite_occ(int(co_types[i]))
            rho_na_co[i, iself, :, :] = 0.0
            for k in range(4):
                rho_na_co[i, iself, k, k] = occ[k]
        
        # Project to grid
        projector = ocl_grid.GridProjector(fdata_dir=fdata_basis, verbosity=0)
        projector.load_basis(sorted(set(co_types.tolist())))
        atoms_dict = {
            'pos': co_pos,
            'Rcut': np.array([RCUT_DEFAULT.get(int(z), 4.5) for z in co_types]),
            'type': co_types,
        }
        
        print("  Projecting CO total density...")
        rho_total = projector.project(rho_sparse_co, neighs, atoms_dict, grid_spec, nMaxAtom=64, use_tiled=True)
        
        print("  Projecting CO neutral-atom density...")
        rho_na_grid = projector.project(rho_na_co, neighs, atoms_dict, grid_spec, nMaxAtom=64, use_tiled=True)
        
        os.chdir(orig_cwd)
    
    elif backend == 'dftb':
        # Use DFTB backend with dftb_utils functions
        from pyBall import dftb_utils as du
        from pyBall.DFTB.DFTBplusParser import parse_wfc_hsd, convert_wfc_to_species_list_ang
        
        # Use mio-1-1 basis for CO (C and O)
        basis_name = 'mio-1-1'
        sk_prefix = du.get_sk_path(basis_name)
        basis_hsd_path = du.WFC_HSD_PATHS.get(basis_name)
        if basis_hsd_path is None:
            raise RuntimeError(f"Basis HSD file not found for {basis_name}. Check WFC_HSD_PATHS in dftb_utils.py")
        
        # Setup work directory and run DFTB SCF using utility function
        work_dir = os.path.join(out_dir, 'dftb_work')
        enames = ['C' if z == 6 else 'O' if z == 8 else 'H' for z in co_types]
        
        print(f"[CO Tip Subprocess] Running DFTB SCF on CO (nscf={nscf})...")
        geo, evecs = du.run_dftb_for_density(work_dir, enames, co_pos, sk_prefix, xyz_fname='geom.xyz')
        
        # Get density matrix from eigenvectors (rho = sum_k f_k * C_k * C_k^T)
        occupations = geo['occupations']
        norb = geo['norb']
        dm_dense = np.zeros((norb, norb), dtype=np.float32)
        for k in range(len(occupations)):
            if occupations[k] > 0:
                dm_dense += occupations[k] * np.outer(evecs[:, k], evecs[:, k])
        
        # Load basis and setup projector
        basis_data = parse_wfc_hsd(basis_hsd_path)
        basis_ang = convert_wfc_to_species_list_ang(basis_data, resolution_bohr=0.04)
        
        # Build orbital layout
        norb_per_atom = []
        orb_offsets = [0]
        for name in enames:
            sp = basis_data[name]
            norb = sum(2 * orb['AngularMomentum'] + 1 for orb in sp['orbitals'])
            norb_per_atom.append(norb)
            orb_offsets.append(orb_offsets[-1] + norb)
        norb_per_atom = np.array(norb_per_atom, dtype=np.int32)
        orb_offsets = np.array(orb_offsets, dtype=np.int32)
        
        coords_bohr = co_pos * 1.8897259886
        species_per_atom = list(range(len(enames)))
        dftb_data = {
            'coords_bohr': coords_bohr,
            'species_per_atom': species_per_atom,
            'species_names': enames
        }
        
        projector, atoms_dict = dg.setup_gridprojector_from_dftb(dftb_data, basis_ang, verbosity=0, max_shells=2)
        
        # Project SCF density
        print("  Projecting CO total density (DFTB)...")
        rho_total = projector.project_density_dense(dm_dense.astype(np.float32), norb_per_atom, orb_offsets, atoms_dict, grid_spec)
        
        # Project neutral atom density
        print("  Projecting CO neutral-atom density (DFTB)...")
        geo_dict = {
            'natoms': len(enames),
            'species_per_atom': species_per_atom,
            'species_names': enames,
            'coords_bohr': coords_bohr
        }
        rho_na_grid = dg.project_neutral_density(geo_dict, projector, atoms_dict, grid_spec, basis_ang)
    
    rho_delta = (rho_total - rho_na_grid).astype(np.float32)
    
    # NOTE: CO density is kept at grid center (not shifted to index 0).
    # The FFT circular shift for convolution is applied in step3/step4.
    nx, ny, nz = rho_total.shape
    dV = step**3
    q_total = rho_total.sum() * dV
    q_delta = rho_delta.sum() * dV
    print(f"  CO rho_total: shape={rho_total.shape} range=[{rho_total.min():.4f},{rho_total.max():.4f}] q={q_total:.3f}")
    print(f"  CO rho_delta: range=[{rho_delta.min():.4f},{rho_delta.max():.4f}] q={q_delta:.3f}")
    print(f"  CO centered at grid center ({nx//2},{ny//2},{nz//2})")
    
    np.save(os.path.join(out_dir, 'co_rho_total.npy'), rho_total)
    np.save(os.path.join(out_dir, 'co_rho_delta.npy'), rho_delta)
    print(f"  Saved CO densities to {out_dir}")

if __name__ == '__main__':
    main()
