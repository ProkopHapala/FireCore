#!/usr/bin/env python3
"""
test_3d_grid_density.py

Rigorous 3D grid density test comparing:
- Method A: WavePlot.orb2grid (libwaveplot.so reference)
- Method B: Grid_dftb.project_density_dense (OpenCL)

Tests:
1. Electron count integration vs expected valence electrons
2. Cube file export for visual inspection (VESTA)
3. Grid subtraction and maximum error reporting

Usage:
    cd tests/dftb
    python test_3d_grid_density.py --xyz tests/tAFM/pyocl_fdbm/TBTAP_3mols_c3h.xyz --basis 3ob-3-1 --step 0.2
"""

import sys, os, argparse, tempfile, shutil
import numpy as np
from pathlib import Path

REPO_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

from pyBall import atomicUtils as au
from pyBall.DFTB.DFTBcore import DFTBcore
from pyBall.DFTB.DFTBplusParser import parse_wfc_hsd, convert_wfc_to_species_list_ang, build_wp_basis
from pyBall.DFTB.Grid_dftb import setup_gridprojector_from_dftb
from pyBall.DFTB.TestUtils import write_cube

BOHR2ANG = 0.5291772109
ANG2BOHR = 1.0 / BOHR2ANG
ELEM_Z = {'H': 1, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'P': 15, 'S': 16, 'Br': 35}
VALENCE = {'H': 1, 'C': 4, 'N': 5, 'O': 6, 'F': 7, 'P': 5, 'S': 6, 'Br': 7}
SK_PREFIX = os.environ.get('DFTB_SK_PATH', '')

# Default libwaveplot.so path (override via --libpath or WAVEPLOT_LIB env)
_DEFAULT_WAVEPLOT_LIB = os.environ.get(
    'WAVEPLOT_LIB',
    '/home/prokop/git/dftbplus/_build/app/waveplot/libwaveplot.so'
)
B3_FACTOR = 1.0 / (BOHR2ANG ** 3)


def parse_args():
    p = argparse.ArgumentParser(
        description='Rigorous 3D grid density test: WavePlot vs Grid_dftb',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    p.add_argument('--xyz', type=str, required=True, help='Input XYZ geometry file')
    p.add_argument('--basis', type=str, default='3ob-3-1', choices=['mio-1-1', '3ob-3-1'])
    p.add_argument('--step', type=float, default=0.2, help='Grid step in Angstrom')
    p.add_argument('--padding', type=float, default=3.0, help='Padding around molecule in Angstrom')
    p.add_argument('--output-dir', type=str, default=None, help='Output directory for cube files')
    p.add_argument('--libpath', type=str, default=_DEFAULT_WAVEPLOT_LIB, help='Path to libwaveplot.so')
    p.add_argument('--skip-waveplot', action='store_true', help='Skip WavePlot reference (run Grid_dftb only)')
    return p.parse_args()


def write_dftb_input(work_dir, enames, coords_ang, basis_data):
    """Write minimal DFTB+ input file for SCF."""
    xyz_path = os.path.join(work_dir, 'geom.xyz')
    hsd_path = os.path.join(work_dir, 'dftb_in.hsd')
    au.save_xyz(xyz_path, enames, coords_ang)

    species = sorted(set(enames))
    max_am = {}
    for elem in species:
        sp = basis_data[elem]
        max_l = max(orb['AngularMomentum'] for orb in sp['orbitals'])
        max_am[elem] = {0: 's', 1: 'p', 2: 'd'}[max_l]

    sk_dir = os.path.join(SK_PREFIX, args.basis)
    with open(hsd_path, 'w') as f:
        f.write(f'Geometry = xyzFormat {{\n')
        f.write(f'    <<< "geom.xyz"\n')
        f.write(f'}}\n\n')
        f.write(f'Hamiltonian = DFTB {{\n')
        f.write(f'    SCC = Yes\n')
        f.write(f'    SCCTolerance = 1e-7\n')
        f.write(f'    MaxSCCIterations = 200\n')
        f.write(f'    SlaterKosterFiles = Type2FileNames {{\n')
        f.write(f'        Prefix = "{sk_dir}/"\n')
        f.write(f'        Separator = "-"\n')
        f.write(f'        Suffix = ".skf"\n')
        f.write(f'        LowerCaseTypeName = No\n')
        f.write(f'    }}\n')
        f.write(f'    MaxAngularMomentum = {{\n')
        for elem in species:
            f.write(f'        {elem} = "{max_am[elem]}"\n')
        f.write(f'    }}\n')
        f.write(f'}}\n')
    return hsd_path


def copy_sk_files(work_dir, enames):
    """Copy required SK files to work directory."""
    sk_dir = os.path.join(SK_PREFIX, args.basis)
    if not os.path.exists(sk_dir):
        raise RuntimeError(f"SK directory not found: {sk_dir}. Set DFTB_SK_PATH env var.")
    species = sorted(set(enames))
    for i, elem1 in enumerate(species):
        for elem2 in species[i:]:
            for sk_file in [f"{elem1}-{elem2}.skf", f"{elem2}-{elem1}.skf"]:
                src = os.path.join(sk_dir, sk_file)
                if os.path.exists(src):
                    shutil.copy(src, work_dir)


def build_orbital_layout(basis_data, enames):
    """Build norb_per_atom and orb_offsets from basis."""
    norb_per_atom = []
    orb_offsets = [0]
    max_l = 0
    for name in enames:
        sp = basis_data[name]
        norb = sum(2 * orb['AngularMomentum'] + 1 for orb in sp['orbitals'])
        for orb in sp['orbitals']:
            max_l = max(max_l, orb['AngularMomentum'])
        norb_per_atom.append(norb)
        orb_offsets.append(orb_offsets[-1] + norb)
    return (np.array(norb_per_atom, dtype=np.int32),
            np.array(orb_offsets, dtype=np.int32),
            max_l)


def setup_waveplot_reference(coords_bohr, enames, species_list_ang, eigvecs, libpath):
    """Configure WavePlot instance from parsed DFTB+ data."""
    from pyBall.DFTB.WavePlot import WavePlot

    unique_species = sorted(set(enames))
    sp_name_to_idx = {name: i + 1 for i, name in enumerate(unique_species)}
    species_wp = np.array([sp_name_to_idx[name] for name in enames], dtype=np.int32)

    wp_basis, resoln_b = build_wp_basis(species_list_ang, unique_species)

    wp = WavePlot(libpath)
    wp.set_geometry(coords_bohr, species_wp, is_periodic=False)
    wp.set_basis(wp_basis, resolution=resoln_b)
    wp.set_eigenvectors(eigvecs)
    return wp


def compute_waveplot_density(wp, nocc, origin_b, gridVecs_b, nPoints):
    """Accumulate electron density from occupied MOs using WavePlot.orb2grid.

    Returns density in atomic units (a0^-3).
    """
    rho = np.zeros(tuple(nPoints), dtype=np.float64)
    for imo in range(nocc):
        psi = wp.orb2grid(imo + 1, origin_b, gridVecs_b, nPoints)
        rho += 2.0 * (psi ** 2)
    return rho.astype(np.float32)


def main():
    global args
    args = parse_args()

    system_name = Path(args.xyz).stem
    if args.output_dir:
        output_dir = Path(args.output_dir)
    else:
        output_dir = Path('.') / f'cube_{system_name}_{args.basis}_s{args.step:.2f}'
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("3D Grid Density Test: WavePlot vs Grid_dftb")
    print("=" * 70)

    # 1. Load geometry
    print(f"\nLoading {args.xyz}")
    coords_ang, _, enames, _, _ = au.load_xyz(args.xyz)
    coords_ang = np.array(coords_ang, dtype=np.float64)
    coords_bohr = coords_ang * ANG2BOHR
    natoms = len(enames)
    print(f"  {natoms} atoms: {set(enames)}")

    # 2. Expected valence electrons
    n_electrons_expected = sum(VALENCE.get(n, 0) for n in enames)
    nocc = n_electrons_expected // 2
    print(f"  Expected valence electrons: {n_electrons_expected}")
    print(f"  Occupied MOs: {nocc}")

    # 3. Load basis
    basis_path = REPO_ROOT / 'pyBall/DFTB/data' / f'wfc.{args.basis}.hsd'
    print(f"\nLoading basis from {basis_path}")
    basis_data = parse_wfc_hsd(str(basis_path))
    species_list_ang = convert_wfc_to_species_list_ang(basis_data, resolution_bohr=0.04)
    print(f"  Species: {[sp['name'] for sp in species_list_ang]}")

    # 4. Build orbital layout
    norb_per_atom, orb_offsets, max_l = build_orbital_layout(basis_data, enames)
    max_shells = 3 if max_l >= 2 else 2
    norb_total = orb_offsets[-1]
    print(f"  Orbital layout: max_l={max_l}, max_shells={max_shells}, norb_total={norb_total}")

    # 5. Create work directory and setup DFTB
    work_dir = tempfile.mkdtemp(prefix='dftb_dense_')
    print(f"\nWork dir: {work_dir}")
    copy_sk_files(work_dir, enames)
    hsd_path = write_dftb_input(work_dir, enames, coords_ang, basis_data)

    # 6. Run DFTBcore SCF
    print("\nRunning DFTB SCF via DFTBcore...")
    old_cwd = os.getcwd()
    os.chdir(work_dir)
    try:
        dftb = DFTBcore()
        dftb.init('dftb_in.hsd')
        dftb.enable_matrix_collection(dm=True, h=False, s=False)
        energy = dftb.run_scf()
        dm = dftb.get_dm_dense()
        eigvecs, eigvals = dftb.get_eigvecs_dense()
        basis_size = dftb.get_basis_size()
        dftb.finalize()
    finally:
        os.chdir(old_cwd)
    print(f"  Energy: {energy:.6f} Ha")
    print(f"  Basis size: {basis_size}")
    print(f"  DM shape: {dm.shape}")

    # 7. Define 3D grid
    padding = args.padding
    rmin = coords_ang.min(axis=0) - padding
    rmax = coords_ang.max(axis=0) + padding
    ngrid = np.ceil((rmax - rmin) / args.step).astype(np.int32)
    ngrid = np.maximum(ngrid, 2)  # ensure at least 2 points per axis
    nx, ny, nz = ngrid
    print(f"\nGrid: {nx}x{ny}x{nz} = {nx*ny*nz:,} points")
    print(f"  Step: {args.step:.3f} Å")
    print(f"  Range x: [{rmin[0]:.2f}, {rmax[0]:.2f}] Å")
    print(f"  Range y: [{rmin[1]:.2f}, {rmax[1]:.2f}] Å")
    print(f"  Range z: [{rmin[2]:.2f}, {rmax[2]:.2f}] Å")

    dV_ang = args.step ** 3
    dV_bohr = (args.step * ANG2BOHR) ** 3

    # Grid spec for Grid_dftb (Angstrom)
    grid_spec = {
        'origin': rmin.tolist(),
        'dA': [args.step, 0.0, 0.0],
        'dB': [0.0, args.step, 0.0],
        'dC': [0.0, 0.0, args.step],
        'ngrid': ngrid.tolist()
    }

    # Grid spec for WavePlot (Bohr)
    origin_b = rmin * ANG2BOHR
    step_b = args.step * ANG2BOHR
    gridVecs_b = np.diag([step_b, step_b, step_b])
    nPoints = ngrid.astype(np.int32)

    # Atoms for cube file (Bohr)
    atoms_bohr = [
        (ELEM_Z[name], coords_bohr[i, 0], coords_bohr[i, 1], coords_bohr[i, 2])
        for i, name in enumerate(enames)
    ]

    # ================================================================
    # Method A: WavePlot (libwaveplot.so)
    # ================================================================
    rho_wp_bohr = None
    rho_wp_ang = None
    if not args.skip_waveplot:
        if not os.path.exists(args.libpath):
            print(f"\nWARNING: libwaveplot.so not found at {args.libpath}")
            print("  Set --libpath or WAVEPLOT_LIB env var. Skipping WavePlot reference.")
        else:
            print("\n" + "=" * 70)
            print("Method A: WavePlot.orb2grid")
            print("=" * 70)
            wp = setup_waveplot_reference(coords_bohr, enames, species_list_ang, eigvecs, args.libpath)
            rho_wp_bohr = compute_waveplot_density(wp, nocc, origin_b, gridVecs_b, nPoints)
            # Convert to e/Ang^3 for comparison with Grid_dftb
            rho_wp_ang = (rho_wp_bohr * B3_FACTOR).astype(np.float32)
            n_wp_bohr = float(np.sum(rho_wp_bohr) * dV_bohr)
            n_wp_ang = float(np.sum(rho_wp_ang) * dV_ang)
            print(f"  Integrated electrons (Bohr): {n_wp_bohr:.4f}")
            print(f"  Integrated electrons (Ang):  {n_wp_ang:.4f}")
            print(f"  Max density: {rho_wp_bohr.max():.4e} a0^-3")
            print(f"  Max density: {rho_wp_ang.max():.4e} e/Ang^3")

    # ================================================================
    # Method B: Grid_dftb.project_density_dense
    # ================================================================
    print("\n" + "=" * 70)
    print("Method B: Grid_dftb.project_density_dense")
    print("=" * 70)

    dftb_data = {
        'coords_bohr': coords_bohr,
        'species_per_atom': np.array([sorted(set(enames)).index(n) for n in enames], dtype=np.int32),
        'species_names': sorted(set(enames))
    }
    projector, atoms_dict = setup_gridprojector_from_dftb(
        dftb_data, species_list_ang, verbosity=1, max_shells=max_shells
    )

    dm32 = dm.astype(np.float32)
    rho_ocl = projector.project_density_dense(
        dm32, norb_per_atom, orb_offsets, atoms_dict, grid_spec, nMaxAtom=64
    )
    # project_density_dense now returns e/Ang^3 (B3_FACTOR applied internally)
    n_ocl_ang = float(np.sum(rho_ocl) * dV_ang)
    print(f"  Integrated electrons (Ang): {n_ocl_ang:.4f}")
    print(f"  Max density: {rho_ocl.max():.4e} e/Ang^3")

    # ================================================================
    # Comparison
    # ================================================================
    print("\n" + "=" * 70)
    print("Comparison")
    print("=" * 70)

    print(f"\n{'Quantity':<30} {'Value':>15}")
    print("-" * 50)
    print(f"{'Expected valence electrons':<30} {n_electrons_expected:>15.2f}")
    if rho_wp_ang is not None:
        print(f"{'WavePlot electrons (Ang)':<30} {n_wp_ang:>15.4f}")
        rel_err_wp = abs(n_wp_ang - n_electrons_expected) / n_electrons_expected * 100.0
        print(f"{'WavePlot relative error (%)':<30} {rel_err_wp:>15.2f}")
    print(f"{'Grid_dftb electrons (Ang)':<30} {n_ocl_ang:>15.4f}")
    rel_err_ocl = abs(n_ocl_ang - n_electrons_expected) / n_electrons_expected * 100.0
    print(f"{'Grid_dftb relative error (%)':<30} {rel_err_ocl:>15.2f}")

    if rho_wp_ang is not None:
        diff = rho_wp_ang - rho_ocl
        max_err = float(np.max(np.abs(diff)))
        rms_err = float(np.sqrt(np.mean(diff ** 2)))
        ref_max = float(np.max(np.abs(rho_wp_ang)))
        rel_rms = (rms_err / ref_max * 100.0) if ref_max > 1e-20 else 0.0
        print(f"\nWavePlot vs Grid_dftb:")
        print(f"  Max absolute error: {max_err:.4e} e/Ang^3")
        print(f"  RMS error:            {rms_err:.4e} e/Ang^3")
        print(f"  Relative RMS:         {rel_rms:.4f} %")

    # ================================================================
    # Cube file export
    # ================================================================
    print("\n" + "=" * 70)
    print("Cube File Export")
    print("=" * 70)

    # Write in Bohr units (standard quantum chemistry convention)
    if rho_wp_bohr is not None:
        cube_wp = output_dir / f"rho_{system_name}_{args.basis}_waveplot.cube"
        write_cube(
            str(cube_wp), rho_wp_bohr, atoms_bohr, origin_b, gridVecs_b,
            title=f"WavePlot density {system_name} {args.basis}"
        )
        print(f"  Saved: {cube_wp}")

    # Grid_dftb density: convert e/Ang^3 back to a0^-3 for cube consistency
    rho_ocl_bohr = (rho_ocl / B3_FACTOR).astype(np.float32)
    cube_ocl = output_dir / f"rho_{system_name}_{args.basis}_grid_dftb.cube"
    write_cube(
        str(cube_ocl), rho_ocl_bohr, atoms_bohr, origin_b, gridVecs_b,
        title=f"Grid_dftb density {system_name} {args.basis}"
    )
    print(f"  Saved: {cube_ocl}")

    if rho_wp_bohr is not None:
        diff_bohr = rho_wp_bohr - rho_ocl_bohr
        cube_diff = output_dir / f"diff_{system_name}_{args.basis}_wp_minus_ocl.cube"
        write_cube(
            str(cube_diff), diff_bohr, atoms_bohr, origin_b, gridVecs_b,
            title=f"Density difference {system_name} {args.basis}"
        )
        print(f"  Saved: {cube_diff}")

    # Cleanup
    shutil.rmtree(work_dir)
    print(f"\nCleaned up work dir")
    print(f"Output directory: {output_dir.absolute()}")
    print("=" * 70)
    print("Done!")
    print("=" * 70)
    return 0


if __name__ == '__main__':
    sys.exit(main())
