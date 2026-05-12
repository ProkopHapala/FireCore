#!/usr/bin/env python3
"""
test_dense_projection.py

End-to-end DFTB dense projection test.

Workflow: XYZ -> DFTBcore SCF -> density matrix -> Grid_dftb projection -> plots

Uses direct library access (DFTBcore) instead of file-based DFTB+ execution.
No intermediate files (detailed.xml, eigenvec.bin) required.

Usage:
    python test_dense_projection.py --xyz tests/tAFM/pyocl_fdbm/TBTAP_3mols_c3h.xyz --basis 3ob-3-1 --z-offset 2.0
    python test_dense_projection.py --xyz tests/tAFM/pyocl_fdbm/TBTAP_3mols_c3h.xyz --basis 3ob-3-1 --step 0.2 --z-offset 2.5 --no-plot
"""

import sys, os, argparse, tempfile, shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

REPO_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

from pyBall import atomicUtils as au
from pyBall.DFTB.DFTBcore import DFTBcore
from pyBall.DFTB.DFTBplusParser import parse_wfc_hsd, convert_wfc_to_species_list_ang
from pyBall.DFTB.Grid_dftb import setup_gridprojector_from_dftb
from pyBall.DFTB.TestUtils import generate_2d_point_grid

BOHR2ANG = 0.5291772109
ANG2BOHR = 1.0 / BOHR2ANG
ELEM_Z = {'H': 1, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'P': 15, 'S': 16, 'Br': 35}
VALENCE = {'H': 1, 'C': 4, 'N': 5, 'O': 6, 'F': 7, 'P': 5, 'S': 6, 'Br': 7}
SK_PREFIX = os.environ.get('DFTB_SK_PATH', '')


def parse_args():
    p = argparse.ArgumentParser(description='End-to-end DFTB dense projection', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('--xyz', type=str, required=True, help='Input XYZ geometry file')
    p.add_argument('--basis', type=str, default='3ob-3-1', choices=['mio-1-1', '3ob-3-1'])
    p.add_argument('--step', type=float, default=0.1, help='Grid step in Angstrom')
    p.add_argument('--z-offset', type=float, default=2.0, help='Z offset for XY plane')
    p.add_argument('--no-plot', action='store_true', help='Disable plot generation')
    p.add_argument('--output-dir', type=str, default=None, help='Output directory')
    p.add_argument('--dpi', type=int, default=150)
    return p.parse_args()


def write_dftb_input(work_dir, enames, coords_ang, basis_data):
    """Write minimal DFTB+ input file for SCF."""
    xyz_path = os.path.join(work_dir, 'geom.xyz')
    hsd_path = os.path.join(work_dir, 'dftb_in.hsd')
    au.save_xyz(xyz_path, enames, coords_ang)

    species = sorted(set(enames))
    max_am = {}
    for elem in species:
        sp = basis_data[elem]  # basis_data is dict with species names as keys
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
        sp = basis_data[name]  # basis_data is dict with species names as keys
        norb = sum(2 * orb['AngularMomentum'] + 1 for orb in sp['orbitals'])
        for orb in sp['orbitals']:
            max_l = max(max_l, orb['AngularMomentum'])
        norb_per_atom.append(norb)
        orb_offsets.append(orb_offsets[-1] + norb)
    return (np.array(norb_per_atom, dtype=np.int32),
            np.array(orb_offsets, dtype=np.int32),
            max_l)


def main():
    global args
    args = parse_args()

    # Extract system name from XYZ path for identification
    system_name = Path(args.xyz).stem
    
    # Create output directory with descriptive name
    if args.output_dir:
        output_dir = Path(args.output_dir)
    else:
        output_dir = Path('.') / f'outputs_{system_name}_{args.basis}_z{args.z_offset}'
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("End-to-End DFTB Dense Projection")
    print("=" * 70)

    # 1. Load geometry
    print(f"\nLoading {args.xyz}")
    coords_ang, _, enames, _, _ = au.load_xyz(args.xyz)
    coords_ang = np.array(coords_ang, dtype=np.float64)
    coords_bohr = coords_ang * ANG2BOHR
    natoms = len(enames)
    print(f"  {natoms} atoms: {set(enames)}")

    # 2. Load basis
    basis_path = REPO_ROOT / 'pyBall/DFTB/data' / f'wfc.{args.basis}.hsd'
    print(f"\nLoading basis from {basis_path}")
    basis_data = parse_wfc_hsd(str(basis_path))
    basis = convert_wfc_to_species_list_ang(basis_data, resolution_bohr=0.04)
    print(f"  Species: {[sp['name'] for sp in basis]}")

    # 3. Create work directory and setup DFTB
    work_dir = tempfile.mkdtemp(prefix='dftb_dense_')
    print(f"\nWork dir: {work_dir}")
    copy_sk_files(work_dir, enames)
    hsd_path = write_dftb_input(work_dir, enames, coords_ang, basis_data)

    # 4. Run DFTBcore SCF
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

    # 5. Build orbital layout
    norb_per_atom, orb_offsets, max_l = build_orbital_layout(basis_data, enames)
    max_shells = 3 if max_l >= 2 else 2
    print(f"\nOrbital layout: max_l={max_l}, max_shells={max_shells}, norb_total={orb_offsets[-1]}")

    # 6. Setup projector
    dftb_data = {
        'coords_bohr': coords_bohr,
        'species_per_atom': [ELEM_Z[n] for n in enames],
        'species_names': enames
    }
    projector, atoms_dict = setup_gridprojector_from_dftb(dftb_data, basis, verbosity=1, max_shells=max_shells)

    # 7. Generate grid
    rmin = coords_ang.min() - 2.0
    rmax = coords_ang.max() + 2.0
    ngrid = int(np.ceil((rmax - rmin) / args.step))
    print(f"\nGrid: {ngrid}x{ngrid} at z={args.z_offset} Å, step={args.step} Å")
    points, extent = generate_2d_point_grid(plane='xy', npoints=ngrid, z_offset=args.z_offset, xy_range=(rmin, rmax))
    points = points.astype(np.float32)

    # 8. Determine occupied orbitals
    n_electrons = sum(VALENCE.get(n, 0) for n in enames)
    nocc = n_electrons // 2
    print(f"  Electrons: {n_electrons}, Occupied MOs: {nocc}")

    # 9. Project density from DM
    print("\nProjecting density from DM...")
    dm = dm.astype(np.float32)
    rho = projector.project_density_dense_points(points, dm, norb_per_atom, orb_offsets, atoms_dict)
    print(f"  rho: min={rho.min():.4e}, max={rho.max():.4e}, sum={rho.sum():.4e}")

    # 10. Project orbitals around HOMO
    homo_idx = nocc - 1
    mo_indices = [max(0, homo_idx - 1), homo_idx, min(len(eigvals) - 1, homo_idx + 1)]
    print(f"\nProjecting orbitals {mo_indices} (HOMO={homo_idx})")
    orb_results = {}
    for imo in mo_indices:
        coeffs = eigvecs[imo].astype(np.float32)
        psi = projector.project_orbital_dense_points(points, coeffs, norb_per_atom, orb_offsets, atoms_dict)
        orb_results[imo] = psi
        print(f"  MO{imo + 1}: E={eigvals[imo] * 27.2114:.4f} eV")

    # 10b. Project specific orbitals for partial density (LUMO+1,+2,+3)
    lumo_offsets = [1, 2, 3]  # HOMO+1, HOMO+2, HOMO+3
    lumo_indices = [min(len(eigvals) - 1, homo_idx + offset) for offset in lumo_offsets]
    print(f"\nProjecting LUMO orbitals {lumo_indices} for partial density (HOMO+1,+2,+3)")
    lumo_results = {}
    for imo in lumo_indices:
        coeffs = eigvecs[imo].astype(np.float32)
        psi = projector.project_orbital_dense_points(points, coeffs, norb_per_atom, orb_offsets, atoms_dict)
        lumo_results[imo] = psi
        print(f"  MO{imo + 1}: E={eigvals[imo] * 27.2114:.4f} eV")

    # Sum squares of LUMO orbitals for partial density
    rho_lumo = np.zeros_like(rho)
    for imo in lumo_indices:
        rho_lumo += lumo_results[imo] ** 2
    print(f"  LUMO partial density: min={rho_lumo.min():.4e}, max={rho_lumo.max():.4e}, sum={rho_lumo.sum():.4e}")

    # 11. Plot
    if not args.no_plot:
        print("\n" + "=" * 70)
        print("Generating Plots")
        print("=" * 70)

        # Density plot
        rho_grid = rho.reshape(ngrid, ngrid)
        fig, ax = plt.subplots(figsize=(8, 6))
        im = ax.imshow(rho_grid, origin='lower', extent=extent, cmap='viridis')
        ax.set_title(f'Density from DM\n{system_name} | {args.basis} | z={args.z_offset} Å')
        ax.set_xlabel('X (Å)')
        ax.set_ylabel('Y (Å)')
        plt.colorbar(im, ax=ax)
        ax.scatter(coords_ang[:, 0], coords_ang[:, 1], c='black', marker='o', s=15, alpha=0.5, zorder=10)
        plt.tight_layout()
        density_file = output_dir / f'density_{system_name}_{args.basis}_z{args.z_offset}.png'
        plt.savefig(density_file, dpi=args.dpi)
        print(f"  Saved: {density_file}")
        plt.close()

        # LUMO partial density plot
        rho_lumo_grid = rho_lumo.reshape(ngrid, ngrid)
        fig, ax = plt.subplots(figsize=(8, 6))
        im = ax.imshow(rho_lumo_grid, origin='lower', extent=extent, cmap='viridis')
        ax.set_title(f'LUMO Partial Density (HOMO+1,+2,+3)\n{system_name} | {args.basis} | z={args.z_offset} Å')
        ax.set_xlabel('X (Å)')
        ax.set_ylabel('Y (Å)')
        plt.colorbar(im, ax=ax)
        ax.scatter(coords_ang[:, 0], coords_ang[:, 1], c='black', marker='o', s=15, alpha=0.5, zorder=10)
        plt.tight_layout()
        lumo_density_file = output_dir / f'density_lumo_{system_name}_{args.basis}_z{args.z_offset}.png'
        plt.savefig(lumo_density_file, dpi=args.dpi)
        print(f"  Saved: {lumo_density_file}")
        plt.close()

        # Orbital plots
        fig, axes = plt.subplots(1, len(mo_indices), figsize=(5 * len(mo_indices), 4))
        if len(mo_indices) == 1:
            axes = [axes]
        for i, imo in enumerate(mo_indices):
            psi_grid = orb_results[imo].reshape(ngrid, ngrid)
            vmax = np.max(np.abs(psi_grid))
            im = axes[i].imshow(psi_grid, origin='lower', extent=extent, cmap='RdBu_r', vmin=-vmax, vmax=vmax)
            # Calculate relative index to HOMO
            rel_idx = imo - homo_idx
            rel_label = f'HOMO{rel_idx:+d}' if rel_idx != 0 else 'HOMO'
            axes[i].set_title(f'MO{imo + 1} ({rel_label})\nE={eigvals[imo] * 27.2114:.2f} eV')
            axes[i].set_xlabel('X (Å)')
            axes[i].set_ylabel('Y (Å)')
            plt.colorbar(im, ax=axes[i])
            axes[i].scatter(coords_ang[:, 0], coords_ang[:, 1], c='black', marker='o', s=15, alpha=0.5, zorder=10)
        plt.suptitle(f'Orbital Projection\n{system_name} | {args.basis} | z={args.z_offset} Å', fontsize=14)
        plt.tight_layout()
        orbitals_file = output_dir / f'orbitals_{system_name}_{args.basis}_z{args.z_offset}.png'
        plt.savefig(orbitals_file, dpi=args.dpi)
        print(f"  Saved: {orbitals_file}")
        plt.close()

    # Cleanup
    shutil.rmtree(work_dir)
    print(f"\nCleaned up work dir")
    print("=" * 70)
    print("Done!")
    print("=" * 70)
    return 0


if __name__ == '__main__':
    sys.exit(main())
