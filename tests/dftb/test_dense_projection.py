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
    p.add_argument('--dpi', type=int, default=100, help='Plot DPI')
    p.add_argument('--use-exp-basis', action='store_true', help='Use exponential radial decay for STM (LUMO orbitals only)')
    p.add_argument('--exp-beta', type=float, default=1.0, help='Exponential decay constant (Å^-1), higher = steeper decay')
    p.add_argument('--exp-r0', type=float, default=3.0, help='Reference distance (Å) where f=1')
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

    # 7. Generate grid for density (always at z=2.0 with original basis)
    rmin = coords_ang.min() - 2.0
    rmax = coords_ang.max() + 2.0
    ngrid = int(np.ceil((rmax - rmin) / args.step))
    z_offset_density = 2.0  # Always use z=2.0 for density projection
    print(f"\nDensity grid: {ngrid}x{ngrid} at z={z_offset_density} Å, step={args.step} Å")
    points, extent = generate_2d_point_grid(plane='xy', npoints=ngrid, z_offset=z_offset_density, xy_range=(rmin, rmax))
    points = points.astype(np.float32)

    # 7b. Generate grid for STM (LUMO) at user-specified z-offset
    ngrid_stm = int(np.ceil((rmax - rmin) / args.step))
    print(f"\nSTM grid: {ngrid_stm}x{ngrid_stm} at z={args.z_offset} Å, step={args.step} Å")
    points_stm, extent_stm = generate_2d_point_grid(plane='xy', npoints=ngrid_stm, z_offset=args.z_offset, xy_range=(rmin, rmax))
    points_stm = points_stm.astype(np.float32)

    # 8. Determine occupied orbitals
    n_electrons = sum(VALENCE.get(n, 0) for n in enames)
    nocc = n_electrons // 2
    print(f"  Electrons: {n_electrons}, Occupied MOs: {nocc}")

    # 9. Project density from DM
    print("\nProjecting density from DM...")
    dm = dm.astype(np.float32)
    rho = projector.project_density_dense_points(points, dm, norb_per_atom, orb_offsets, atoms_dict)
    print(f"  rho: min={rho.min():.4e}, max={rho.max():.4e}, sum={rho.sum():.4e}")

    # 10. Project orbitals around HOMO (always use original basis at z=2.0)
    homo_idx = nocc - 1
    mo_indices = [max(0, homo_idx - 1), homo_idx, min(len(eigvals) - 1, homo_idx + 1)]
    print(f"\nProjecting orbitals {mo_indices} (HOMO={homo_idx}) at z={z_offset_density} Å")
    orb_results = {}
    for imo in mo_indices:
        coeffs = eigvecs[imo].astype(np.float32)
        psi = projector.project_orbital_dense_points(points, coeffs, norb_per_atom, orb_offsets, atoms_dict)
        orb_results[imo] = psi
        print(f"  MO{imo + 1}: E={eigvals[imo] * 27.2114:.4f} eV")

    # 10b. Project LUMO orbitals for STM (use exponential if requested)
    lumo_offsets = [1, 2, 3]  # HOMO+1, HOMO+2, HOMO+3
    lumo_indices = [min(len(eigvals) - 1, homo_idx + offset) for offset in lumo_offsets]
    print(f"\nProjecting LUMO orbitals {lumo_indices} for STM (HOMO+1,+2,+3) at z={args.z_offset} Å")
    if args.use_exp_basis:
        print(f"  Using exponential radial decay: beta={args.exp_beta} Å^-1, r0={args.exp_r0} Å")
    lumo_results = {}
    for imo in lumo_indices:
        coeffs = eigvecs[imo].astype(np.float32)
        if args.use_exp_basis:
            psi = projector.project_orbital_dense_points_exp(points_stm, coeffs, norb_per_atom, orb_offsets, atoms_dict, beta=args.exp_beta, r0=args.exp_r0)
        else:
            psi = projector.project_orbital_dense_points(points_stm, coeffs, norb_per_atom, orb_offsets, atoms_dict)
        lumo_results[imo] = psi
        print(f"  MO{imo + 1}: E={eigvals[imo] * 27.2114:.4f} eV")

    # Sum squares of LUMO orbitals for partial density
    rho_lumo = np.zeros_like(points_stm[:, 0], dtype=np.float32)
    for imo in lumo_indices:
        rho_lumo += lumo_results[imo] ** 2
    print(f"  LUMO partial density: min={rho_lumo.min():.4e}, max={rho_lumo.max():.4e}, sum={rho_lumo.sum():.4e}")

    # 11. Plot
    if not args.no_plot:
        print("\n" + "=" * 70)
        print("Generating Plots")
        print("=" * 70)

        # Density plot (always at z=2.0 with original basis)
        rho_grid = rho.reshape(ngrid, ngrid)
        fig, ax = plt.subplots(figsize=(8, 6))
        im = ax.imshow(rho_grid, origin='lower', extent=extent, cmap='viridis')
        ax.set_title(f"Electron Density - {system_name} ({args.basis}) z={z_offset_density:.1f} Å")
        ax.set_xlabel("x (Å)")
        ax.set_ylabel("y (Å)")
        plt.colorbar(im, ax=ax, label="rho (e/Å³)")
        density_path = os.path.join(output_dir, f"density_{system_name}_{args.basis}_z{z_offset_density:.1f}.png")
        plt.savefig(density_path, dpi=args.dpi)
        plt.close()
        print(f"  Saved: {density_path}")

        # LUMO partial density plot (at STM z-offset, exponential if requested)
        rho_lumo_grid = rho_lumo.reshape(ngrid_stm, ngrid_stm)
        fig, ax = plt.subplots(figsize=(8, 6))
        im = ax.imshow(rho_lumo_grid, origin='lower', extent=extent_stm, cmap='viridis')
        mode_str = "exp" if args.use_exp_basis else "basis"
        beta_str = f"beta={args.exp_beta:.1f}" if args.use_exp_basis else ""
        title = f"STM (LUMO HOMO+1,+2,+3) - {system_name} ({args.basis}) z={args.z_offset:.1f} Å ({mode_str}"
        if beta_str:
            title += f", {beta_str}"
        title += ")"
        ax.set_title(title)
        ax.set_xlabel("x (Å)")
        ax.set_ylabel("y (Å)")
        plt.colorbar(im, ax=ax, label="rho_LUMO (e/Å³)")
        stm_filename = f"STM_{system_name}_{args.basis}_z{args.z_offset:.1f}"
        if args.use_exp_basis:
            stm_filename += f"_beta{args.exp_beta:.1f}"
        stm_path = os.path.join(output_dir, f"{stm_filename}.png")
        plt.savefig(stm_path, dpi=args.dpi)
        plt.close()
        print(f"  Saved: {stm_path}")

        # Orbitals plot (always at z=2.0 with original basis)
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        for i, imo in enumerate(mo_indices):
            psi_grid = orb_results[imo].reshape(ngrid, ngrid)
            im = axes[i].imshow(psi_grid, origin='lower', extent=extent, cmap='RdBu_r')
            rel_idx = imo - homo_idx
            label = "HOMO" if rel_idx == 0 else f"HOMO{rel_idx:+d}"
            axes[i].set_title(f"MO{imo + 1} ({label})\nE={eigvals[imo] * 27.2114:.4f} eV")
            axes[i].set_xlabel("x (Å)")
            axes[i].set_ylabel("y (Å)")
            plt.colorbar(im, ax=axes[i])
        plt.tight_layout()
        orb_path = os.path.join(output_dir, f"orbitals_{system_name}_{args.basis}_z{z_offset_density:.1f}.png")
        plt.savefig(orb_path, dpi=args.dpi)
        plt.close()
        print(f"  Saved: {orb_path}")

    # Cleanup
    shutil.rmtree(work_dir)
    print(f"\nCleaned up work dir")
    print("=" * 70)
    print("Done!")
    print("=" * 70)
    return 0


if __name__ == '__main__':
    sys.exit(main())
