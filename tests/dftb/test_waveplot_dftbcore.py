#!/usr/bin/env python3
"""
test_waveplot_dftbcore.py

Test orbital projection using eigenvectors exported DIRECTLY from the DFTBcore
library (libdftbcore.so) — no eigenvec.bin file required.

The eigenvector matrix from the library is in (nStates, nOrb) C row-major order,
identical to what parse_eigenvec_bin_custom() produces from eigenvec.bin:
  - Fortran stores eigvecsReal(nOrb, nStates) column-major
  - .bin file is the raw Fortran memory dump → same column-major layout
  - DFTBcore.get_eigvecs_dense() does: reshape(n,n, order='F').T  → (nStates, nOrb)
  - parse_eigenvec_bin_custom() also produces (nStates, nOrb)
So the two paths are equivalent; no extra reordering is needed.

Optionally compares against eigenvec.bin (if present) to validate.

Usage:
    # H2O, 2D XY plane, compare lib vs .bin
    python test_waveplot_dftbcore.py --dftb-dir tests/grid/dftb_h2o --points --plane2d xy --z-offset 0.0

    # PTCDA HOMO-4..LUMO+4, XY plane at z=2 Å
    python test_waveplot_dftbcore.py --dftb-dir tests/grid/dftb_ptcda --points --plane2d xy --z-offset 2.0 --mo-range 66 75 --npoints 64

    # H2O, 3D grid projection
    python test_waveplot_dftbcore.py --dftb-dir tests/grid/dftb_h2o

    # 1D z-scan
    python test_waveplot_dftbcore.py --dftb-dir tests/grid/dftb_h2o --points
"""

import sys, os, argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

REPO_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

from pyBall.DFTB.DFTBplusParser import ( parse_basis_hsd_ang, parse_detailed_xml_custom, parse_eigenvec_bin_custom, build_wp_basis, evec_to_kernel_coeffs)
from pyBall.DFTB.Grid_dftb import GridProjector, setup_gridprojector_from_dftb, evaluate_mos_on_points as ocl_evaluate_mos_on_points
from pyBall.DFTB.TestUtils import print_eigenvecs
from pyBall.DFTB.DFTBcore import DFTBcore

BOHR2ANG = 0.5291772109
OUTPUT_DIR = Path(__file__).parent / 'waveplot_output' / 'dftbcore'


# ================================================================
# CLI
# ================================================================

def parse_args():
    p = argparse.ArgumentParser( description='DFTBcore orbital projection test (no eigenvec.bin needed)', formatter_class=argparse.ArgumentDefaultsHelpFormatter )
    p.add_argument('--dftb-dir',  type=str, default=None, help='DFTB+ run dir containing dftb_in.hsd, waveplot_in.hsd, detailed.xml')
    p.add_argument('--lib-path',  type=str, default=None,  help='Path to libdftbcore.so (auto-detected if absent)')
    p.add_argument('--no-show',   action='store_true')
    p.add_argument('--dpi',       type=int, default=150)
    p.add_argument('--nmo',       type=int, default=6,   help='Number of MOs around HOMO (default range if --mo-range not set)')
    p.add_argument('--mo-range',  type=int, nargs=2, default=None, metavar=('START','END'),   help='1-based inclusive MO index range')
    p.add_argument('--compare-bin', action='store_true',  help='Compare library eigenvectors against eigenvec.bin (validation)')
    # Point evaluation
    p.add_argument('--points',    action='store_true',  help='Evaluate at explicit points instead of full 3D grid')
    p.add_argument('--plane2d',   type=str, choices=['xy','xz','yz'], default=None,  help='2D plane for --points; omit for 1D z-scan')
    p.add_argument('--z-offset',  type=float, default=0.0,  help='Fixed coordinate value for out-of-plane axis (Å)')
    p.add_argument('--xy-range',  type=float, nargs=2, default=None, metavar=('MIN','MAX'),  help='Coordinate range for 2D scan (Å); default: mol extent + 3 Å')
    p.add_argument('--z-range',   type=float, nargs=2, default=[-3.0, 3.0])
    p.add_argument('--npoints',   type=int, default=64)
    # 3D grid
    p.add_argument('--step',      type=float, default=0.2,  help='Grid spacing in Å (3D grid mode)')
    p.add_argument('--margin',    type=float, default=3.0,  help='Grid margin around molecule in Å (3D grid mode)')
    p.add_argument('--nmax-atom', type=int, default=64)
    p.add_argument('--print-eigenvec', action='store_true',  help='Print eigenvectors from eigenvec.bin and exit')
    return p.parse_args()


# ================================================================
# Core: run DFTB+ via library and extract eigenvectors
# ================================================================

def run_dftb_and_get_eigvecs(dftb_dir, lib_path=None):
    """
    Run DFTB+ via libdftbcore.so and return eigenvectors (nStates, nOrb).

    DFTB+ always reads 'dftb_in.hsd' from the current working directory
    (hardcoded in hsdhelpers.F90:parseHsdInput). We chdir into dftb_dir,
    run, then restore the original CWD.

    The Fortran array eigvecsReal(nOrb, nStates) is stored column-major.
    dftbcore_get_eigvecs_dense() flattens it as a Fortran-contiguous buffer
    and DFTBcore.get_eigvecs_dense() reshapes it as:
        np.asfortranarray(buf.reshape(n, n, order='F')).T.copy()
    which gives (nStates, nOrb) in C order — same as parse_eigenvec_bin_custom().
    """
    assert (dftb_dir / 'dftb_in.hsd').exists(), f"dftb_in.hsd not found in {dftb_dir}"

    orig_cwd = os.getcwd()
    try:
        os.chdir(dftb_dir)
        dftb = DFTBcore(libpath=lib_path)
        dftb.init('dftb_in.hsd')
        dftb.enable_matrix_collection(dm=False, h=False, s=False)
        energy = dftb.run_scf()
        eigvecs, eigenvals = dftb.get_eigvecs_dense()   # (nStates, nOrb), (nStates,)
        dftb.finalize()
    finally:
        os.chdir(orig_cwd)

    print(f"  DFTBcore: E={energy:.6f} Ha  nOrb={eigvecs.shape[1]}  nStates={eigvecs.shape[0]}")
    return eigvecs, eigenvals


# ================================================================
# Optional validation: compare library eigenvecs against .bin
# ================================================================

def compare_against_bin(evecs_lib, eigenvals_lib, dftb_dir, nstates, norb):
    """
    Compare library eigenvectors against eigenvec.bin.
    Eigenvectors may differ by a global sign per column (degenerate states).
    Reports max|diff| after sign alignment.
    """
    bin_path = dftb_dir / 'eigenvec.bin'
    if not bin_path.exists():
        print("  [compare_bin] eigenvec.bin not found — skipping")
        return
    evecs_bin = parse_eigenvec_bin_custom(bin_path, nstates, norb)
    # Sign-align: flip sign of library MO if dot product with .bin MO is negative
    for i in range(nstates):
        if np.dot(evecs_lib[i], evecs_bin[i]) < 0:
            evecs_lib[i] *= -1
    diff = np.abs(evecs_lib - evecs_bin)
    print(f"  [compare_bin] max|lib - bin| after sign alignment = {diff.max():.3e}")
    print(f"  [compare_bin] mean|lib - bin| = {diff.mean():.3e}")


# ================================================================
# Main
# ================================================================

def main():
    args = parse_args()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    dftb_dir = Path(args.dftb_dir) if args.dftb_dir else Path(__file__).parent / 'dftb_h2o'
    assert dftb_dir.exists(), f"DFTB+ directory not found: {dftb_dir}"
    system_name = dftb_dir.name
    
    # Print eigenvectors if requested
    if args.print_eigenvec:
        print_eigenvecs(dftb_dir / 'eigenvec.bin', dftb_dir / 'detailed.xml', dftb_dir / 'waveplot_in.hsd', max_orbitals=args.nmo if args.nmo else None)
        return
    
    print("=" * 60)
    print(f"test_waveplot_dftbcore.py — {system_name}")
    print(f"Directory: {dftb_dir}")
    print("=" * 60)

    # ---- 1. Parse detailed.xml for geometry / occupation info ----
    geo = parse_detailed_xml_custom(dftb_dir / 'detailed.xml')
    atom_coords_b    = geo['coords_bohr']          # (natoms, 3) Bohr
    species_per_atom = geo['species_per_atom']     # 0-based
    species_names    = geo['species_names']
    natoms           = geo['natoms']
    nstates_total    = geo['nstates']
    norb_total       = geo['norb']
    occs_full        = geo['occupations']
    occs_flat        = occs_full[:, 0, 0]
    print(f"  natoms={natoms}  species={species_names}  nstates={nstates_total}  norb={norb_total}")

    # ---- 2. Parse STO basis from waveplot_in.hsd ----
    hsd_path = dftb_dir / 'waveplot_in.hsd'
    assert hsd_path.exists(), f"waveplot_in.hsd not found: {hsd_path}"
    species_list_ang = parse_basis_hsd_ang(hsd_path)
    sp_by_name = {sp['name']: sp for sp in species_list_ang}
    norb_per_atom = np.array(
        [sum(2*o['l']+1 for o in sp_by_name[species_names[si]]['orbitals'])
         for si in species_per_atom], dtype=np.int32
    )
    assert norb_per_atom.sum() == norb_total, \
        f"norb mismatch: basis sums to {norb_per_atom.sum()}, xml says {norb_total}"
    print(f"  STO basis loaded: {[sp['name'] for sp in species_list_ang]}")

    # ---- 3. Get eigenvectors from DFTBcore library (not from .bin) ----
    print("\n[DFTBcore] Running DFTB+ via library...")
    evecs_full, eigenvals_lib = run_dftb_and_get_eigvecs(dftb_dir, lib_path=args.lib_path)
    assert evecs_full.shape == (nstates_total, norb_total), \
        f"Eigenvec shape mismatch: got {evecs_full.shape}, expected ({nstates_total},{norb_total})"

    # ---- 4. Optional: validate against eigenvec.bin ----
    if args.compare_bin:
        print("\n[Validation] Comparing library eigvecs against eigenvec.bin...")
        compare_against_bin(evecs_full.copy(), eigenvals_lib, dftb_dir, nstates_total, norb_total)

    # ---- HOMO ----
    homo_idx = np.where(occs_flat > 0.5)[0]
    homo = int(homo_idx[-1]) + 1 if len(homo_idx) else nstates_total // 2

    # ---- Energies from band.out (labels only) ----
    energies_ev = np.zeros(nstates_total)
    band_path = dftb_dir / 'band.out'
    if band_path.exists():
        with open(band_path) as f:
            for line in f:
                parts = line.split()
                if len(parts) == 3:
                    try:
                        idx = int(parts[0])
                        if 1 <= idx <= nstates_total:
                            energies_ev[idx-1] = float(parts[1])
                    except ValueError:
                        pass
    else:
        # Fall back to eigenvalues from library (Hartree -> eV)
        HARTREE2EV = 27.2114
        energies_ev = eigenvals_lib * HARTREE2EV

    # ---- MO selection ----
    if args.mo_range:
        mo_start = max(1, args.mo_range[0])
        mo_end   = min(nstates_total, args.mo_range[1])
    else:
        mo_start = max(1, homo - args.nmo // 2)
        mo_end   = min(nstates_total, homo + args.nmo // 2)
    nstates  = mo_end - mo_start + 1
    evecs    = evecs_full[mo_start-1 : mo_end]     # (nstates, norb) for selected MOs
    energies = energies_ev[mo_start-1 : mo_end]
    occs     = occs_flat[mo_start-1 : mo_end]

    print(f"\n  HOMO = MO{homo} ({energies_ev[homo-1]:.3f} eV)")
    print(f"  Evaluating MOs {mo_start}–{mo_end}  ({nstates} total)")
    for i in range(min(nstates, 15)):
        tag = " <--HOMO" if (mo_start+i)==homo else (" <--LUMO" if (mo_start+i)==homo+1 else "")
        print(f"    MO{mo_start+i:4d}  E={energies[i]:8.3f} eV  occ={occs[i]:.1f}{tag}")

    # ---- Setup OpenCL projector ----
    dftb_data_ocl = {
        'coords_bohr':    atom_coords_b,
        'species_per_atom': species_per_atom,
        'species_names':  species_names,
    }
    projector, atoms_dict = setup_gridprojector_from_dftb(dftb_data_ocl, species_list_ang, verbosity=0)

    mo_indices = list(range(nstates))   # 0-indexed within the selected slice

    # ================================================================
    # Method 2: Point evaluation
    # ================================================================
    if args.points:
        atom_ang = atom_coords_b * BOHR2ANG
        if args.xy_range:
            rmin, rmax_r = args.xy_range
        else:
            rmin   = float(atom_ang.min()) - 3.0
            rmax_r = float(atom_ang.max()) + 3.0

        u = np.linspace(rmin, rmax_r, args.npoints)
        if args.plane2d == 'xy':
            uu, vv = np.meshgrid(u, u, indexing='ij')
            ww = np.full_like(uu, args.z_offset)
            points_ang = np.column_stack([uu.ravel(), vv.ravel(), ww.ravel()])
            ax_labels = ('x (Å)', 'y (Å)')
            plane_desc = f"XY plane  z={args.z_offset:.3f} Å"
        elif args.plane2d == 'xz':
            uu, vv = np.meshgrid(u, u, indexing='ij')
            points_ang = np.column_stack([uu.ravel(), np.full(uu.size, args.z_offset), vv.ravel()])
            ax_labels = ('x (Å)', 'z (Å)')
            plane_desc = f"XZ plane  y={args.z_offset:.3f} Å"
        elif args.plane2d == 'yz':
            uu, vv = np.meshgrid(u, u, indexing='ij')
            points_ang = np.column_stack([np.full(uu.size, args.z_offset), uu.ravel(), vv.ravel()])
            ax_labels = ('y (Å)', 'z (Å)')
            plane_desc = f"YZ plane  x={args.z_offset:.3f} Å"
        else:
            z_vals = np.linspace(args.z_range[0], args.z_range[1], args.npoints)
            points_ang = np.column_stack([np.zeros(args.npoints), np.zeros(args.npoints), z_vals])
            ax_labels = ('z (Å)', 'ψ')
            plane_desc = "1D z-scan  x=0  y=0"

        print(f"\n[Points] {len(points_ang)} pts  {plane_desc}")

        ocl_vals_list = ocl_evaluate_mos_on_points(
            projector, mo_indices, points_ang.astype(np.float32),
            evecs, natoms, species_per_atom, species_names, species_list_ang,
            norb_per_atom, atoms_dict
        )
        ocl_vals = np.array(ocl_vals_list)
        print(f"  OpenCL done. max|ψ|={np.abs(ocl_vals).max():.4e}")

        mo_labels = list(range(mo_start, mo_end + 1))

        if args.plane2d:
            s2 = (args.npoints, args.npoints)
            n_cols = min(4, nstates)
            n_rows = (nstates + n_cols - 1) // n_cols
            fig, axes = plt.subplots(n_rows, n_cols, figsize=(4.5*n_cols, 4*n_rows))
            axes = np.array(axes).flatten()
            extent = [rmin, rmax_r, rmin, rmax_r]
            for idx, imo_abs in enumerate(mo_labels):
                ax = axes[idx]
                dat = ocl_vals[idx].reshape(s2)
                clim = max(abs(dat.min()), abs(dat.max())) or 1e-10
                im = ax.imshow(dat.T, origin='lower', extent=extent,
                               cmap='RdBu_r', vmin=-clim, vmax=clim)
                tag = " HOMO" if imo_abs==homo else (" LUMO" if imo_abs==homo+1 else "")
                ax.set_title(f"MO{imo_abs}{tag}\nE={energies[idx]:.2f}eV  occ={occs[idx]:.0f}",
                             fontsize=8)
                ax.set_xlabel(ax_labels[0]); ax.set_ylabel(ax_labels[1])
                plt.colorbar(im, ax=ax, fraction=0.046)
            for ax in axes[nstates:]:
                ax.set_visible(False)
            fig.suptitle(f"{system_name} [DFTBcore lib]  {plane_desc}", fontsize=10)
            fig.tight_layout()
            out = OUTPUT_DIR / f"{system_name}_{args.plane2d}_z{args.z_offset:.2f}_MO{mo_start}-{mo_end}.png"
            fig.savefig(str(out), dpi=args.dpi)
            print(f"\nSaved: {out}")
            if not args.no_show: plt.show()
        else:
            # 1D z-scan
            n_cols = min(3, nstates)
            n_rows = (nstates + n_cols - 1) // n_cols
            fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 3.5*n_rows))
            axes = np.array(axes).flatten()
            for idx, imo_abs in enumerate(mo_labels):
                ax = axes[idx]
                tag = " HOMO" if imo_abs==homo else (" LUMO" if imo_abs==homo+1 else "")
                ax.plot(z_vals, ocl_vals[idx], 'r-', lw=1.5, label='DFTBcore+OCL')
                ax.axhline(0, c='gray', lw=0.5, ls='--')
                ax.set_xlabel('z (Å)'); ax.set_ylabel('ψ')
                ax.set_title(f"MO{imo_abs}{tag}  E={energies[idx]:.2f}eV", fontsize=8)
                ax.legend(fontsize=7)
            for ax in axes[nstates:]:
                ax.set_visible(False)
            fig.suptitle(f"{system_name} [DFTBcore lib]  {plane_desc}", fontsize=10)
            fig.tight_layout()
            out = OUTPUT_DIR / f"{system_name}_1dscan_MO{mo_start}-{mo_end}.png"
            fig.savefig(str(out), dpi=args.dpi)
            print(f"\nSaved: {out}")
            if not args.no_show: plt.show()
        return

    # ================================================================
    # Method 1: 3D grid projection
    # ================================================================
    coords_ang = atom_coords_b * BOHR2ANG
    pos_min = coords_ang.min(axis=0) - args.margin
    pos_max = coords_ang.max(axis=0) + args.margin
    span    = pos_max - pos_min
    block   = 8
    ngrid   = ((np.ceil(span / args.step).astype(int) + block - 1) // block) * block
    origin  = pos_min
    grid_spec = {
        'origin': origin,
        'dA': np.array([args.step, 0., 0.]),
        'dB': np.array([0., args.step, 0.]),
        'dC': np.array([0., 0., args.step]),
        'ngrid': ngrid,
    }
    print(f"\n[3D grid]  ngrid={ngrid}  step={args.step} Å  margin={args.margin} Å")

    grids = []
    for i, imo_abs in enumerate(range(mo_start, mo_end+1)):
        coeffs_k = evec_to_kernel_coeffs(
            evecs[i], natoms, species_per_atom, species_names, species_list_ang
        )
        grid_3d = projector.project_orbital(coeffs_k, norb_per_atom, atoms_dict, grid_spec,
                                            nMaxAtom=args.nmax_atom)
        grids.append(grid_3d)
        print(f"  MO{imo_abs}  E={energies[i]:.2f}eV  occ={occs[i]:.0f}  |ψ|max={np.abs(grid_3d).max():.4e}")

    # Plot middle XY slice for each MO
    n_cols = min(4, nstates)
    n_rows = (nstates + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4.5*n_cols, 4*n_rows))
    axes = np.array(axes).flatten()
    for i, imo_abs in enumerate(range(mo_start, mo_end+1)):
        ax    = axes[i]
        g3d   = grids[i]
        iz    = g3d.shape[2] // 2
        g2d   = g3d[:, :, iz]
        z_val = origin[2] + iz * args.step
        ext   = [origin[0], origin[0]+ngrid[0]*args.step,
                 origin[1], origin[1]+ngrid[1]*args.step]
        clim  = max(abs(g2d.min()), abs(g2d.max())) or 1e-10
        im    = ax.imshow(g2d.T, origin='lower', extent=ext, cmap='RdBu_r', vmin=-clim, vmax=clim)
        tag   = " HOMO" if imo_abs==homo else (" LUMO" if imo_abs==homo+1 else "")
        ax.set_title(f"MO{imo_abs}{tag}  E={energies[i]:.2f}eV\nz={z_val:.2f}Å  occ={occs[i]:.0f}", fontsize=8)
        ax.set_xlabel('x (Å)'); ax.set_ylabel('y (Å)')
        plt.colorbar(im, ax=ax, fraction=0.046)
    for ax in axes[nstates:]:
        ax.set_visible(False)
    fig.suptitle(f"{system_name} [DFTBcore lib] — XY slice  (MO{mo_start}–{mo_end})", fontsize=10)
    fig.tight_layout()
    out = OUTPUT_DIR / f"{system_name}_grid3d_MO{mo_start}-{mo_end}.png"
    fig.savefig(str(out), dpi=args.dpi)
    print(f"\nSaved: {out}")
    if not args.no_show: plt.show()


if __name__ == '__main__':
    main()
