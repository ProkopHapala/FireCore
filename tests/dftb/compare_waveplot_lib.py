#!/usr/bin/env python3
"""
Compare libwaveplot.so orbital grids vs reference WAVEPLOT cube files, and vs OpenCL.

General-purpose: works for any molecule given a DFTB+ run directory.
All system data (geometry, eigenvectors, basis) is read automatically from:
  <dftb-dir>/detailed.xml   -- geometry, nstates, norb, occupations
  <dftb-dir>/eigenvec.bin   -- eigenvectors (Fortran unformatted binary)
  <dftb-dir>/waveplot_in.hsd -- STO basis parameters

Usage:
  # Method 1: 3D grid comparison (needs wp-1-1-N-real.cube files)
  python compare_waveplot_lib.py --dftb-dir tests/grid/dftb_h2o --nmo 6

  # Method 2: 2D point evaluation on XY plane at z=2 Å
  python compare_waveplot_lib.py --dftb-dir tests/grid/dftb_h2o --points --plane2d xy --z-offset 2.0 --npoints 51

  # PTCDA, HOMO-4 to LUMO+4
  python compare_waveplot_lib.py --dftb-dir tests/grid/dftb_ptcda --points --plane2d xy --z-offset 2.0 --mo-range 66 75 --npoints 64
"""

import sys, os, struct, argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
REPO_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

from pyBall.DFTB.WavePlot import WavePlot, setup_waveplot_from_dftb, evaluate_mos_on_points as wp_evaluate_mos_on_points
from pyBall.DFTB.DFTBplusParser import ( parse_basis_hsd_ang, parse_detailed_xml_custom, parse_eigenvec_bin_custom, read_cube, build_wp_basis,evec_to_kernel_coeffs)
from pyBall.DFTB.Grid_dftb import GridProjector, setup_gridprojector_from_dftb, evaluate_mos_on_points as ocl_evaluate_mos_on_points
from pyBall.DFTB.TestUtils import ( compare_point_evaluations, print_comparison_results, generate_1d_z_scan, print_eigenvecs, generate_2d_point_grid)
from pyBall.plotUtils import plot_comparison_2d

LIB_PATH   = str(Path('/home/prokop/git/dftbplus/_build/app/waveplot/libwaveplot.so'))
OUTPUT_DIR = Path(__file__).parent / 'waveplot_output' / 'comparison'
BOHR2ANG   = 0.5291772109


# ================================================================
# CLI
# ================================================================

def parse_args():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('--dftb-dir',  type=str, default=None, help='DFTB+ run directory containing detailed.xml, eigenvec.bin, waveplot_in.hsd')
    p.add_argument('--no-show',   action='store_true')
    p.add_argument('--dpi',       type=int, default=150)
    p.add_argument('--nmo',       type=int, default=6, help='Number of MOs (Method 1, or default range if --mo-range not set)')
    p.add_argument('--mo-range',  type=int, nargs=2, default=None, metavar=('START','END'), help='1-based inclusive MO index range')
    p.add_argument('--points',    action='store_true', help='Method 2: evaluate at explicit points (libwaveplot + OpenCL)')
    p.add_argument('--plane2d',   type=str, choices=['xy','xz','yz'], default=None, help='2D plane for --points; omit for 1D scan')
    p.add_argument('--line-scan', type=str, choices=['z','bond'], default=None, help='1D line scan: z (along z-axis at origin) or bond (through two atoms)')
    p.add_argument('--atoms',     type=int, nargs=2, default=None, metavar=('I','J'), help='Atom indices (0-based) for bond line scan')
    p.add_argument('--z-offset',  type=float, default=0.0,  help='Fixed coordinate value for out-of-plane axis (Å)')
    p.add_argument('--xy-range',  type=float, nargs=2, default=None, metavar=('MIN','MAX'), help='Coordinate range for 2D plane (Å); default: mol extent + 3 Å')
    p.add_argument('--z-range',   type=float, nargs=2, default=[-3.0, 3.0],  metavar=('ZMIN','ZMAX'), help='Range for 1D z-scan (Å)')
    p.add_argument('--line-range', type=float, nargs=2, default=None, metavar=('RMIN','RMAX'), help='Range for bond line scan (Å); default: bond length + 3 Å each side')
    p.add_argument('--npoints',   type=int, default=64, help='Grid points per axis')
    p.add_argument('--print-eigenvec', action='store_true', help='Print eigenvectors from eigenvec.bin and exit')
    return p.parse_args()


# ================================================================
# MAIN
# ================================================================

def main():
    args = parse_args()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    dftb_dir = Path(args.dftb_dir) if args.dftb_dir else Path(__file__).parent / 'dftb_h2o'
    assert dftb_dir.exists(), f"Not found: {dftb_dir}"
    system_name = dftb_dir.name
    
    # Print eigenvectors if requested
    if args.print_eigenvec:
        print_eigenvecs(dftb_dir / 'eigenvec.bin', dftb_dir / 'detailed.xml', dftb_dir / 'waveplot_in.hsd', max_orbitals=args.nmo if args.nmo else None)
        return
    
    print("=" * 60)
    print(f"System: {system_name}  ({dftb_dir})")
    print("=" * 60)

    # ---- 1. Parse detailed.xml ----
    geo = parse_detailed_xml_custom(dftb_dir / 'detailed.xml')
    atom_coords_b   = geo['coords_bohr']        # (natoms,3) Bohr
    species_per_atom = geo['species_per_atom']  # 0-based
    species_names    = geo['species_names']
    natoms  = geo['natoms']
    nstates_total = geo['nstates']
    norb_total    = geo['norb']
    nkpoints      = geo['nkpoints']
    nspin         = geo['nspin']
    occs_full     = geo['occupations']          # (nstates,nkpoints,nspin)
    occs_flat     = occs_full[:, 0, 0]          # k=1, spin=1
    print(f"  natoms={natoms}, species={species_names}, nstates={nstates_total}, norb={norb_total}")

    # ---- 2. Parse STO basis ----
    hsd_path = dftb_dir / 'waveplot_in.hsd'
    assert hsd_path.exists(), f"waveplot_in.hsd not found: {hsd_path}"
    species_list_ang = parse_basis_hsd_ang(hsd_path)
    sp_names_hsd = [sp['name'] for sp in species_list_ang]
    print(f"  HSD species: {sp_names_hsd}")

    # species array for libwaveplot: 1-based index into sp_names_hsd
    sp_name_to_idx = {name: i+1 for i, name in enumerate(sp_names_hsd)}
    species_wp = np.array([sp_name_to_idx[species_names[si]] for si in species_per_atom],  dtype=np.int32)

    wp_basis, resoln_b = build_wp_basis(species_list_ang, sp_names_hsd)
    print(f"  Basis: {len(wp_basis)} species, resoln={resoln_b:.4f} Bohr")

    # norb per atom (from basis)
    sp_by_name = {sp['name']: sp for sp in species_list_ang}
    norb_per_atom = np.array([sum(2*o['l']+1 for o in sp_by_name[species_names[si]]['orbitals'])
                               for si in species_per_atom], dtype=np.int32)
    assert norb_per_atom.sum() == norb_total, f"norb mismatch: basis sums to {norb_per_atom.sum()}, XML says {norb_total}"

    # ---- 3. Parse eigenvectors ----
    # parse_eigenvec_bin returns (nstates, norb) directly
    evecs_full = parse_eigenvec_bin_custom(dftb_dir / 'eigenvec.bin',  nstates_total, norb_total)

    # ---- HOMO ----
    homo_idx = np.where(occs_flat > 0.5)[0]
    homo = int(homo_idx[-1]) + 1 if len(homo_idx) else nstates_total // 2

    # Read energies from band.out (optional, for labels)
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

    # ---- MO range ----
    if args.mo_range:
        mo_start = max(1, args.mo_range[0])
        mo_end   = min(nstates_total, args.mo_range[1])
    else:
        mo_start = 1
        mo_end   = min(nstates_total, args.nmo)
    nstates  = mo_end - mo_start + 1
    evecs    = evecs_full[mo_start-1 : mo_end]
    energies = energies_ev[mo_start-1 : mo_end]
    occs     = occs_flat[mo_start-1 : mo_end]

    print(f"  HOMO = MO{homo} ({energies_ev[homo-1]:.3f} eV)")
    print(f"  Evaluating MOs {mo_start}–{mo_end}  ({nstates} total)")
    for i in range(min(nstates, 15)):
        tag = " <--HOMO" if (mo_start+i)==homo else (" <--LUMO" if (mo_start+i)==homo+1 else "")
        print(f"    MO{mo_start+i:4d}  E={energies[i]:8.3f} eV  occ={occs[i]:.1f}{tag}")

    # ================================================================
    # Method 2: point evaluation
    # ================================================================
    if args.points:
        print("\n" + "=" * 60)
        scan_type = args.plane2d.upper() if args.plane2d else (args.line_scan.upper() if args.line_scan else 'z-scan')
        print(f"Method 2: {scan_type}")
        print("=" * 60)

        atom_ang = atom_coords_b * BOHR2ANG
        if args.xy_range:
            rmin, rmax_r = args.xy_range
        else:
            rmin  = float(atom_ang.min()) - 3.0
            rmax_r = float(atom_ang.max()) + 3.0

        if args.plane2d:
            points_ang, extent = generate_2d_point_grid(args.plane2d, args.npoints, args.z_offset, (rmin, rmax_r))
            if args.plane2d == 'xy':
                plane_desc = f"XY plane  z={args.z_offset:.3f} Å"
            elif args.plane2d == 'xz':
                plane_desc = f"XZ plane  y={args.z_offset:.3f} Å"
            elif args.plane2d == 'yz':
                plane_desc = f"YZ plane  x={args.z_offset:.3f} Å"
        elif args.line_scan == 'bond':
            # Line scan through two atoms
            if args.atoms is None:
                # Default: use first two atoms
                i_atom, j_atom = 0, 1
            else:
                i_atom, j_atom = args.atoms
            print(f"  Bond line scan: atom {i_atom} -> atom {j_atom}")
            pos_i = atom_ang[i_atom]
            pos_j = atom_ang[j_atom]
            bond_vec = pos_j - pos_i
            bond_len = np.linalg.norm(bond_vec)
            bond_dir = bond_vec / bond_len
            print(f"  Bond length: {bond_len:.3f} Å")
            print(f"  Bond direction: ({bond_dir[0]:.3f}, {bond_dir[1]:.3f}, {bond_dir[2]:.3f})")
            
            if args.line_range:
                rmin, rmax_r = args.line_range
            else:
                rmin = -3.0
                rmax_r = bond_len + 3.0
            
            t_vals = np.linspace(rmin, rmax_r, args.npoints)  # distance along bond from atom i
            points_ang = np.array([pos_i + t * bond_dir for t in t_vals])
            ax_labels = ('distance along bond (Å)', 'ψ')
            plane_desc = f"Bond line scan  atom{i_atom}→atom{j_atom}  bond_len={bond_len:.3f}Å"
        else:  # 1D z-scan (default)
            z_vals = np.linspace(args.z_range[0], args.z_range[1], args.npoints)
            points_ang = np.column_stack([np.zeros(args.npoints), np.zeros(args.npoints), z_vals])
            ax_labels = ('z (Å)', 'ψ')
            plane_desc = "1D z-scan  x=0  y=0"

        npts = len(points_ang)
        print(f"  {npts} points, coord range [{rmin:.2f}, {rmax_r:.2f}] Å")
        print(f"  {plane_desc}")

        # --- libwaveplot ---
        points_bohr = points_ang / BOHR2ANG
        dftb_data_wp = { 'coords_bohr': atom_coords_b, 'species_wp': species_wp, 'basis': wp_basis, 'resolution': resoln_b, 'evecs': evecs }
        wp = setup_waveplot_from_dftb(dftb_data_wp, LIB_PATH)
        mo_indices_wp = list(range(1, nstates + 1))  # 1-indexed for libwaveplot
        wp_vals_list = wp_evaluate_mos_on_points(wp, mo_indices_wp, points_bohr)
        wp_vals = np.array(wp_vals_list)
        print(f"  libwaveplot done. max|ψ|={np.abs(wp_vals).max():.4e}")

        # --- OpenCL ---
        dftb_data_ocl = { 'coords_bohr': atom_coords_b, 'species_per_atom': species_per_atom,  'species_names': species_names}
        projector, atoms_dict = setup_gridprojector_from_dftb(dftb_data_ocl, species_list_ang, verbosity=0)
        mo_indices_ocl = list(range(nstates))  # 0-indexed for OpenCL
        ocl_vals_list = ocl_evaluate_mos_on_points(  projector, mo_indices_ocl, points_ang.astype(np.float32), evecs, natoms, species_per_atom, species_names, species_list_ang, norb_per_atom, atoms_dict )
        ocl_vals = np.array(ocl_vals_list)
        print(f"  OpenCL done. max|ψ|={np.abs(ocl_vals).max():.4e}")

        # --- comparison ---
        diff_vals = wp_vals - ocl_vals
        mo_indices = list(range(mo_start, mo_end + 1))
        results = compare_point_evaluations(wp_vals, ocl_vals, mo_indices, energies, homo)
        print_comparison_results(results, "libwaveplot vs OpenCL")

        # --- plot ---
        method_tag = "orb2points"
        if args.plane2d:
            s2 = (args.npoints, args.npoints)
            wp_vals_2d = np.array([v.reshape(s2) for v in wp_vals])
            ocl_vals_2d = np.array([v.reshape(s2) for v in ocl_vals])
            diff_vals_2d = np.array([v.reshape(s2) for v in diff_vals])
            suffix = f"{system_name}_{method_tag}_{args.plane2d}_z{args.z_offset:.2f}A_n{args.npoints}_MO{mo_start}-{mo_end}"
            out_file = OUTPUT_DIR / f'comparison_points_{suffix}.png'
            # Convert atom coordinates to Angstrom for overlay
            atom_coords_ang = atom_coords_b * BOHR2ANG
            plot_comparison_2d(wp_vals_2d, ocl_vals_2d, diff_vals_2d, extent, system_name, plane_desc, method_tag, mo_indices, energies, homo, out_file, args.dpi, atom_coords_ang)
        else:
            wp_vals_1d = wp_vals
            ocl_vals_1d = ocl_vals
            if args.line_scan == 'bond':
                x_vals = t_vals
            else:
                x_vals = z_vals
            suffix = f"{system_name}_{method_tag}_{args.line_scan if args.line_scan else 'zscan'}_n{args.npoints}_MO{mo_start}-{mo_end}"
            out_file = OUTPUT_DIR / f'comparison_points_{suffix}.png'
            plot_comparison_1d(x_vals, wp_vals_1d, ocl_vals_1d, system_name, plane_desc, method_tag, mo_indices, energies, homo, out_file, args.dpi)
        print(f"\nSaved: {out_file}")
        if not args.no_show: plt.show()
        return

    # ================================================================
    # Method 1: 3D grid comparison vs cube files
    # ================================================================
    print("\n" + "=" * 60)
    print("Method 1: 3D grid (libwaveplot vs WAVEPLOT cube)")
    print("=" * 60)

    cube_grids, origin_b, step_b, nPoints = [], None, None, None
    for imo in range(mo_start, mo_end+1):
        cube_path = dftb_dir / f'wp-1-1-{imo}-real.cube'
        if not cube_path.exists():
            print(f"  WARNING: {cube_path} not found — skipping Method 1")
            return
        data, _ori, _stp, _ac, _az = read_cube(cube_path)
        cube_grids.append(data)
        if origin_b is None:
            origin_b = _ori; step_b = _stp
            nPoints  = np.array(data.shape, dtype=np.int32)
        print(f"  MO{imo}: shape={data.shape}  |ψ|max={np.abs(data).max():.4e}")

    wp = WavePlot(LIB_PATH)
    wp.set_geometry(atom_coords_b, species_wp, is_periodic=False)
    wp.set_basis(wp_basis, resolution=resoln_b)
    wp.set_eigenvectors(evecs)

    gridVecs = np.diag(step_b)
    lib_grids = []
    for i in range(nstates):
        g = wp.orb2grid(i+1, origin_b, gridVecs, nPoints)
        lib_grids.append(g)
        print(f"  libwaveplot MO{mo_start+i}: |ψ|max={np.abs(g).max():.4e}")

    iz_mol   = int(np.clip(round(-origin_b[2]/step_b[2]), 0, nPoints[2]-1))
    z_slice_ang = (origin_b[2] + iz_mol * step_b[2]) * BOHR2ANG
    ext_xy = [origin_b[0]*BOHR2ANG, (origin_b[0]+nPoints[0]*step_b[0])*BOHR2ANG, origin_b[1]*BOHR2ANG, (origin_b[1]+nPoints[1]*step_b[1])*BOHR2ANG]
    method1_tag = "orb2grid"
    slice_desc  = f"XY slice  iz={iz_mol}  z={z_slice_ang:.3f} Å  [{method1_tag}]"
    print(f"  {slice_desc}")

    fig, axes = plt.subplots(nstates, 3, figsize=(13, 4*nstates))
    if nstates == 1: axes = axes[np.newaxis, :]
    all_stats = []
    for i in range(nstates):
        ref  = cube_grids[i]; lib  = lib_grids[i]; diff = lib - ref
        rs   = ref[:,:,iz_mol]; ls  = lib[:,:,iz_mol]; ds = diff[:,:,iz_mol]
        clim = max(np.abs(rs).max(), np.abs(ls).max()) or 1.0
        tag  = " [HOMO]" if (mo_start+i)==homo else (" [LUMO]" if (mo_start+i)==homo+1 else "")
        for ax, dat, ttl, cm, vl, vh in [
            (axes[i,0], rs, f"WAVEPLOT cube MO{mo_start+i}{tag}\nE={energies[i]:.2f}eV  {slice_desc}", 'RdBu_r', -clim, clim),
            (axes[i,1], ls, f"libwaveplot.so MO{mo_start+i}{tag}\n{slice_desc}",                        'RdBu_r', -clim, clim),
            (axes[i,2], ds, f"diff (lib−cube)\nRMS={np.sqrt(np.mean(diff**2)):.2e}",                    'bwr',    -clim*0.05, clim*0.05),
        ]:
            im = ax.imshow(dat.T, origin='lower', cmap=cm, vmin=vl, vmax=vh,  extent=ext_xy, interpolation='bilinear')
            ax.set_xlabel('x (Å)'); ax.set_ylabel('y (Å)')
            ax.set_title(ttl, fontsize=7)
            plt.colorbar(im, ax=ax, fraction=0.046)
        rms = np.sqrt(np.mean(diff**2)); rmax = np.abs(ref).max()
        all_stats.append((mo_start+i, rms, rms/rmax if rmax>1e-10 else float('nan'), rmax))

    fig.suptitle(f"{system_name}: libwaveplot vs WAVEPLOT cube  [{method1_tag}]  {slice_desc}  (MO{mo_start}–{mo_end})", fontsize=9)
    fig.tight_layout()
    out1 = OUTPUT_DIR / f'comparison_lib_vs_cube_{system_name}_MO{mo_start}-{mo_end}.png'
    fig.savefig(str(out1), dpi=args.dpi)
    print(f"\nSaved: {out1}")
    if not args.no_show: plt.show()
    plt.close(fig)

    print(f"\n{'MO':>5}  {'E(eV)':>8}  {'RMS':>12}  {'RelRMS':>10}  {'|ref|max':>12}")
    for mo, rms, rel, rmax in all_stats:
        i = mo - mo_start
        print(f"  {mo:3d}   {energies[i]:8.3f}   {rms:12.4e}   {rel:10.4e}   {rmax:12.4e}")
    print(f"\nOutputs in: {OUTPUT_DIR}")


if __name__ == '__main__':
    main()
