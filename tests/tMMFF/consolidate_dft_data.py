#!/usr/bin/env python3
"""
CLI script to consolidate DFT reference data for FDBM fitting.

Organizes .xyz geometries and .dat energies, generates NPZ files and XYZ movies.

Usage:
    python consolidate_dft_data.py consolidate --molecule H2O-H
    python consolidate_dft_data.py match --molecule H2O-H
    python consolidate_dft_data.py movie --molecule H2O-H --movie-type scan_z --site Na --tilt 0
"""
import sys
import os
import argparse
import numpy as np

sys.path.insert(0, '/home/prokop/git/FireCore')
from pyBall import atomicUtils as au

# Default data directory - all paths derived from this
DEFAULT_DATA_DIR = '../tSurf/small_mols_NaCl_New'


def cmd_consolidate(args):
    """Consolidate all XYZ files into single NPZ file."""
    data_dir = args.data_dir or DEFAULT_DATA_DIR
    base_dir = args.base_dir or os.path.join(data_dir, '1-inputs')
    output_file = os.path.join(base_dir, f'{args.molecule}_consolidated.npz')
    result = au.consolidate_xyz_to_npz(
        base_dir=base_dir,
        molecule=args.molecule,
        n_substrate=args.n_substrate,
        output_file=output_file
    )
    print(f"Consolidated {args.molecule} to: {result}")
    # Print summary
    data = np.load(result)
    print(f"  Coords shape: {data['coords'].shape}")
    print(f"  Element names: {list(data['enames'])}")
    print(f"  Unique ix: {np.unique(data['ix'])}")
    print(f"  Unique iy: {np.unique(data['iy'])}")
    print(f"  Sites: {np.unique(data['orientxy'])}")
    print(f"  Tilts: {np.unique(data['orientz'])}")
    print(f"  Z range: [{data['zdist'].min():.2f}, {data['zdist'].max():.2f}]")


def cmd_match(args):
    """Match XYZ files with DFT energies."""
    data_dir = args.data_dir or DEFAULT_DATA_DIR
    base_dir = args.base_dir or os.path.join(data_dir, '1-inputs')
    results_dir = args.results_dir or os.path.join(data_dir, '3-results')
    matched = au.match_xyz_with_energies(
        base_dir=base_dir,
        results_dir=results_dir,
        molecule=args.molecule,
        n_substrate=args.n_substrate
    )
    print(f"Matched {len(matched)} configurations for {args.molecule}")
    # Save matched data
    output_file = f'{args.molecule}_matched.npz'
    coords = np.array([m['apos'] for m in matched])
    energies = np.array([m['energy'] for m in matched])
    metadata = {
        'ix': np.array([m['ix'] for m in matched]),
        'iy': np.array([m['iy'] for m in matched]),
        'orientxy': np.array([m['orientxy'] for m in matched]),
        'orientz': np.array([m['orientz'] for m in matched]),
        'zdist': np.array([m['zdist'] for m in matched]),
    }
    np.savez(output_file,
             coords=coords,
             energies=energies,
             enames=matched[0]['enames'] if matched else [],
             **metadata)
    print(f"Saved matched data to: {output_file}")
    print(f"  Coords shape: {coords.shape}")
    print(f"  Energy range: [{energies.min():.6f}, {energies.max():.6f}] eV")


def cmd_movie(args):
    """Generate XYZ movie from matched data."""
    data_dir = args.data_dir or DEFAULT_DATA_DIR
    base_dir = args.base_dir or os.path.join(data_dir, '1-inputs')
    results_dir = args.results_dir or os.path.join(data_dir, '3-results')
    matched = au.match_xyz_with_energies(
        base_dir=base_dir,
        results_dir=results_dir,
        molecule=args.molecule,
        n_substrate=args.n_substrate
    )
    # Build filter function based on movie type
    if args.movie_type == 'scan_z':
        # Fix site and tilt, scan z distances
        def filter_func(d):
            return (d['orientxy'] == args.site and
                    d['orientz'] == args.tilt and
                    d['ix'] == args.ix and d['iy'] == args.iy)
        output_file = f'{args.molecule}_scan_z_ix{args.ix}_iy{args.iy}_{args.site}_tilt{args.tilt}.xyz'
    elif args.movie_type == 'scan_orientation':
        # Fix z and site, scan tilts
        def filter_func(d):
            return (d['zdist'] == args.zdist and
                    d['orientxy'] == args.site and
                    d['ix'] == args.ix and d['iy'] == args.iy)
        output_file = f'{args.molecule}_scan_orientation_ix{args.ix}_iy{args.iy}_{args.site}_z{args.zdist}.xyz'
    elif args.movie_type == 'scan_site':
        # Fix z and tilt, scan sites
        def filter_func(d):
            return (d['zdist'] == args.zdist and
                    d['orientz'] == args.tilt and
                    d['ix'] == args.ix and d['iy'] == args.iy)
        output_file = f'{args.molecule}_scan_site_ix{args.ix}_iy{args.iy}_z{args.zdist}_tilt{args.tilt}.xyz'
    elif args.movie_type == 'scan_xy':
        # Fix z, site, tilt, scan (ix, iy) positions
        def filter_func(d):
            return (d['zdist'] == args.zdist and
                    d['orientxy'] == args.site and
                    d['orientz'] == args.tilt)
        output_file = f'{args.molecule}_scan_xy_{args.site}_z{args.zdist}_tilt{args.tilt}.xyz'
    else:
        raise ValueError(f"Unknown movie type: {args.movie_type}")
    substrate_xyz = os.path.join(data_dir, '0-geoms/NaCl.xyz') if args.include_substrate else None
    au.generate_xyz_movie(
        matched_data=matched,
        output_file=output_file,
        include_substrate=args.include_substrate,
        substrate_xyz=substrate_xyz,
        n_substrate=args.n_substrate,
        filter_func=filter_func
    )
    print(f"Generated movie: {output_file}")


def cmd_extract(args):
    """Extract XYZ with z-min frame(s) for given (ix, iy, site, tilt) for alignment diagnostics."""
    data_dir = args.data_dir or DEFAULT_DATA_DIR
    base_dir = args.base_dir or os.path.join(data_dir, '1-inputs')
    results_dir = args.results_dir or os.path.join(data_dir, '3-results')
    matched = au.match_xyz_with_energies(
        base_dir=base_dir,
        results_dir=results_dir,
        molecule=args.molecule,
        n_substrate=args.n_substrate
    )
    sel = [m for m in matched if (m['ix'] == args.ix and m['iy'] == args.iy and m['orientxy'] == args.site and m['orientz'] == args.tilt)]
    if len(sel) == 0:
        raise RuntimeError(f"No matched frames for molecule={args.molecule} ix={args.ix} iy={args.iy} site={args.site} tilt={args.tilt}")
    zmax = max(m['zdist'] for m in sel)
    sel_zmax = [m for m in sel if m['zdist'] == zmax]
    if len(sel_zmax) == 0:
        raise RuntimeError('Internal error: empty zmax selection')
    E0 = float(np.mean([m['energy'] for m in sel_zmax]))
    zmin = min(m['zdist'] for m in sel)
    sel_zmin = [m for m in sel if m['zdist'] == zmin]
    if len(sel_zmin) == 0:
        raise RuntimeError('Internal error: empty zmin selection')
    sel_zmin = sorted(sel_zmin, key=lambda d: d['energy'])
    print(f"Selected {len(sel_zmin)} frame(s) at zmin={zmin:.6f} for {args.molecule} ix={args.ix} iy={args.iy} site={args.site} tilt={args.tilt}")
    E = np.array([m['energy'] for m in sel_zmin], dtype=float)
    print(f"  Energy range: [{E.min():.6f}, {E.max():.6f}] eV")
    print(f"  Baseline: zmax={zmax:.6f} E0=<E(zmax)>={E0:.6f} eV  (n={len(sel_zmax)})")
    Ei = E - E0
    print(f"  Interaction E_int=E-E0 range: [{Ei.min():.6f}, {Ei.max():.6f}] eV")

    if args.include_substrate:
        substrate_xyz = os.path.join(data_dir, '0-geoms/NaCl.xyz')
        def _load_xyz_simple(fname):
            with open(fname, 'r') as f:
                lines = f.readlines()
            n = int(lines[0].strip())
            names = []
            pos = np.zeros((n, 3), dtype=float)
            for i in range(n):
                w = lines[i + 2].split()
                names.append(w[0])
                pos[i, :] = (float(w[1]), float(w[2]), float(w[3]))
            return np.array(names), pos
        sub_names, sub_pos = _load_xyz_simple(substrate_xyz)
        z_top = np.max(sub_pos[:, 2])
        top_mask = np.isclose(sub_pos[:, 2], z_top, rtol=0.0, atol=1e-6)
        top_names = sub_names[top_mask]
        top_pos = sub_pos[top_mask]
        # use O as reference atom for lateral site
        ref = sel_zmin[0]['apos'][0]
        dxy2 = (top_pos[:, 0] - ref[0])**2 + (top_pos[:, 1] - ref[1])**2
        i0 = int(np.argmin(dxy2))
        print(f"  Top-layer z={z_top:.6f} A; nearest top atom under O: {top_names[i0]} at ({top_pos[i0,0]:.6f},{top_pos[i0,1]:.6f},{top_pos[i0,2]:.6f}) dxy={np.sqrt(dxy2[i0]):.6f} A")

    output_file = args.output
    if output_file is None:
        output_file = f"{args.molecule}_zmin_ix{args.ix}_iy{args.iy}_{args.site}_tilt{args.tilt}.xyz"
    substrate_xyz = os.path.join(data_dir, '0-geoms/NaCl.xyz') if args.include_substrate else None

    zmin_val = zmin
    def filter_func(d):
        return (d['ix'] == args.ix and d['iy'] == args.iy and d['orientxy'] == args.site and d['orientz'] == args.tilt and d['zdist'] == zmin_val)
    au.generate_xyz_movie(
        matched_data=matched,
        output_file=output_file,
        include_substrate=args.include_substrate,
        substrate_xyz=substrate_xyz,
        n_substrate=args.n_substrate,
        filter_func=filter_func
    )
    print(f"Wrote: {output_file}")


def cmd_summary(args):
    """Print summary of matched data."""
    data_dir = args.data_dir or DEFAULT_DATA_DIR
    base_dir = args.base_dir or os.path.join(data_dir, '1-inputs')
    results_dir = args.results_dir or os.path.join(data_dir, '3-results')
    matched = au.match_xyz_with_energies(
        base_dir=base_dir,
        results_dir=results_dir,
        molecule=args.molecule,
        n_substrate=args.n_substrate
    )
    print(f"\n=== Summary for {args.molecule} ===")
    print(f"Total matched configurations: {len(matched)}")
    # Group by parameters
    from collections import defaultdict
    by_site = defaultdict(int)
    by_tilt = defaultdict(int)
    by_ixiy = defaultdict(int)
    for m in matched:
        by_site[m['orientxy']] += 1
        by_tilt[m['orientz']] += 1
        by_ixiy[(m['ix'], m['iy'])] += 1
    print(f"\nBy surface site:")
    for site, count in sorted(by_site.items()):
        print(f"  {site}: {count}")
    print(f"\nBy tilt angle:")
    for tilt, count in sorted(by_tilt.items()):
        print(f"  {tilt}°: {count}")
    print(f"\nBy (ix, iy) position: {len(by_ixiy)} unique positions")
    # Energy stats
    energies = np.array([m['energy'] for m in matched])
    print(f"\nEnergy statistics:")
    print(f"  Min: {energies.min():.6f} eV")
    print(f"  Max: {energies.max():.6f} eV")
    print(f"  Mean: {energies.mean():.6f} eV")
    print(f"  Range: {energies.max() - energies.min():.6f} eV")


def main():
    parser = argparse.ArgumentParser(
        description='Consolidate DFT reference data for FDBM fitting',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Consolidate XYZ files to NPZ
  python consolidate_dft_data.py consolidate --molecule H2O-H

  # Match geometries with DFT energies
  python consolidate_dft_data.py match --molecule H2O-H

  # Generate movie scanning z distances
  python consolidate_dft_data.py movie --molecule H2O-H --movie-type scan_z --site Na --tilt 0 --ix 0 --iy 0

  # Generate movie with substrate atoms
  python consolidate_dft_data.py movie --molecule H2O-H --movie-type scan_z --site Na --tilt 0 --ix 0 --iy 0 --include-substrate

  # Print summary statistics
  python consolidate_dft_data.py summary --molecule H2O-H
""")
    parser.add_argument('--data-dir', default=DEFAULT_DATA_DIR, help='Data directory (default: ../tSurf/small_mols_NaCl_New)')
    parser.add_argument('--molecule', default='H2O-H', help='Molecule name (default: H2O-H)')
    parser.add_argument('--base-dir', default=None, help='Base directory with XYZ files (default: {data_dir}/1-inputs)')
    parser.add_argument('--results-dir', default=None, help='Directory with .dat files (default: {data_dir}/3-results)')
    parser.add_argument('--n-substrate', type=int, default=32, help='Number of substrate atoms (default: 32)')
    subparsers = parser.add_subparsers(dest='command', help='Command to run')
    # consolidate command
    p_consolidate = subparsers.add_parser('consolidate', help='Consolidate XYZ to NPZ')
    p_consolidate.set_defaults(func=cmd_consolidate)
    # match command
    p_match = subparsers.add_parser('match', help='Match XYZ with DFT energies')
    p_match.set_defaults(func=cmd_match)
    # movie command
    p_movie = subparsers.add_parser('movie', help='Generate XYZ movie')
    p_movie.add_argument('--movie-type', choices=['scan_z', 'scan_orientation', 'scan_site', 'scan_xy'], default='scan_z', help='Type of movie to generate')
    p_movie.add_argument('--site', choices=['Na', 'Cl', 'hollow'], default='Na',  help='Surface site for filtering')
    p_movie.add_argument('--tilt', type=int, choices=[-45, 0, 45], default=0,  help='Tilt angle for filtering')
    p_movie.add_argument('--ix', type=int, default=0, help='ix index')
    p_movie.add_argument('--iy', type=int, default=0, help='iy index')
    p_movie.add_argument('--zdist', type=float, default=2.0, help='z distance')
    p_movie.add_argument('--include-substrate', action='store_true',  help='Include substrate atoms in movie')
    p_movie.set_defaults(func=cmd_movie)
    # extract command
    p_extract = subparsers.add_parser('extract', help='Extract z-min XYZ for given (ix, iy, site, tilt)')
    p_extract.add_argument('--ix', type=int, default=0, help='ix index')
    p_extract.add_argument('--iy', type=int, default=0, help='iy index')
    p_extract.add_argument('--site', choices=['Na', 'Cl', 'hollow'], default='Na', help='Surface site for filtering')
    p_extract.add_argument('--tilt', type=int, choices=[-45, 0, 45], default=0, help='Tilt angle for filtering')
    p_extract.add_argument('--include-substrate', action='store_true', help='Include substrate atoms')
    p_extract.add_argument('--output', default=None, help='Output XYZ filename (default: auto)')
    p_extract.set_defaults(func=cmd_extract)
    # summary command
    p_summary = subparsers.add_parser('summary', help='Print summary statistics')
    p_summary.set_defaults(func=cmd_summary)
    args = parser.parse_args()
    if args.command is None:
        parser.print_help()
        return
    args.func(args)


if __name__ == '__main__':
    main()
