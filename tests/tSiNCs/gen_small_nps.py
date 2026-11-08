#!/usr/bin/env python3
"""
gen_small_nps.py - Generate small regular Si/C nanoparticles for PySCF vibration calculations.

Uses build_spherical_nanoparticle from pyBall.nanocrystal_gen to create H-passivated
spherical nanoparticles from diamond/sphalerite primitive cells.

Targets: ~30-50 atoms (a bit bigger than adamantane C10H16 / Si10H16 at 26 atoms).

Usage:
    python3 gen_small_nps.py                    # generate all default sizes
    python3 gen_small_nps.py --element Si       # only Si
    python3 gen_small_nps.py --element C        # only C (diamond)
    python3 gen_small_nps.py --outdir geometry  # output directory
"""
import argparse
import os
import sys
import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, REPO)

from pyBall.nanocrystal_gen import build_spherical_nanoparticle, save_xyz


# Target: a bit bigger than adamantane (26 atoms) / Si10H16 (26 atoms)
# (R, nrep) pairs to try — will pick those giving 28-60 atoms
SI_RADII = [4.2, 5.0, 6.0]
C_RADII  = [2.7, 3.5, 4.0]

# Also include reference structures from existing files
REFERENCES = {
    'adamantane':  os.path.join(REPO, 'cpp/common_resources/mol/adamantane.mol2'),
    'Si10H16':     os.path.join(REPO, 'cpp/common_resources/xyz/Si10_H.xyz'),
}
# Fallback: adamantane as XYZ from geometry/g4.xyz
ADAMANTANE_XYZ_FALLBACK = os.path.join(os.path.dirname(__file__), 'geometry', 'g4.xyz')


def load_reference_xyz_or_mol2(path):
    """Load reference structure from .xyz or .mol2 into ASE Atoms."""
    from ase.io import read
    if path.endswith('.mol2'):
        try:
            return read(path, format='mol2')
        except Exception:
            # ASE mol2 reader may fail on non-standard files; try generic read
            return read(path)
    return read(path)


def try_radii(element, radii, nrep, outdir, min_atoms=28, max_atoms=60):
    """Generate nanoparticles at several radii, save those in target size range."""
    heavy_z = 14 if element.lower() == 'si' else 6
    results = []
    for R in radii:
        try:
            elems, apos = build_spherical_nanoparticle(
                R=R, nrep=nrep, heavy_z=heavy_z,
                resolve_clashes=True,
            )
            nat = len(elems)
            n_heavy = int(np.sum(elems != 'H'))
            n_h = int(np.sum(elems == 'H'))
            tag = f'{element}_sphere_R{R:.1f}_nat{nat}'
            xyz_path = os.path.join(outdir, f'{tag}.xyz')
            save_xyz(xyz_path, elems, apos, comment=f'{element} sphere R={R} nrep={nrep} natoms={nat} nHeavy={n_heavy} nH={n_h}')
            in_range = min_atoms <= nat <= max_atoms
            print(f"  R={R:.1f}  natoms={nat:3d}  (heavy={n_heavy}, H={n_h})  {'<-- in range' if in_range else ''}")
            results.append({'tag': tag, 'path': xyz_path, 'natoms': nat, 'n_heavy': n_heavy, 'n_h': n_h, 'R': R, 'in_range': in_range})
        except Exception as e:
            print(f"  R={R:.1f}  FAILED: {e}")
    return results


def main():
    ap = argparse.ArgumentParser(description='Generate small regular nanoparticles for PySCF vibrations')
    ap.add_argument('--element', choices=['Si', 'C', 'both'], default='both')
    ap.add_argument('--outdir', default=os.path.join(os.path.dirname(__file__), 'geometry', 'small_nps'))
    ap.add_argument('--min-atoms', type=int, default=28, help='Min total atoms (default 28, a bit > 26)')
    ap.add_argument('--max-atoms', type=int, default=60, help='Max total atoms for PySCF feasibility')
    ap.add_argument('--nrep', type=int, default=5, help='Lattice repetition for sphere cut')
    ap.add_argument('--include-references', action='store_true', default=True, help='Copy adamantane/Si10H16 as references')
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    all_results = {}

    # Copy reference structures
    if args.include_references:
        print("\n=== Reference structures ===")
        for name, path in REFERENCES.items():
            if not os.path.exists(path):
                print(f"  {name}: file not found at {path}")
                continue
            try:
                atoms = load_reference_xyz_or_mol2(path)
                nat = len(atoms)
                ref_xyz = os.path.join(args.outdir, f'{name}.xyz')
                from ase.io import write
                write(ref_xyz, atoms, format='xyz')
                print(f"  {name}: {nat} atoms -> {ref_xyz}")
                all_results[name] = {'tag': name, 'path': ref_xyz, 'natoms': nat}
            except Exception as e:
                # Try XYZ fallback for adamantane
                if name == 'adamantane' and os.path.exists(ADAMANTANE_XYZ_FALLBACK):
                    atoms = load_reference_xyz_or_mol2(ADAMANTANE_XYZ_FALLBACK)
                    nat = len(atoms)
                    ref_xyz = os.path.join(args.outdir, f'{name}.xyz')
                    from ase.io import write
                    write(ref_xyz, atoms, format='xyz')
                    print(f"  {name}: {nat} atoms -> {ref_xyz} (from XYZ fallback)")
                    all_results[name] = {'tag': name, 'path': ref_xyz, 'natoms': nat}
                else:
                    print(f"  {name}: FAILED to load: {e}")

    # Generate Si nanoparticles
    if args.element in ('Si', 'both'):
        print(f"\n=== Si nanoparticles (nrep={args.nrep}) ===")
        si_results = try_radii('Si', SI_RADII, args.nrep, args.outdir, args.min_atoms, args.max_atoms)
        all_results['Si'] = si_results

    # Generate C nanoparticles
    if args.element in ('C', 'both'):
        print(f"\n=== C nanoparticles (nrep={args.nrep}) ===")
        c_results = try_radii('C', C_RADII, args.nrep, args.outdir, args.min_atoms, args.max_atoms)
        all_results['C'] = c_results

    # Summary
    print(f"\n=== Summary ===")
    print(f"Output directory: {args.outdir}")
    for key, val in all_results.items():
        if isinstance(val, list):
            in_range = [r for r in val if r['in_range']]
            print(f"  {key}: {len(in_range)} structures in range [{args.min_atoms}, {args.max_atoms}] atoms")
            for r in in_range:
                print(f"    {r['tag']}  ({r['natoms']} atoms, {r['n_heavy']} heavy + {r['n_h']} H)")
        else:
            print(f"  {key}: {val['natoms']} atoms -> {val['path']}")

    # Write index JSON
    import json
    index_path = os.path.join(args.outdir, 'index.json')
    index = {}
    for key, val in all_results.items():
        if isinstance(val, list):
            index[key] = [{'tag': r['tag'], 'path': r['path'], 'natoms': r['natoms'],
                           'n_heavy': r['n_heavy'], 'n_h': r['n_h'], 'R': r['R']}
                          for r in val if r['in_range']]
        else:
            index[key] = [{'tag': val['tag'], 'path': val['path'], 'natoms': val['natoms']}]
    with open(index_path, 'w') as f:
        json.dump(index, f, indent=2)
    print(f"\nIndex: {index_path}")


if __name__ == '__main__':
    main()
