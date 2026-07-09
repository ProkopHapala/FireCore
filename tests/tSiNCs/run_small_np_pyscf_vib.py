#!/usr/bin/env python3
"""
run_small_np_pyscf_vib.py - Full PySCF relax + vibration pipeline for small nanoparticles.

Generates small regular Si/C nanoparticles (or loads existing XYZ),
relaxes them with PySCF, computes Hessian + harmonic analysis,
and saves hessian, eigenmodes, frequencies, geometry to .npz.

Usage:
    python3 run_small_np_pyscf_vib.py                          # run all from geometry/small_nps/
    python3 run_small_np_pyscf_vib.py --names Si_sphere_R5.0_nat42  # specific structure
    python3 run_small_np_pyscf_vib.py --method hf --basis sto-3g    # PySCF settings
    python3 run_small_np_pyscf_vib.py --gen-only                    # just generate, don't run PySCF
    python3 run_small_np_pyscf_vib.py --export-modes                # also export modes to XYZ
"""
import argparse
import os
import sys
import json
import numpy as np
from pathlib import Path

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.dirname(__file__))

from vib_utils import run_pyscf_vib_full, export_modes_to_xyz, vib_npz_path, relaxed_xyz_path


def load_xyz_to_ase(xyz_path):
    """Load XYZ file into ASE Atoms."""
    from ase.io import read
    return read(xyz_path)


def discover_structures(geom_dir):
    """Find all .xyz files in geometry directory (excluding *_relaxed.xyz, *_all_modes.xyz)."""
    structures = []
    for f in sorted(os.listdir(geom_dir)):
        if not f.endswith('.xyz'):
            continue
        if '_relaxed' in f or '_all_modes' in f:
            continue
        xyz_path = os.path.join(geom_dir, f)
        tag = f.replace('.xyz', '')
        structures.append({'tag': tag, 'path': xyz_path})
    return structures


def run_pipeline(structures, method, basis, workdir, export_modes=False, threshold=10.0):
    """Run PySCF relax + vib pipeline for each structure."""
    results = []
    for s in structures:
        tag = s['tag']
        xyz_path = s['path']
        print(f"\n{'='*60}")
        print(f"  {tag}")
        print(f"{'='*60}")

        try:
            atoms = load_xyz_to_ase(xyz_path)
            nat = len(atoms)
            print(f"  Loaded {nat} atoms from {xyz_path}")
        except Exception as e:
            print(f"  FAILED to load: {e}")
            results.append({'tag': tag, 'status': 'load_failed', 'error': str(e)})
            continue

        try:
            data, method_tag = run_pyscf_vib_full(
                tag, atoms, method=method, basis=basis, workdir=workdir,
            )
            npz = vib_npz_path(tag, method_tag, workdir)
            freqs = data['freqs']
            n_real = np.sum((freqs.real > threshold) if np.iscomplexobj(freqs) else (freqs > threshold))
            print(f"  OK: {len(freqs)} modes ({n_real} real > {threshold} cm^-1)")
            print(f"  NPZ: {npz}")

            if export_modes:
                export_modes_to_xyz(atoms, tag, method_tag, workdir=workdir, threshold=threshold)
                modes_xyz = Path(workdir) / f'{tag}_{method_tag}_all_modes.xyz'
                print(f"  Modes XYZ: {modes_xyz}")

            results.append({
                'tag': tag, 'status': 'ok', 'natoms': nat,
                'npz': str(npz), 'nmodes': len(freqs), 'n_real': int(n_real),
                'energy': float(data['energy']),
            })
        except Exception as e:
            import traceback
            traceback.print_exc()
            print(f"  FAILED: {e}")
            results.append({'tag': tag, 'status': 'failed', 'error': str(e)})

    return results


def main():
    ap = argparse.ArgumentParser(description='PySCF relax + vibration pipeline for small nanoparticles')
    ap.add_argument('--names', nargs='+', default=None, help='Specific structure tags to run (default: all in geometry dir)')
    ap.add_argument('--method', default='hf', help='PySCF method (hf, b3lyp, pbe, ...)')
    ap.add_argument('--basis', default='sto-3g', help='Basis set (sto-3g, 6-31g(d), ...)')
    ap.add_argument('--workdir', default=os.path.dirname(__file__), help='Working directory for cached files')
    ap.add_argument('--geom-dir', default=os.path.join(os.path.dirname(__file__), 'geometry', 'small_nps'),
                    help='Directory with input XYZ structures')
    ap.add_argument('--gen-only', action='store_true', help='Only generate nanoparticles, skip PySCF')
    ap.add_argument('--export-modes', action='store_true', help='Export vibration modes to multi-frame XYZ')
    ap.add_argument('--threshold', type=float, default=10.0, help='Frequency threshold for real modes (cm^-1)')
    args = ap.parse_args()

    workdir = args.workdir
    geom_dir = args.geom_dir

    # Step 1: Generate nanoparticles if needed
    if not os.path.exists(geom_dir) or not any(f.endswith('.xyz') for f in os.listdir(geom_dir)):
        print(f"[gen] No structures found in {geom_dir}, generating...")
        gen_script = os.path.join(os.path.dirname(__file__), 'gen_small_nps.py')
        import subprocess
        subprocess.run([sys.executable, gen_script, '--outdir', geom_dir], check=True)

    # Step 2: Discover structures
    all_structures = discover_structures(geom_dir)
    if not all_structures:
        print(f"No .xyz structures found in {geom_dir}")
        return

    # Filter by names if specified
    if args.names:
        name_set = set(args.names)
        structures = [s for s in all_structures if s['tag'] in name_set]
        missing = name_set - {s['tag'] for s in structures}
        if missing:
            print(f"WARNING: structures not found: {missing}")
            print(f"Available: {[s['tag'] for s in all_structures]}")
    else:
        structures = all_structures

    print(f"\n[run] {len(structures)} structures to process")
    print(f"  method={args.method}, basis={args.basis}")
    print(f"  workdir={workdir}")
    for s in structures:
        print(f"  - {s['tag']}  ({s['path']})")

    if args.gen_only:
        print("\n[gen-only] Skipping PySCF calculations")
        return

    # Step 3: Run PySCF pipeline
    results = run_pipeline(structures, args.method, args.basis, workdir,
                           export_modes=args.export_modes, threshold=args.threshold)

    # Step 4: Summary
    print(f"\n{'='*60}")
    print(f"  SUMMARY")
    print(f"{'='*60}")
    ok = [r for r in results if r['status'] == 'ok']
    failed = [r for r in results if r['status'] != 'ok']
    print(f"  OK: {len(ok)} / {len(results)}")
    for r in ok:
        print(f"    {r['tag']:40s}  natoms={r['natoms']:3d}  nmodes={r['nmodes']}  E={r['energy']:.8f} Ha")
        print(f"      -> {r['npz']}")
    if failed:
        print(f"  FAILED: {len(failed)}")
        for r in failed:
            print(f"    {r['tag']:40s}  {r['status']}: {r.get('error', '?')}")

    # Write results JSON
    results_path = os.path.join(workdir, 'small_np_pyscf_vib_results.json')
    with open(results_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults: {results_path}")


if __name__ == '__main__':
    main()
