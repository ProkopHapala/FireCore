#!/usr/bin/env python3
"""
run_dftb_zscan.py - Run DFTB rigid z-scan for CO tip approaching molecule.

Produces zscan_z.npy and zscan_energy_eV.npy in the output directory.
Can use any DFTB+ basis set (mio-1-1, 3ob-3-1, etc.).

Usage:
    cd tests/tAFM/pyocl_fdbm
    PYTHONPATH=/path/to/FireCore:$PYTHONPATH python run_dftb_zscan.py --basis mio-1-1
    PYTHONPATH=/path/to/FireCore:$PYTHONPATH python run_dftb_zscan.py --basis 3ob-3-1 --xyz mymol.xyz
"""

import os, sys, argparse, time
import numpy as np

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.realpath(os.path.join(_THIS_DIR, '..', '..', '..'))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from pyBall import atomicUtils as au
from pyBall import dftb_utils as du

HAU2EV = 27.211386245988

# ── Main z-scan ───────────────────────────────────────────────────────────
def run_zscan_for_atom(target_idx, mol_names, mol_pos, tip_names, sk_prefix, z_distances, out_dir):
    """Run DFTB z-scan for a single target atom. Returns z_vals, e_vals."""
    target_name = mol_names[target_idx]
    target_pos = mol_pos[target_idx]
    atom_dir = os.path.join(out_dir, f'atom_{target_idx}')
    os.makedirs(atom_dir, exist_ok=True)

    print(f"\n{'='*60}")
    print(f"Target atom {target_idx}: {target_name} at [{target_pos[0]:.4f}, {target_pos[1]:.4f}, {target_pos[2]:.4f}]")
    print(f"Output: {atom_dir}")
    print(f"{'='*60}")

    cache_path = os.path.join(atom_dir, 'zscan_results_cache.npz')
    results = []
    if os.path.exists(cache_path):
        print(f"Loading cache: {cache_path}")
        cache = np.load(cache_path, allow_pickle=True)
        cached_z = cache['z_distances']
        if np.allclose(cached_z, z_distances):
            results = cache['results'].tolist()
            print(f"  Using {len(results)} cached results")
        else:
            print("  Cache z-range mismatch, recomputing")

    if len(results) != len(z_distances):
        combined_names = list(mol_names) + list(tip_names)
        for iz, z in enumerate(z_distances):
            print(f"\n[z-scan {iz+1}/{len(z_distances)}] z = {z:.2f} Å")
            o_pos = np.array([target_pos[0], target_pos[1], target_pos[2] + z])
            c_pos = np.array([target_pos[0], target_pos[1], target_pos[2] + z + 1.13])
            co_pos_shifted = np.array([o_pos, c_pos])
            combined_pos = np.vstack([mol_pos, co_pos_shifted])

            work_dir = os.path.join(atom_dir, f'zscan_z{z:.2f}')
            t_start = time.time()
            try:
                energy_ha = du.run_dftb_sp(work_dir, combined_names, combined_pos, sk_prefix)
                energy_ev = energy_ha * HAU2EV
                t_elapsed = time.time() - t_start
                print(f"  Energy: {energy_ha:.8f} Ha = {energy_ev:.6f} eV  ({t_elapsed:.1f}s)")
                results.append({'z': float(z), 'energy_Ha': float(energy_ha), 'energy_eV': float(energy_ev)})
            except Exception as e:
                print(f"  ERROR: {e}")
                raise

        np.savez(cache_path, results=results, z_distances=z_distances)
        print(f"\nSaved cache to {cache_path}")

    z_vals = np.array([r['z'] for r in results])
    e_vals = np.array([r['energy_eV'] for r in results])
    np.save(os.path.join(atom_dir, 'zscan_z.npy'), z_vals)
    np.save(os.path.join(atom_dir, 'zscan_energy_eV.npy'), e_vals)

    with open(os.path.join(atom_dir, 'zscan_results.txt'), 'w') as f:
        f.write("DFTB Z-Scan Results\n")
        f.write("="*70 + "\n")
        f.write(f"Target: {target_name} [{target_pos[0]:.4f}, {target_pos[1]:.4f}, {target_pos[2]:.4f}]\n")
        f.write(f"CO bond: 1.13 Å\n\n")
        f.write(f"{'z [Å]':>10}  {'E [Ha]':>16}  {'E [eV]':>16}\n")
        f.write("-"*70 + "\n")
        for r in results:
            f.write(f"{r['z']:10.2f}  {r['energy_Ha']:16.8f}  {r['energy_eV']:16.6f}\n")

    e_rel = e_vals - e_vals[-1]
    print(f"\nAtom {target_idx} complete: {len(z_vals)} points")
    print(f"  Rel energy at contact: {e_rel[0]:.4f} eV")
    return z_vals, e_vals

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--basis', type=str, default='mio-1-1', choices=['mio-1-1', '3ob-3-1'])
    parser.add_argument('--xyz', type=str, default='pentacene.xyz', help='Molecule XYZ file')
    parser.add_argument('--tip_xyz', type=str, default='CO.xyz', help='Tip molecule XYZ file')
    parser.add_argument('--output_dir', type=str, default=None)
    parser.add_argument('--target_indices', type=str, default='0',
                        help='Comma-separated list of target atom indices (e.g., "0,1,20,21")')
    parser.add_argument('--z_min', type=float, default=2.0, help='Minimum z-distance [Å]')
    parser.add_argument('--z_max', type=float, default=10.0, help='Maximum z-distance [Å]')
    parser.add_argument('--z_step', type=float, default=0.15, help='z step [Å]')
    args = parser.parse_args()

    basis = args.basis
    sk_prefix = du.SK_PATHS[basis]
    out_dir = os.path.join(_THIS_DIR, args.output_dir) if args.output_dir else os.path.join(_THIS_DIR, f'zscan_{basis.replace("-", "_")}')
    os.makedirs(out_dir, exist_ok=True)

    target_indices = [int(x.strip()) for x in args.target_indices.split(',')]

    print(f"="*60)
    print(f"DFTB Z-Scan: basis={basis}, atoms={target_indices}")
    print(f"Output: {out_dir}")
    print(f"="*60)

    xyz_path = os.path.join(_THIS_DIR, args.xyz)
    tip_path = os.path.join(_THIS_DIR, args.tip_xyz)
    mol_pos, _, mol_names, _, _ = au.load_xyz(xyz_path)
    mol_pos = np.array(mol_pos, dtype=np.float64)
    tip_pos, _, tip_names, _, _ = au.load_xyz(tip_path)
    print(f"Molecule: {len(mol_names)} atoms ({args.xyz})")
    print(f"Tip:      {len(tip_names)} atoms ({args.tip_xyz})")

    z_distances = np.arange(args.z_min, args.z_max + args.z_step*0.5, args.z_step)
    print(f"Z-scan: {len(z_distances)} points from {z_distances.min():.2f} to {z_distances.max():.2f} Å")

    for target_idx in target_indices:
        run_zscan_for_atom(target_idx, mol_names, mol_pos, tip_names, sk_prefix, z_distances, out_dir)

    print(f"\nAll atoms done. Outputs in: {out_dir}/")

if __name__ == "__main__":
    main()
