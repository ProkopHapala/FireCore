#!/usr/bin/env python3
"""
run_dftb_zscan.py - Run DFTB rigid z-scan for CO tip approaching pentacene.

Produces zscan_z.npy and zscan_energy_eV.npy in the output directory.
Can use any DFTB+ basis set (mio-1-1, 3ob-3-1, etc.).

Usage:
    cd tests/tAFM/pyocl_fdbm
    PYTHONPATH=/path/to/FireCore:$PYTHONPATH python run_dftb_zscan.py --basis mio-1-1
    PYTHONPATH=/path/to/FireCore:$PYTHONPATH python run_dftb_zscan.py --basis 3ob-3-1
"""

import os, sys, argparse, time
import numpy as np

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_PENTACENE_XYZ = os.path.join(_THIS_DIR, 'pentacene.xyz')
_CO_XYZ        = os.path.join(_THIS_DIR, 'CO.xyz')

BOHR2ANG = 0.5291772109
HAU2EV   = 27.211386245988

SLAKO_PREFIXES = {
    'mio-1-1': '/home/prokop/SIMULATIONS/dftbplus/slakos/mio-1-1/',
    '3ob-3-1': '/home/prokop/SIMULATIONS/dftbplus/slakos/3ob-3-1/',
}

# ── XYZ helpers ────────────────────────────────────────────────────────────
def load_xyz(fname):
    with open(fname, 'r') as f:
        lines = f.readlines()
    natoms = int(lines[0].strip())
    enames = []; apos = []
    for line in lines[2:2+natoms]:
        parts = line.split()
        enames.append(parts[0])
        apos.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return np.array(enames), np.array(apos, dtype=np.float64)

def write_xyz(fname, enames, apos, comment="Generated"):
    with open(fname, 'w') as f:
        f.write(f"{len(enames)}\n{comment}\n")
        for e, p in zip(enames, apos):
            f.write(f"{e} {p[0]:.10f} {p[1]:.10f} {p[2]:.10f}\n")

# ── DFTB+ helpers ─────────────────────────────────────────────────────────
def write_dftb_input(enames, xyz_path, out_path, sk_prefix):
    species = sorted(set(enames))
    max_ang = {}
    for s in species:
        max_ang[s] = '"s"' if s == 'H' else '"p"'
    max_ang_str = '\n    '.join([f'{s} = {max_ang[s]}' for s in species])

    hsd = f"""Geometry = xyzFormat {{
  <<< "{os.path.basename(xyz_path)}"
}}

Hamiltonian = DFTB {{
  SCC = Yes
  SlaterKosterFiles = Type2FileNames {{
    Prefix = "{sk_prefix}"
    Separator = "-"
    Suffix = ".skf"
  }}
  MaxAngularMomentum {{
    {max_ang_str}
  }}
  SCCTolerance = 1e-7
  MaxSccIterations = 200
}}

Analysis {{
  CalculateForces = Yes
}}
"""
    with open(out_path, 'w') as f:
        f.write(hsd)

def run_dftb_single_point(work_dir, enames, apos, sk_prefix):
    os.makedirs(work_dir, exist_ok=True)
    xyz_path = os.path.join(work_dir, "geom.xyz")
    hsd_path = os.path.join(work_dir, "dftb_in.hsd")
    write_xyz(xyz_path, enames, apos)
    write_dftb_input(enames, xyz_path, hsd_path, sk_prefix)

    cwd = os.getcwd()
    os.chdir(work_dir)
    try:
        ret = os.system('dftb+ > OUT 2> ERR')
        if ret != 0:
            raise RuntimeError(f"DFTB+ failed in {work_dir}")

        energy = None
        with open('OUT', 'r') as f:
            for line in f:
                if "Total Energy" in line:
                    energy = float(line[51:70].strip())
        if energy is None:
            raise RuntimeError(f"Could not parse energy from {work_dir}/OUT")
        return energy
    finally:
        os.chdir(cwd)

# ── Main z-scan ───────────────────────────────────────────────────────────
def run_zscan_for_atom(target_idx, pen_names, pen_pos, co_names, sk_prefix, z_distances, out_dir):
    """Run DFTB z-scan for a single target atom. Returns z_vals, e_vals."""
    target_name = pen_names[target_idx]
    target_pos = pen_pos[target_idx]
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
        combined_names = list(pen_names) + list(co_names)
        for iz, z in enumerate(z_distances):
            print(f"\n[z-scan {iz+1}/{len(z_distances)}] z = {z:.2f} Å")
            o_pos = np.array([target_pos[0], target_pos[1], target_pos[2] + z])
            c_pos = np.array([target_pos[0], target_pos[1], target_pos[2] + z + 1.13])
            co_pos_shifted = np.array([o_pos, c_pos])
            combined_pos = np.vstack([pen_pos, co_pos_shifted])

            work_dir = os.path.join(atom_dir, f'zscan_z{z:.2f}')
            t_start = time.time()
            try:
                energy_ha = run_dftb_single_point(work_dir, combined_names, combined_pos, sk_prefix)
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
    parser.add_argument('--output_dir', type=str, default=None)
    parser.add_argument('--target_indices', type=str, default='0',
                        help='Comma-separated list of target atom indices (e.g., "0,1,20,21")')
    parser.add_argument('--z_min', type=float, default=2.0, help='Minimum z-distance [Å]')
    parser.add_argument('--z_max', type=float, default=10.0, help='Maximum z-distance [Å]')
    parser.add_argument('--z_step', type=float, default=0.15, help='z step [Å]')
    args = parser.parse_args()

    basis = args.basis
    sk_prefix = SLAKO_PREFIXES[basis]
    out_dir = os.path.join(_THIS_DIR, args.output_dir) if args.output_dir else os.path.join(_THIS_DIR, f'zscan_{basis.replace("-", "_")}')
    os.makedirs(out_dir, exist_ok=True)

    target_indices = [int(x.strip()) for x in args.target_indices.split(',')]

    print(f"="*60)
    print(f"DFTB Z-Scan: basis={basis}, atoms={target_indices}")
    print(f"Output: {out_dir}")
    print(f"="*60)

    pen_names, pen_pos = load_xyz(_PENTACENE_XYZ)
    co_names,   _       = load_xyz(_CO_XYZ)
    print(f"Pentacene: {len(pen_names)} atoms")
    print(f"CO tip:    {len(co_names)} atoms")

    z_distances = np.arange(args.z_min, args.z_max + args.z_step*0.5, args.z_step)
    print(f"Z-scan: {len(z_distances)} points from {z_distances.min():.2f} to {z_distances.max():.2f} Å")

    for target_idx in target_indices:
        run_zscan_for_atom(target_idx, pen_names, pen_pos, co_names, sk_prefix, z_distances, out_dir)

    print(f"\nAll atoms done. Outputs in: {out_dir}/")

if __name__ == "__main__":
    main()
