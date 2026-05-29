#!/usr/bin/env python3
"""
run_vib_spectra.py - Compute vibrational spectra for a molecule using multiple methods.

Usage:
    python run_vib_spectra.py adamantane
    python run_vib_spectra.py sila_adamantane
    python run_vib_spectra.py adamantane --methods dftb_mio pyscf_hf

Results are cached on disk (relaxed geometry + frequencies) so re-runs skip costly steps.
"""

import argparse
import numpy as np
from vib_utils import (
    get_structure, make_sila,
    run_dftb, run_pyscf, run_psi4, run_gpaw, run_cp2k,
    freq_npy_path, export_modes_to_xyz,
    compute_or_load_hessian_pyscf,
    relaxed_xyz_path, make_pyscf_calc,
    plot_hessians
)

MOLECULES = {
    'adamantane':      lambda: get_structure('adamantane'),
    'sila_adamantane': lambda: make_sila(get_structure('adamantane')),
}

METHODS = {
    'dftb_mio':    lambda atoms, mol, wd: run_dftb(mol, atoms, sk_set='mio-1-1',    workdir=wd, with_ir=True),
    'dftb_3ob':    lambda atoms, mol, wd: run_dftb(mol, atoms, sk_set='3ob-3-1',    workdir=wd, with_ir=True),
    'dftb_matsci': lambda atoms, mol, wd: run_dftb(mol, atoms, sk_set='matsci-0-3', workdir=wd, with_ir=True),
    'dftb_pbc':    lambda atoms, mol, wd: run_dftb(mol, atoms, sk_set='pbc-0-3',    workdir=wd, with_ir=True),
    'cp2k_pbe':    lambda atoms, mol, wd: run_cp2k(mol, atoms, xc='PBE', basis_set='SZV-MOLOPT-SR-GTH', workdir=wd),
    'gpaw_pbe':    lambda atoms, mol, wd: run_gpaw(mol, atoms, mode='lcao', basis='dzp', xc='PBE', workdir=wd),
    'psi4_b3lyp':  lambda atoms, mol, wd: run_psi4(mol, atoms, method='b3lyp', basis='cc-pvdz', workdir=wd),
    'pyscf_hf':    lambda atoms, mol, wd: run_pyscf(mol, atoms, method='hf',    basis='sto-3g', workdir=wd),
    'pyscf_b3lyp': lambda atoms, mol, wd: run_pyscf(mol, atoms, method='b3lyp', basis='sto-3g', workdir=wd),
}

DEFAULT_METHODS = {
    'adamantane':      ['dftb_mio', 'cp2k_pbe'],
    'sila_adamantane': ['dftb_matsci', 'cp2k_pbe'],
}


def main():
    parser = argparse.ArgumentParser(description='Compute vibrational spectra for a molecule')
    parser.add_argument('molecule', choices=list(MOLECULES.keys()), help='Molecule to compute')
    parser.add_argument('--methods', nargs='+', choices=list(METHODS.keys()),
                        help='Methods to use (default: per-molecule defaults)')
    parser.add_argument('--workdir', default='.', help='Directory for cached files and output')
    parser.add_argument('--export-modes', action='store_true',
                        help='Export vibration modes to multi-frame XYZ with displacement vectors')
    parser.add_argument('--cache-hessians', action='store_true',
                        help='Compute and cache Hessian matrices (PySCF only)')
    parser.add_argument('--plot-hessians', action='store_true',
                        help='Plot cached Hessian matrices side-by-side (PySCF only)')
    args = parser.parse_args()

    mol_name = args.molecule
    methods  = args.methods or DEFAULT_METHODS[mol_name]
    workdir  = args.workdir

    # Plot only mode - no computation
    if args.plot_hessians and not args.cache_hessians:
        print(f"\n=== Plotting Hessian comparison for {mol_name} ===")
        # Collect method tags for available Hessians
        method_tags = []
        for mkey in methods:
            if mkey.startswith('dftb_'):
                sk_set = mkey.replace('dftb_', '') if '_' in mkey else 'mio-1-1'
                tag = f'dftb_{sk_set}'
            else:
                tag = f'{mkey}_sto-3g'.replace('pyscf_', 'pyscf_')
            method_tags.append(tag)
        plot_hessians(mol_name, method_tags, workdir=workdir)
        return

    print(f"\n=== {mol_name} | methods: {methods} ===\n")
    atoms_in = MOLECULES[mol_name]()

    results = {}
    for mkey in methods:
        print(f"\n--- Method: {mkey} ---")
        try:
            freqs, method_tag = METHODS[mkey](atoms_in, mol_name, workdir)
            results[mkey] = (freqs, method_tag)
            print(f"  {len(freqs)} modes, tag: {method_tag}")
        except Exception as e:
            print(f"  FAILED: {e}")
            continue

        # Compute and cache Hessian if requested (PySCF only - ASE Vibrations doesn't expose Hessian)
        if args.cache_hessians and mkey in results and not mkey.startswith('dftb_'):
            print(f"  [hess] Caching Hessian for {mkey}...")
            try:
                # PySCF - need to rebuild mol/mf at cached geometry
                from ase.io import read
                xyz_path = relaxed_xyz_path(mol_name, method_tag, workdir)
                if xyz_path.exists():
                    atoms_opt = read(str(xyz_path))
                else:
                    atoms_opt = atoms_in.copy()
                method = 'hf' if 'hf' in mkey else 'b3lyp'
                mol, mf = make_pyscf_calc(atoms_opt, method=method, basis='sto-3g')
                mf.kernel()
                compute_or_load_hessian_pyscf(mol, mf, mol_name, method_tag, workdir=workdir)
            except Exception as e:
                print(f"    FAILED: {e}")

    print(f"\n=== Finished: {mol_name} ===")
    print("Cached frequency files:")
    for mkey, (freqs, tag) in results.items():
        p = freq_npy_path(mol_name, tag, workdir)
        print(f"  {mkey:15s}  ->  {p}")

    if args.export_modes:
        print(f"\n=== Exporting vibration modes ===")
        for mkey, (freqs, tag) in results.items():
            print(f"\n--- {mkey} ({tag}) ---")
            try:
                export_modes_to_xyz(atoms_in, mol_name, tag, workdir=workdir)
            except Exception as e:
                print(f"  FAILED: {e}")

    if args.cache_hessians and results:
        print(f"\n=== Plotting Hessian comparison ===")
        method_tags = [tag for _, tag in results.values()]
        plot_hessians(mol_name, method_tags, workdir=workdir)


if __name__ == '__main__':
    main()
