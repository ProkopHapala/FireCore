#!/usr/bin/env python3
"""
vib_utils.py - Shared utilities for vibration spectra computation via ASE.

Handles:
 - Structure setup from PubChem or XYZ
 - Geometry caching: relaxed XYZ saved per method/basis to avoid recomputation
 - Vibration calculation dispatch for DFTB+ and PySCF backends
 - Frequency loading for plotting
"""

import numpy as np
import os
from pathlib import Path
from ase.io import read, write
from ase.optimize import BFGS
from ase.vibrations import Vibrations, Infrared

SK_PATHS = {
    'mio-1-1':   '/home/prokop/SIMULATIONS/dftbplus/slakos/mio-1-1/',
    'matsci-0-3': '/home/prokop/SIMULATIONS/dftbplus/slakos/matsci-0-3/',
    '3ob-3-1':   '/home/prokop/SIMULATIONS/dftbplus/slakos/3ob-3-1/',
    'pbc-0-3':   '/home/prokop/SIMULATIONS/dftbplus/slakos/pbc-0-3/',
}

# Max angular momentum for DFTB+
MAX_ANG = {'H': 's', 'C': 'p', 'N': 'p', 'O': 'p', 'S': 'p', 'P': 'p', 'Si': 'd', 'Ge': 'd'}


# ============================================================
# Structure setup
# ============================================================

def get_structure(name, cache_dir='.'):
    """
    Get molecule by name. Loads from local XYZ cache if available, else fetches from PubChem and caches.
    Cache file: <cache_dir>/<name>_pubchem.xyz
    """
    from ase.io import read, write
    xyz_cache = Path(cache_dir) / f'{name}_pubchem.xyz'
    if xyz_cache.exists():
        print(f"  [cache] Loading structure from {xyz_cache}")
        return read(str(xyz_cache))
    from ase.data.pubchem import pubchem_atoms_search
    print(f"Fetching '{name}' from PubChem...")
    atoms = pubchem_atoms_search(name)
    write(str(xyz_cache), atoms)
    print(f"  Got {len(atoms)} atoms: {set(atoms.get_chemical_symbols())}. Cached to {xyz_cache}")
    return atoms


def make_sila(adamantane):
    """C->Si substitution + scale positions for longer Si-Si bonds."""
    sila = adamantane.copy()
    for atom in sila:
        if atom.symbol == 'C':
            atom.symbol = 'Si'
    sila.positions *= 1.5
    return sila


# ============================================================
# Calculator factory
# ============================================================

def make_dftb_calc(atoms, sk_set='mio-1-1'):
    """Build ASE Dftb calculator for given atom types and SK set."""
    from ase.calculators.dftb import Dftb
    assert sk_set in SK_PATHS, f"Unknown SK set '{sk_set}'. Available: {list(SK_PATHS)}"
    sk_path = SK_PATHS[sk_set]
    
    atom_types = set(atoms.get_chemical_symbols())
    params = {
        'Hamiltonian_': 'DFTB',
        'Hamiltonian_SlaterKosterFiles_Type': 'Directory',
        'Hamiltonian_SlaterKosterFiles_Prefix': sk_path,
        'Hamiltonian_MaxAngularMomentum_': '',
        'Options_': '',
        'Options_WriteCharges': 'Yes',
        'Options_WriteBandOut': 'No',
    }
    for sym in atom_types:
        assert sym in MAX_ANG, f"No MaxAngularMomentum defined for element {sym}"
        params[f'Hamiltonian_MaxAngularMomentum_{sym}'] = MAX_ANG[sym]
    
    return Dftb(**params)


def make_pyscf_calc(atoms, method='b3lyp', basis='sto-3g'):
    """Build PySCF mean-field object for given atoms."""
    from pyscf import gto, dft, scf
    syms = atoms.get_chemical_symbols()
    pos  = atoms.get_positions()
    atom_str = '; '.join(f'{s} {x:.6f} {y:.6f} {z:.6f}' for s, (x, y, z) in zip(syms, pos))
    mol = gto.M(atom=atom_str, basis=basis, charge=0, spin=0, verbose=0)
    if method.lower() == 'hf':
        mf = mol.RHF()
    else:
        mf = dft.RKS(mol); mf.xc = method
    return mol, mf


def make_psi4_calc(atoms, method='b3lyp', basis='cc-pvdz'):
    """Build ASE Psi4 calculator for given atoms."""
    from ase.calculators.psi4 import Psi4
    return Psi4(atoms=atoms, method=method, basis=basis, charge=0, multiplicity=1)


def make_gpaw_calc(atoms, mode='lcao', basis='dzp', xc='PBE'):
    """Build ASE GPAW calculator for given atoms."""
    from gpaw import GPAW
    # GPAW requires periodic boundary conditions - add vacuum box for molecules
    if atoms.pbc.sum() == 0:
        # Add 10 Å vacuum in all directions
        atoms.center(vacuum=10.0)
        atoms.pbc = True
    return GPAW(mode=mode, basis=basis, xc=xc, txt='gpaw.txt')


def make_cp2k_calc(atoms, xc='PBE', basis_set='SZV-MOLOPT-SR-GTH'):
    """Build ASE CP2K calculator for given atoms."""
    from ase.calculators.cp2k import CP2K
    # CP2K requires periodic boundary conditions - add vacuum box for molecules
    if atoms.pbc.sum() == 0:
        # Add 10 Å vacuum in all directions
        atoms.center(vacuum=10.0)
        atoms.pbc = True
    return CP2K(
        xc=xc,
        basis_set=basis_set,
        basis_set_file='BASIS_MOLOPT',
        potential_file='POTENTIAL',
        command='cp2k_shell',
        label='cp2k',
        cutoff=300,
        max_scf=200,
        inp='''
&FORCE_EVAL
  &DFT
    &SCF
      EPS_SCF 1.0E-4
      &OT
        MINIMIZER DIIS
      &END OT
      &OUTER_SCF
        MAX_SCF 10
        EPS_SCF 1.0E-4
      &END OUTER_SCF
    &END SCF
  &END DFT
&END FORCE_EVAL
'''
    )


# ============================================================
# Geometry optimization with caching
# ============================================================

def relaxed_xyz_path(mol_name, method_tag, workdir='.'):
    """Return path for cached relaxed geometry XYZ file."""
    return Path(workdir) / f'{mol_name}_{method_tag}_relaxed.xyz'


def optimize_and_cache(atoms, mol_name, method_tag, fmax=0.001, steps=500, workdir='.'):
    """
    Optimize geometry. Load from cache if available, else run and save.

    Returns optimized ASE Atoms.
    """
    xyz_path = relaxed_xyz_path(mol_name, method_tag, workdir)
    
    if xyz_path.exists():
        print(f"  [cache] Loading relaxed geometry from {xyz_path}")
        return read(str(xyz_path))
    
    print(f"  [opt] Optimizing {mol_name} / {method_tag} (fmax={fmax})...")
    traj_path = Path(workdir) / f'{mol_name}_{method_tag}_opt.traj'
    opt = BFGS(atoms, trajectory=str(traj_path))
    opt.run(fmax=fmax, steps=steps)
    
    write(str(xyz_path), atoms)
    print(f"  [opt] Done. Saved to {xyz_path}")
    return atoms


# ============================================================
# Vibration calculation with caching
# ============================================================

def vib_cache_dir(mol_name, method_tag, workdir='.'):
    """Return name prefix for ASE Vibrations cache directory."""
    return str(Path(workdir) / f'{mol_name}_{method_tag}_vib')


def freq_npy_path(mol_name, method_tag, workdir='.'):
    """Path for saved frequency numpy array."""
    return Path(workdir) / f'{mol_name}_{method_tag}_freq.npy'


def hessian_npy_path(mol_name, method_tag, workdir='.'):
    """Path for saved Hessian numpy array (n_atoms x 3 x n_atoms x 3)."""
    return Path(workdir) / f'{mol_name}_{method_tag}_hessian.npy'


def compute_or_load_vibrations(atoms, mol_name, method_tag, workdir='.', with_ir=False):
    """
    Compute vibrational frequencies. Load from cache if already computed.

    Returns: complex frequency array in cm^-1
    """
    npy_path = freq_npy_path(mol_name, method_tag, workdir)
    
    if npy_path.exists():
        print(f"  [cache] Loading frequencies from {npy_path}")
        return np.load(str(npy_path))
    
    print(f"  [vib] Computing vibrations for {mol_name} / {method_tag}...")
    name_prefix = vib_cache_dir(mol_name, method_tag, workdir)
    vib = Vibrations(atoms, name=name_prefix)
    vib.run()
    vib.summary()
    
    freqs = vib.get_frequencies()
    np.save(str(npy_path), freqs)
    print(f"  [vib] Done. Saved to {npy_path}")
    
    if with_ir:
        try:
            ir_prefix = str(Path(workdir) / f'{mol_name}_{method_tag}_ir')
            ir = Infrared(atoms, name=ir_prefix)
            ir.run()
            ir.summary()
        except Exception as e:
            print(f"  [ir] IR failed: {e}")
    
    return freqs


# ============================================================
# Hessian caching and export
# ============================================================

def compute_or_load_hessian_pyscf(mol, mf, mol_name, method_tag, workdir='.'):
    """
    Compute Hessian matrix for PySCF.
    Load from cache if already computed.
    Returns: Hessian array (n_atoms x 3 x n_atoms x 3) in Hartree/Bohr^2
    """
    hess_path = hessian_npy_path(mol_name, method_tag, workdir)
    if hess_path.exists():
        print(f"  [cache] Loading Hessian from {hess_path}")
        return np.load(str(hess_path))

    print(f"  [hess] Computing Hessian for {mol_name}/{method_tag}...")
    hess_obj = mf.Hessian()
    hess = hess_obj.kernel()
    np.save(str(hess_path), hess)
    print(f"  [hess] Done. Saved to {hess_path}")
    return hess


def plot_hessians(mol_name, method_tags, workdir='.', figsize=(14, 5)):
    """
    Plot Hessian matrices side-by-side for comparison using imshow.

    Args:
        mol_name: Molecule name
        method_tags: List of method tags to compare (e.g., ['dftb_mio-1-1', 'pyscf_hf_sto-3g'])
        workdir: Directory with cached Hessian files
        figsize: Figure size for the plot
    """
    import matplotlib.pyplot as plt

    n_methods = len(method_tags)
    if n_methods == 0:
        print("No methods provided for comparison")
        return None

    fig, axes = plt.subplots(1, n_methods, figsize=figsize)
    if n_methods == 1:
        axes = [axes]

    for ax, tag in zip(axes, method_tags):
        hess_path = hessian_npy_path(mol_name, tag, workdir)
        if not hess_path.exists():
            print(f"  Missing Hessian: {hess_path}")
            ax.set_title(f'{tag}\n(NOT FOUND)')
            ax.axis('off')
            continue

        hess = np.load(str(hess_path))
        # Reshape to 2D for visualization: (n_atoms*3, n_atoms*3)
        n = hess.shape[0] * hess.shape[1]
        hess_2d = hess.reshape(n, n)

        # Use symmetric log scale for better visualization
        vmax = np.abs(hess_2d).max()
        im = ax.imshow(hess_2d, cmap='RdBu_r', vmin=-vmax, vmax=vmax, aspect='equal')
        ax.set_title(f'{tag}\nRange: [{hess_2d.min():.2e}, {hess_2d.max():.2e}]')
        ax.set_xlabel('Atom*3 + coord')
        ax.set_ylabel('Atom*3 + coord')
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    plt.suptitle(f'Hessian Matrix Comparison - {mol_name}', fontsize=14)
    plt.tight_layout()

    out_path = Path(workdir) / f'{mol_name}_hessian_comparison.png'
    fig.savefig(str(out_path), dpi=150)
    print(f"Saved: {out_path}")
    return fig


# ============================================================
# DFTB+ full pipeline
# ============================================================

def run_dftb(mol_name, atoms_in, sk_set='mio-1-1', fmax=0.001, workdir='.', with_ir=False):
    """
    Full DFTB+ pipeline: optimize (cached) -> vibrations (cached).
    Returns frequencies in cm^-1.
    """
    method_tag = f'dftb_{sk_set}'
    atoms = atoms_in.copy()
    atoms.calc = make_dftb_calc(atoms, sk_set)
    
    atoms = optimize_and_cache(atoms, mol_name, method_tag, fmax=fmax, workdir=workdir)
    atoms.calc = make_dftb_calc(atoms, sk_set)  # reattach after read from xyz
    
    freqs = compute_or_load_vibrations(atoms, mol_name, method_tag, workdir=workdir, with_ir=with_ir)
    return freqs, method_tag


# ============================================================
# Psi4 full pipeline
# ============================================================

def run_psi4(mol_name, atoms_in, method='b3lyp', basis='cc-pvdz', fmax=0.001, workdir='.'):
    """
    Full Psi4 pipeline: optimize (cached) -> vibrations (cached).
    Returns frequencies in cm^-1.
    """
    method_tag = f'psi4_{method}_{basis}'.replace('(', '').replace(')', '').replace('*', 's').replace(',', '')
    atoms = atoms_in.copy()
    atoms.calc = make_psi4_calc(atoms, method=method, basis=basis)

    atoms = optimize_and_cache(atoms, mol_name, method_tag, fmax=fmax, workdir=workdir)
    atoms.calc = make_psi4_calc(atoms, method=method, basis=basis)  # reattach after read from xyz

    freqs = compute_or_load_vibrations(atoms, mol_name, method_tag, workdir=workdir, with_ir=False)
    return freqs, method_tag


# ============================================================
# GPAW full pipeline
# ============================================================

def run_gpaw(mol_name, atoms_in, mode='lcao', basis='dzp', xc='PBE', fmax=0.001, workdir='.'):
    """
    Full GPAW pipeline: optimize (cached) -> vibrations (cached).
    Returns frequencies in cm^-1.
    """
    method_tag = f'gpaw_{mode}_{basis}_{xc}'.lower()
    atoms = atoms_in.copy()
    atoms.calc = make_gpaw_calc(atoms, mode=mode, basis=basis, xc=xc)

    atoms = optimize_and_cache(atoms, mol_name, method_tag, fmax=fmax, workdir=workdir)
    atoms.calc = make_gpaw_calc(atoms, mode=mode, basis=basis, xc=xc)  # reattach after read from xyz

    freqs = compute_or_load_vibrations(atoms, mol_name, method_tag, workdir=workdir, with_ir=False)
    return freqs, method_tag


# ============================================================
# CP2K full pipeline
# ============================================================

def run_cp2k(mol_name, atoms_in, xc='PBE', basis_set='DZVP-MOLOPT-SR-GTH', fmax=0.01, workdir='.'):
    """
    Full CP2K pipeline: optimize (cached) -> vibrations (cached).
    Returns frequencies in cm^-1.
    Uses looser fmax for CP2K due to computational cost.
    """
    method_tag = f'cp2k_{xc}_{basis_set}'.replace('-', '_').lower()
    atoms = atoms_in.copy()
    atoms.calc = make_cp2k_calc(atoms, xc=xc, basis_set=basis_set)

    atoms = optimize_and_cache(atoms, mol_name, method_tag, fmax=fmax, workdir=workdir)
    atoms.calc = make_cp2k_calc(atoms, xc=xc, basis_set=basis_set)  # reattach after read from xyz

    freqs = compute_or_load_vibrations(atoms, mol_name, method_tag, workdir=workdir, with_ir=False)
    return freqs, method_tag


# ============================================================
# PySCF full pipeline
# ============================================================

def run_pyscf(mol_name, atoms_in, method='b3lyp', basis='sto-3g', fmax=0.001, workdir='.'):
    """
    Full PySCF pipeline: optimize (cached) -> vibrations (cached).
    Returns frequencies in cm^-1.
    """
    from pyscf.geomopt import optimize as pyscf_optimize
    from pyscf.hessian import thermo
    from pyscf import gto, dft
    
    method_tag = f'pyscf_{method}_{basis}'.replace('(', '').replace(')', '').replace('*', 's').replace(',', '')
    xyz_path = relaxed_xyz_path(mol_name, method_tag, workdir)
    npy_path = freq_npy_path(mol_name, method_tag, workdir)
    
    if npy_path.exists():
        print(f"  [cache] Loading PySCF frequencies from {npy_path}")
        return np.load(str(npy_path)), method_tag
    
    # Build PySCF molecule
    atoms = atoms_in.copy()
    mol, mf = make_pyscf_calc(atoms, method=method, basis=basis)
    
    # Optimize (or load)
    if xyz_path.exists():
        print(f"  [cache] Loading PySCF relaxed geometry from {xyz_path}")
        atoms_opt = read(str(xyz_path))
        syms = atoms_opt.get_chemical_symbols()
        pos  = atoms_opt.get_positions()
        atom_str = '; '.join(f'{s} {x:.6f} {y:.6f} {z:.6f}' for s, (x, y, z) in zip(syms, pos))
        mol = gto.M(atom=atom_str, basis=basis, charge=0, spin=0, verbose=0)
        if method.lower() == 'hf':
            mf = mol.RHF()
        else:
            mf = dft.RKS(mol); mf.xc = method
    else:
        print(f"  [opt] Optimizing {mol_name} via PySCF {method}/{basis}...")
        mf.kernel()
        mol_eq = pyscf_optimize(mf)
        # Save relaxed geometry
        coords_bohr = mol_eq.atom_coords()
        from ase import Atoms as ASEAtoms
        from ase.units import Bohr
        syms = [mol_eq.atom_symbol(i) for i in range(mol_eq.natm)]
        pos  = coords_bohr * Bohr
        atoms_opt = ASEAtoms(symbols=syms, positions=pos)
        write(str(xyz_path), atoms_opt)
        # Rebuild mf at optimized geometry
        atom_str = '; '.join(f'{s} {x:.6f} {y:.6f} {z:.6f}' for s, (x, y, z) in zip(syms, pos))
        mol = gto.M(atom=atom_str, basis=basis, charge=0, spin=0, verbose=0)
        if method.lower() == 'hf':
            mf = mol.RHF()
        else:
            mf = dft.RKS(mol); mf.xc = method
    
    # Hessian + frequencies
    print(f"  [vib] Computing PySCF Hessian {mol_name} / {method_tag}...")
    mf.kernel()
    hess_obj = mf.Hessian()
    hess = hess_obj.kernel()
    freq_result = thermo.harmonic_analysis(mol, hess)
    freqs = freq_result['freq_wavenumber']
    np.save(str(npy_path), freqs)
    print(f"  [vib] Done. Saved to {npy_path}")
    
    return freqs, method_tag


# ============================================================
# PySCF full pipeline with Hessian + eigenmodes export
# ============================================================

def vib_npz_path(mol_name, method_tag, workdir='.'):
    """Path for saved vibration NPZ (hessian, frequencies, eigenmodes, geometry)."""
    return Path(workdir) / f'{mol_name}_{method_tag}_vib.npz'


def run_pyscf_vib_full(mol_name, atoms_in, method='b3lyp', basis='sto-3g', fmax=0.001, workdir='.'):
    """
    Full PySCF pipeline: optimize (cached) -> Hessian -> harmonic analysis.
    Saves hessian, frequencies, eigenmodes, geometry to a single .npz file.

    Returns dict with keys: freqs, hessian, modes, symbols, positions, method_tag.
    """
    from pyscf.geomopt import optimize as pyscf_optimize
    from pyscf.hessian import thermo
    from pyscf import gto, dft, scf

    method_tag = f'pyscf_{method}_{basis}'.replace('(', '').replace(')', '').replace('*', 's').replace(',', '')
    xyz_path = relaxed_xyz_path(mol_name, method_tag, workdir)
    npz_path = vib_npz_path(mol_name, method_tag, workdir)

    if npz_path.exists():
        print(f"  [cache] Loading PySCF vibration data from {npz_path}")
        data = np.load(str(npz_path), allow_pickle=True)
        return {k: data[k] for k in data.files}, method_tag

    # Build PySCF molecule
    atoms = atoms_in.copy()
    mol, mf = make_pyscf_calc(atoms, method=method, basis=basis)

    # Optimize (or load from cache)
    if xyz_path.exists():
        print(f"  [cache] Loading PySCF relaxed geometry from {xyz_path}")
        atoms_opt = read(str(xyz_path))
    else:
        print(f"  [opt] Optimizing {mol_name} via PySCF {method}/{basis}...")
        mf.kernel()
        mol_eq = pyscf_optimize(mf)
        from ase import Atoms as ASEAtoms
        from ase.units import Bohr
        syms = [mol_eq.atom_symbol(i) for i in range(mol_eq.natm)]
        pos = mol_eq.atom_coords() * Bohr
        atoms_opt = ASEAtoms(symbols=syms, positions=pos)
        write(str(xyz_path), atoms_opt)

    # Rebuild mol/mf at optimized geometry
    syms = atoms_opt.get_chemical_symbols()
    pos = atoms_opt.get_positions()
    atom_str = '; '.join(f'{s} {x:.6f} {y:.6f} {z:.6f}' for s, (x, y, z) in zip(syms, pos))
    mol = gto.M(atom=atom_str, basis=basis, charge=0, spin=0, verbose=0)
    if method.lower() == 'hf':
        mf = mol.RHF()
    else:
        mf = dft.RKS(mol); mf.xc = method

    # SCF + Hessian
    print(f"  [scf] Running SCF for {mol_name} / {method_tag}...")
    mf.kernel()
    print(f"  [hess] Computing PySCF Hessian {mol_name} / {method_tag}...")
    hess_obj = mf.Hessian()
    hess = hess_obj.kernel()  # shape (natm, 3, natm, 3) in Hartree/Bohr^2

    # Harmonic analysis -> frequencies + normal modes
    print(f"  [vib] Harmonic analysis...")
    freq_result = thermo.harmonic_analysis(mol, hess)
    freqs = freq_result['freq_wavenumber']  # (3N,) cm^-1
    # Normal mode matrix: (3N, 3N) — columns are mode vectors
    modes = freq_result.get('C_ij', None)
    if modes is None:
        # Fallback: some PySCF versions use different key
        modes = freq_result.get('mode', None)
    if modes is None:
        # Compute from hessian directly: mass-weighted eigenvectors
        from pyscf.data import elements
        mass = np.array([elements.masses[mol.atom_symbol(i)] for i in range(mol.natm)])
        mass_3n = np.repeat(mass, 3)
        hess_2d = hess.reshape(mol.natm * 3, mol.natm * 3)
        # Mass-weighted Hessian
        mw_hess = hess_2d / np.sqrt(mass_3n[:, None] * mass_3n[None, :])
        evals, evecs = np.linalg.eigh(mw_hess)
        modes = evecs  # (3N, 3N)

    # IR intensities if available
    ir_int = freq_result.get('ir_int', None)

    # Save everything to .npz
    save_dict = {
        'freqs': freqs,
        'hessian': hess,
        'modes': modes,
        'symbols': np.array(syms),
        'positions': pos,  # Angstrom
        'method': method,
        'basis': basis,
        'energy': float(mf.e_tot),
    }
    if ir_int is not None:
        save_dict['ir_int'] = ir_int

    np.savez(str(npz_path), **save_dict)
    print(f"  [vib] Saved to {npz_path}")
    print(f"    natoms={len(syms)}, nmodes={len(freqs)}, energy={mf.e_tot:.8f} Ha")

    # Also save freqs as .npy for compatibility with existing plot code
    npy_path = freq_npy_path(mol_name, method_tag, workdir)
    if not npy_path.exists():
        np.save(str(npy_path), freqs)

    return save_dict, method_tag


# ============================================================
# Frequency filtering
# ============================================================

def filter_real_freqs(freqs, threshold=10.0):
    """Return only real positive frequencies above threshold (cm^-1)."""
    if np.iscomplexobj(freqs):
        mask = (freqs.imag == 0) & (freqs.real > threshold)
        return freqs[mask].real
    return freqs[freqs > threshold]


def export_modes_to_xyz(atoms, mol_name, method_tag, workdir='.', threshold=10.0):
    """
    Export all vibration modes to a single multi-frame XYZ file with displacement vectors.
    Each frame represents one vibration mode with frequency in comment line.
    Extended XYZ format: atom_symbol x y z dx dy dz (displacement vector)

    Output: <mol>_<method>_all_modes.xyz (single multi-frame file)
    """
    from ase.io import write
    npy_path = freq_npy_path(mol_name, method_tag, workdir)
    if not npy_path.exists():
        print(f"  [export] No frequency file found: {npy_path}")
        return

    freqs = np.load(str(npy_path))
    real_freqs = filter_real_freqs(freqs, threshold=threshold)
    n_modes = len(real_freqs)
    print(f"  [export] Exporting {n_modes} modes for {mol_name}/{method_tag}")

    # Output single multi-frame XYZ file
    fname = Path(workdir) / f'{mol_name}_{method_tag}_all_modes.xyz'

    # For DFTB+, set up atoms with calculator and use cached Vibrations
    if method_tag.startswith('dftb_'):
        sk_set = method_tag.replace('dftb_', '')
        # Load relaxed geometry
        xyz_path = relaxed_xyz_path(mol_name, method_tag, workdir)
        if xyz_path.exists():
            atoms_calc = read(str(xyz_path))
        else:
            atoms_calc = atoms.copy()
        atoms_calc.calc = make_dftb_calc(atoms_calc, sk_set)
        # Use legacy cache name if it exists, otherwise use new naming
        legacy_cache = str(Path(workdir) / f'{mol_name}_dftb_vib')
        new_cache = vib_cache_dir(mol_name, method_tag, workdir)
        cache_name = legacy_cache if Path(legacy_cache).exists() else new_cache
        # Create Vibrations object with cache name
        vib = Vibrations(atoms_calc, name=cache_name)
        # Write all modes to single multi-frame XYZ file
        with open(fname, 'w') as fp:
            for i, f in enumerate(freqs):
                if f.imag == 0 and f.real > threshold:
                    # Get displacement vectors for this mode
                    mode = vib.get_mode(i)
                    # Write frame header
                    fp.write(f'{len(atoms_calc)}\n')
                    fp.write(f'Freq: {f.real:.2f} cm-1\n')
                    # Write atoms with displacement vectors
                    for j, atom in enumerate(atoms_calc):
                        x, y, z = atom.position
                        dx, dy, dz = mode[j]
                        fp.write(f'{atom.symbol:2s} {x:12.6f} {y:12.6f} {z:12.6f} {dx:12.6f} {dy:12.6f} {dz:12.6f}\n')
        print(f"  [export] Done: {fname} ({n_modes} modes)")
        return

    # PySCF: load from .npz saved by run_pyscf_vib_full
    npz_path = vib_npz_path(mol_name, method_tag, workdir)
    if not npz_path.exists():
        print(f"  [export] No NPZ file found: {npz_path}")
        print(f"  [export] Run run_pyscf_vib_full() first to generate vibration data")
        return

    data = np.load(str(npz_path), allow_pickle=True)
    symbols = data['symbols']
    positions = data['positions']
    modes = data['modes']  # (3N, 3N) — columns are mode vectors
    freqs_npz = data['freqs']
    natoms = len(symbols)

    with open(fname, 'w') as fp:
        for i, f in enumerate(freqs_npz):
            f_real = f.real if np.iscomplexobj(freqs_npz) else f
            if f_real > threshold:
                mode_vec = modes[:, i].reshape(natoms, 3)
                fp.write(f'{natoms}\n')
                fp.write(f'Freq: {f_real:.2f} cm-1\n')
                for j in range(natoms):
                    x, y, z = positions[j]
                    dx, dy, dz = mode_vec[j]
                    fp.write(f'{symbols[j]:2s} {x:12.6f} {y:12.6f} {z:12.6f} {dx:12.6f} {dy:12.6f} {dz:12.6f}\n')
    print(f"  [export] Done: {fname}")
