#!/usr/bin/env python3
"""List available basis sets in PySCF."""

import os
import pyscf
from pyscf import gto
from pyscf.gto.basis import parse_nwchem

# --- 1. Basis set .dat files shipped with PySCF ---
basis_dir = os.path.join(os.path.dirname(pyscf.__file__), "gto", "basis")
dat_files = sorted([f for f in os.listdir(basis_dir) if f.endswith(".dat")])
print("=== Built-in .dat basis sets (first 60) ===")
for f in dat_files[:60]:
    print(f"  {f}")
print(f"  ... total {len(dat_files)} .dat files\n")

# --- 2. Minimal-basis to large-basis quick reference ---
quick_ref = [
    ("STO-3G"     , "minimal (single-zeta, 3-primitive Gaussian per orbital)"),
    ("3-21G"      , "split-valence minimal"),
    ("6-31G"      , "split-valence Pople"),
    ("6-31G(d)"   , "split-valence + polarization on heavy atoms"),
    ("6-31G(d,p)" , "split-valence + polarization on all atoms"),
    ("6-311G(d,p)", "triple-zeta Pople"),
    ("cc-pVDZ"    , "Dunning correlation-consistent double-zeta"),
    ("cc-pVTZ"    , "Dunning correlation-consistent triple-zeta"),
    ("cc-pVQZ"    , "Dunning quadruple-zeta"),
    ("def2-SVP"   , "Karlsruhe split-valence + polarization"),
    ("def2-TZVP"  , "Karlsruhe triple-zeta valence + polarization"),
    ("def2-TZVPP" , "triple-zeta with extra polarization"),
    ("def2-QZVP"  , "quadruple-zeta"),
    ("aug-cc-pVTZ", "augmented triple-zeta (diffuse functions)"),
    ("pc-2"       , "Jensen polarization-consistent triple-zeta"),
    ("pc-3"       , "Jensen polarization-consistent quadruple-zeta"),
    ("lanl2dz"    , "effective core potential / pseudopotential basis"),
    ("cc-pwCVTZ"  , "weighted core-valence triple-zeta"),
]

print("=== Quick-reference commonly used basis sets ===")
for name, desc in quick_ref:
    print(f"  {name:<18s}  {desc}")

# --- 3. Verify a few exist by trying to build a molecule ---
print("\n=== Spot-check: can PySCF load these? ===")
for bs in ["sto-3g", "6-31g(d)", "cc-pvtz", "def2-tzvp", "lanl2dz"]:
    try:
        mol = gto.M(atom="O 0 0 0; H 0 0 1; H 0 1 0", basis=bs, verbose=0)
        print(f"  {bs:<15s}  OK  (nao={mol.nao})")
    except Exception as e:
        print(f"  {bs:<15s}  FAIL: {e}")
