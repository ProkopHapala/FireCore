"""
DOS computation interface: Fireball SCF + spectral function with Γ broadening.

PURPOSE:
  Thin interface to pyBall.FireballOCL.STM.compute_dos() for running Fireball SCF
  and computing the spectral function (projected DOS) with lead coupling broadening.
  Intended to be run as a subprocess to avoid Fireball reinitialization issues.

PHYSICAL BACKGROUND:

1) Fireball SCF:
  Self-consistent field calculation using Fireball DFTB (Density Functional Tight Binding).
  Produces the Hamiltonian H and overlap S matrices in the LCAO (Linear Combination of
  Atomic Orbitals) basis with numerical atomic orbitals.

2) Spectral Function with Γ Broadening:
  The retarded Green's function for a system coupled to leads is:
    G^r(E) = [(E + iη)S - H - Σ]^{-1}
  where Σ = -iΓ is the wide-band limit self-energy applied to contact atoms.
  The spectral function (projected DOS) is:
    A(E) = (1/π) Im[G^r(E)]

  In the molecular orbital basis, the spectral function is a sum of Lorentzians:
    A_μν(E) = Σ_n C_μn C_νn * (Γ_n/π) / [(E - ε_n)² + Γ_n²]
  where:
    - ε_n: MO eigenvalues
    - C_μn: MO eigenvectors (coefficients in AO basis)
    - Γ_n: orbital-resolved broadening, Γ_n = γ Σ_{μ∈contact} C_μn²

3) Generalized Eigenvalue Problem:
  Fireball solves: H C = S C ε
  where C is the matrix of eigenvectors and ε is the diagonal matrix of eigenvalues.
  The overlap S is not identity due to non-orthogonal atomic orbitals.

4) Output Data:
  The computed .npz file contains:
    - E_grid: Energy grid around Fermi level
    - A: Spectral function A_μ(E) [nE, norb]
    - E_fermi: Fermi level
    - eigen: Fireball eigenvalues (for reference)
    - atomTypes, atomPos: Molecular geometry
    - H, S: Dense Hamiltonian and overlap matrices
    - orb2atom: Orbital-to-atom mapping
    - norb_per: Orbitals per atom
    - contact_orbs: Indices of contact orbitals
    - gamma: Broadening strength
    - eps: MO eigenvalues (from generalized eigensolve)
    - C: MO eigenvectors (from generalized eigensolve)

  Note: eps and C are required for MO overlap diagnostics.

USAGE EXAMPLES:

  Compute DOS for H2 tip with Γ=0.05 eV on lowest-z atom:
    cd tests/pyFireball
    python stm_compute_dos.py --mol H2 --out export/stm_phase1/dos_H2.npz --gamma 0.05 --contact lowest_z

  Compute DOS for PTCDA sample with Γ=0.1 eV on all atoms:
    python stm_compute_dos.py --xyz PTCDA.xyz --out export/stm_phase1/dos_PTCDA.npz --gamma 0.1 --contact all

  Compute DOS with custom energy window:
    python stm_compute_dos.py --xyz PTCDA.xyz --out dos_PTCDA.npz --E_range 8.0 --dE 0.01

KEY CLI OPTIONS:
  --mol NAME              Use built-in molecule (H2, CO, etc.)
  --xyz FILE              Load geometry from XYZ file
  --out PATH              Output .npz file path (required)
  --gamma VALUE           Broadening strength Γ (eV), default 0.1
  --contact {all,lowest_z,1,2,3}  Contact atoms for Γ broadening
  --nscf N                Maximum SCF iterations, default 200
  --E_range VALUE         Energy window around Fermi (eV), default 6.0
  --dE VALUE              Energy grid step (eV), default 0.02

OUTPUT:
  Saves .npz file with all spectral data, geometry, and matrices.
  Also prints diagnostic information (norb, eigenvalues, sum-rule check).

NOTES:
  - Must run from tests/pyFireball directory (requires Fdata symlink)
  - Fdata symlink is automatically created to tests/Fireball/Fdata_HCNOS
  - Energies are in eV
  - Positions are in Å
"""
import argparse, os, sys
import numpy as np

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.normpath(os.path.join(_THIS_DIR, "..", ".."))
if _REPO_ROOT not in sys.path: sys.path.insert(0, _REPO_ROOT)
os.chdir(_THIS_DIR)

# Fdata symlink
_FDATA_HCNOS = os.path.join(_REPO_ROOT, "tests", "Fireball", "Fdata_HCNOS")
_FDATA_LOCAL  = os.path.join(_THIS_DIR, "Fdata")
if os.path.realpath(_FDATA_LOCAL) if os.path.exists(_FDATA_LOCAL) else "" != os.path.realpath(_FDATA_HCNOS):
    if os.path.lexists(_FDATA_LOCAL): os.unlink(_FDATA_LOCAL)
    os.symlink(_FDATA_HCNOS, _FDATA_LOCAL)

from pyBall.FireballOCL.STM import compute_dos, load_xyz, make_mol

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--mol",     default=None)
    p.add_argument("--xyz",     default=None)
    p.add_argument("--out",     required=True)
    p.add_argument("--gamma",   type=float, default=0.1)
    p.add_argument("--contact", default="all")
    p.add_argument("--nscf",    type=int,   default=200)
    p.add_argument("--E_range", type=float, default=6.0)
    p.add_argument("--dE",      type=float, default=0.02)
    args = p.parse_args()

    if args.mol:
        atomTypes, atomPos = make_mol(args.mol)
    elif args.xyz:
        atomTypes, atomPos = load_xyz(args.xyz)
    else:
        raise ValueError("Need --mol or --xyz")

    print(f"Molecule: {len(atomTypes)} atoms {[{'H':1,'C':6,'N':7,'O':8,'S':16}.get(z,'?') for z in atomTypes]}")
    result = compute_dos(atomTypes, atomPos, gamma=args.gamma, contact=args.contact,
                         nmax_scf=args.nscf, E_range=args.E_range, dE=args.dE,
                         out_path=args.out)
    return result

if __name__=="__main__": main()
