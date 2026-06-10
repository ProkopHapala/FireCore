#!/usr/bin/env python3
import os, sys, argparse, tempfile
import numpy as np

sys.path.append("../../")
from pyBall import MMFF

_MMFF_BOND_K = {
    ("H","H",1): 10.000,
    ("C","C",1): 10.082,
}

_MASS_AMU = {"H":1.00784, "C":12.0107}

def thz_from_eig_massweighted(eig_eVA2_per_amu):
    eV=1.602176634e-19
    Ang=1e-10
    amu=1.66053906660e-27
    conv = (eV/(Ang*Ang))/amu
    w = np.sign(eig_eVA2_per_amu)*np.sqrt(np.abs(eig_eVA2_per_amu)*conv)
    f = w/(2*np.pi)*1e-12
    return f

def make_mol_diatomic(sym, l0, fname, syma_label=None, symb_label=None):
    """Create MDL .mol file with explicit bond. Labels like 'C_3' force atom type."""
    syma,symb = sym
    la = syma_label if syma_label else syma
    lb = symb_label if symb_label else symb
    with open(fname, "w") as f:
        f.write("\n  Avogadro\n\n")
        f.write("  2  1  0  0  0  0  0  0  0  0999 V2000\n")
        f.write(f"    0.0000    0.0000    0.0000 {la:2s}  0  0  0  0  0  0  0  0  0  0  0  0\n")
        f.write(f"    {l0:7.4f}    0.0000    0.0000 {lb:2s}  0  0  0  0  0  0  0  0  0  0  0  0\n")
        f.write("  1  2  1  0  0  0  0\n")
        f.write("M  END\n")

def run_case(sym=("H","H"), bond_order=1, l0=0.740, bUFF=False, dx=1e-4, syma_label=None, symb_label=None):
    k = _MMFF_BOND_K.get((sym[0],sym[1],bond_order), None)
    if k is None: k = _MMFF_BOND_K.get((sym[1],sym[0],bond_order), None)
    if k is None: raise ValueError(f"No reference k for {sym} order={bond_order}")

    with tempfile.TemporaryDirectory() as td:
        mol = os.path.join(td, "diatomic.mol")
        make_mol_diatomic(sym, l0, mol, syma_label, symb_label)
        MMFF.init(xyz_name=mol, nPBC=(0,0,0), bEpairs=False, bMMFF=True, bUFF=bUFF)
        MMFF.setSwitches(PBC=-1, NonBonded=-1, Angles=-1, PiSigma=-1, PiPiI=-1)
        inds = np.arange(2, dtype=np.int32)
        H = MMFF.getHessian3Nx3N(inds, dx=dx)
        H = 0.5*(H+H.T)
        lam = np.linalg.eigvalsh(H)
        lam = np.sort(lam)

        m = np.array([_MASS_AMU[sym[0]], _MASS_AMU[sym[1]]], dtype=float)
        w = np.repeat(1.0/np.sqrt(m), 3)
        Hw = (H*w[None,:])*w[:,None]
        lamw = np.linalg.eigvalsh(0.5*(Hw+Hw.T))
        lamw = np.sort(lamw)
        fTHz = np.sort(thz_from_eig_massweighted(lamw))

        print("\n==============================")
        print(f"Diatomic {sym[0]}-{sym[1]} l0={l0} order={bond_order} bUFF={int(bUFF)} dx={dx}")
        print("==============================")
        print(f"reference bond k [eV/Å^2] = {k}")
        print(f"expected λ_stretch ≈ 4*k [eV/Å^2] = {4*k}")
        print("eig λ [eV/Å^2] sorted:", lam)
        print("eig μ [eV/Å^2/amu] sorted:", lamw)
        print("freq [THz] sorted:", fTHz)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument('--case', choices=['H2','C2','C2_sp3'], default='H2')
    ap.add_argument('--l0', type=float, default=None)
    ap.add_argument('--uff', action='store_true')
    ap.add_argument('--dx', type=float, default=1e-4)
    args = ap.parse_args()

    if args.case=='H2':
        l0 = 0.740 if args.l0 is None else args.l0
        run_case(("H","H"), bond_order=1, l0=l0, bUFF=args.uff, dx=args.dx)
    elif args.case=='C2':
        # Uses automatic atom type detection (will pick C_1 sp for short bond)
        l0 = 1.538 if args.l0 is None else args.l0
        run_case(("C","C"), bond_order=1, l0=l0, bUFF=args.uff, dx=args.dx)
    elif args.case=='C2_sp3':
        # Force C_3 type to test true C-C single bond from BondTypes.dat (k=10.082)
        l0 = 1.538 if args.l0 is None else args.l0
        run_case(("C","C"), bond_order=1, l0=l0, bUFF=args.uff, dx=args.dx, syma_label='C_3', symb_label='C_3')
