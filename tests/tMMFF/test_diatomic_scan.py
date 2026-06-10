#!/usr/bin/env python3
"""Scan energy/force along bond direction for diatomic molecule."""
import os, sys, argparse, tempfile
import numpy as np

sys.path.append("../../")
from pyBall import MMFF

_MASS_AMU = {"H":1.00784, "C":12.0107}

def make_mol_diatomic(sym, l0, fname, syma_label=None, symb_label=None):
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

def run_scan(sym=("H","H"), bond_order=1, l0=0.74, l_range=0.3, nstep=61, syma_label=None, symb_label=None):
    with tempfile.TemporaryDirectory() as td:
        mol = os.path.join(td, "diatomic.mol")
        make_mol_diatomic(sym, l0, mol, syma_label, symb_label)
        MMFF.init(xyz_name=mol, nPBC=(0,0,0), bEpairs=False, bMMFF=True, bUFF=False)
        # Disable everything except bonds
        MMFF.setSwitches(PBC=-1, NonBonded=-1, Angles=-1, PiSigma=-1, PiPiI=-1)
        
        # Get bond params
        MMFF.print_debugs(bParams=True, bNeighs=True)
        
        poss = []
        Escan = []
        
        for i in range(nstep):
            # Reset positions to equilibrium first
            MMFF.init(xyz_name=mol, nPBC=(0,0,0), bEpairs=False, bMMFF=True, bUFF=False)
            MMFF.setSwitches(PBC=-1, NonBonded=-1, Angles=-1, PiSigma=-1, PiPiI=-1)
            
            disp = -l_range/2 + (i/(nstep-1))*l_range
            # Displace atom 1 along x
            sel = np.array([1], dtype=np.int32)
            shift = np.array([disp, 0.0, 0.0])
            MMFF.shift_atoms_ax(shift, sel=sel)
            E = MMFF.eval()
            poss.append(l0 + disp)
            Escan.append(E)
        
        poss = np.array(poss)
        Escan = np.array(Escan)
        
        # Find minimum
        imin = np.argmin(Escan)
        l_min = poss[imin]
        E_min = Escan[imin]
        
        print(f"\n{'='*50}")
        print(f"Diatomic {sym[0]}-{sym[1]} bond scan")
        print(f"{'='*50}")
        print(f"Equilibrium bond length from scan: {l_min:.4f} Å")
        print(f"Minimum energy: {E_min:.6f} eV")
        print(f"\nBond length (Å) | Energy (eV)")
        print("-" * 35)
        for i in range(0, nstep, 5):
            print(f"{poss[i]:13.4f} | {Escan[i]:11.6f}")
        
        # Quadratic fit around minimum
        mask = np.abs(poss - l_min) < 0.1  # Fit within ±0.1 Å of minimum
        if np.sum(mask) >= 3:
            coeffs = np.polyfit(poss[mask] - l_min, Escan[mask], 2)
            k_fit = 2 * coeffs[0]  # E = 1/2*k*(Δx)² => d²E/dx² = k
            print(f"\nQuadratic fit around minimum (±0.1 Å):")
            print(f"  d²E/dx² = {k_fit:.3f} eV/Å²")
            print(f"  Compare to bond table k (expected ~10 eV/Å² for C-C)")
        
        return poss, Escan

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument('--case', choices=['H2','C2','C2_sp3','C2_R'], default='C2_sp3')
    ap.add_argument('--l0', type=float, default=None)
    ap.add_argument('--range', type=float, default=0.3, help='Scan range (Å)')
    ap.add_argument('--nstep', type=int, default=61)
    args = ap.parse_args()
    
    if args.case == 'H2':
        l0 = 0.74 if args.l0 is None else args.l0
        run_scan(("H","H"), 1, l0, args.range, args.nstep)
    elif args.case == 'C2_sp3':
        l0 = 1.538 if args.l0 is None else args.l0
        run_scan(("C","C"), 1, l0, args.range, args.nstep, 'C_3', 'C_3')
    elif args.case == 'C2_R':
        l0 = 1.538 if args.l0 is None else args.l0
        run_scan(("C","C"), 1, l0, args.range, args.nstep, 'C_R', 'C_R')
