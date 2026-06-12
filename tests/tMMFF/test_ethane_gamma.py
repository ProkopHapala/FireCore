#!/usr/bin/env python3
import os, sys, argparse
import numpy as np

sys.path.append("../../")
from pyBall import MMFF

xyz = os.path.join(os.path.dirname(__file__), "..", "..", "cpp", "common_resources", "xyz", "C2H6.xyz")

with open(xyz) as f:
    natoms_xyz = int(f.readline().strip())

def load_xyz_symbols(fname, natoms):
    syms=[]
    with open(fname) as f:
        nat = int(f.readline().strip()); assert nat==natoms
        f.readline()
        for _ in range(natoms):
            syms.append(f.readline().split()[0])
    return syms

_MASS_AMU = {"H":1.00784, "C":12.0107, "N":14.0067, "O":15.999}

def masses_amu_from_symbols(syms):
    m=np.zeros(len(syms),dtype=float)
    for i,s in enumerate(syms):
        if s not in _MASS_AMU: raise ValueError(f"Unknown symbol '{s}' (add to _MASS_AMU)")
        m[i]=_MASS_AMU[s]
    return m

def thz_from_eig_massweighted(eig_eVA2_per_amu):
    eV=1.602176634e-19
    Ang=1e-10
    amu=1.66053906660e-27
    conv = (eV/(Ang*Ang))/amu
    w = np.sign(eig_eVA2_per_amu)*np.sqrt(np.abs(eig_eVA2_per_amu)*conv)
    f = w/(2*np.pi)*1e-12
    return f

def run_case(bUFF, bBondsOnly=False):
    print("\n==============================")
    print(f"Ethane gamma Hessian bUFF={int(bUFF)} bBondsOnly={int(bBondsOnly)}")
    print("==============================")
    MMFF.init(xyz_name=xyz, nPBC=(0,0,0), bEpairs=False, bMMFF=True, bUFF=bUFF)
    if bBondsOnly:
        MMFF.setSwitches(PBC=-1, NonBonded=-1, Angles=-1, PiSigma=-1, PiPiI=-1)
    else:
        MMFF.setSwitches(NonBonded=-1, PBC=-1)
    inds = np.arange(natoms_xyz, dtype=np.int32)
    H = MMFF.getHessian3Nx3N(inds, dx=1e-4)
    H = 0.5*(H+H.T)

    lam = np.linalg.eigvalsh(H)
    lam = np.sort(lam)
    print("stiffness eig λ  [eV/Å^2]  min..max:", lam[0], lam[-1])
    print("sqrt(λ)           [sqrt(eV)/Å] last 12:", np.sign(lam[-12:])*np.sqrt(np.abs(lam[-12:])))

    syms = load_xyz_symbols(xyz, natoms_xyz)
    m = masses_amu_from_symbols(syms)
    w = np.repeat(1.0/np.sqrt(m), 3)
    Hw = (H*w[None,:])*w[:,None]
    lamw = np.linalg.eigvalsh(0.5*(Hw+Hw.T))
    lamw = np.sort(lamw)
    fTHz = thz_from_eig_massweighted(lamw)
    print("mass-weighted eig μ [eV/Å^2/amu] min..max:", lamw[0], lamw[-1])
    print("freq               [THz]         last 12:", fTHz[-12:])

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--bondsOnly', action='store_true')
    args = parser.parse_args()
    run_case(False, bBondsOnly=args.bondsOnly)
    run_case(True,  bBondsOnly=args.bondsOnly)
