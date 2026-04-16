"""
Thin interface to STM.compute_dos — runs Fireball SCF + spectral function.
Run as subprocess to avoid Fireball reinit issues.

Usage:
  python stm_compute_dos.py --mol H2   --out export/stm_phase1/dos_H2.npz
  python stm_compute_dos.py --xyz PTCDA.xyz --out export/stm_phase1/dos_PTCDA.npz
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
