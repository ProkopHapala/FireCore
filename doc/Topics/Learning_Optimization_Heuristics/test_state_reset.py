#!/usr/bin/env python3
"""§11.2B A->B->A contamination test: a mode-0 run must not inherit adaptive invMasses.
Sequence in ONE worker: A(autoMass=2, im=3) -> B(autoMass=0, im=3) vs B_fresh(autoMass=0 first).
Before the fix B != B_fresh (inherited adaptive D). After: B == B_fresh exactly (deterministic)."""
import os, sys, numpy as np
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
os.chdir(os.path.join(_HERE, '..', '..', '..'))   # repo root: force-field resource paths
import scan_all_ff as S

NMAX, FCONV = 10000, 1e-6

def gh_evals(mmff, am, label, p0, DOFs):
    DOFs[:] = p0
    mmff.setAutoMass(am, 100.0)
    oE = np.zeros(NMAX+1); oF = np.zeros(NMAX+1)
    n = mmff.run_GH(nstepMax=NMAX, Fconv=FCONV, iMethod=3, mHist=8, eta=0.02, a_init=-1, outE=oE, outF=oF)
    F = np.nanmin(oF[:n])
    print(f"{label:12s} evals={n} minF={F:.3e}", flush=True)
    return n, F

def main():
    mmff, p0, E0 = S.init_mmff('propandiol', False)
    DOFs = mmff.DOFs
    nB, fB   = gh_evals(mmff, 0, "B_fresh", p0, DOFs)   # autoMass=0 first: unit masses
    nA, _    = gh_evals(mmff, 2, "A_autoM2", p0, DOFs)   # learns adaptive masses
    nB2, fB2 = gh_evals(mmff, 0, "B_after_A", p0, DOFs)  # must equal B_fresh if masses restored
    ok = (nB == nB2) and (fB == fB2)
    print(("PASS" if ok else "FAIL") + f" A->B->A reset: B_fresh={nB} B_after_A={nB2}", flush=True)
    sys.exit(0 if ok else 1)

if __name__ == '__main__': main()
