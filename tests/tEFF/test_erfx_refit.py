#!/usr/bin/env python3
"""Regression test for the 2026-11-08 `erfx_e6` constrained refit + `evalAE` fsi/fsj fix.

Pre-refit, `erfx_e6` (approx. of erf(kx)/x) had slope -2.6e-4 at r=0 instead of 0
=> a linear |r|-cusp in the electron-nucleus Coulomb. Paired core electrons that
relax onto nuclei then carried a persistent |F|~0.02-0.05 which no optimizer
could remove (the CH4-bPaired convergence floor). The refit enforces
y(0)=2/sqrt(pi) and dy(0)=0 exactly (e0=(sqrt(pi)/2)^(1/8), o0=-1/(8*e0^7)),
max|err| improved 3.7e-7 -> 9.1e-8. See doc/topical_audit/eFF_optimizer_speedup.md sec 4.3.

Checks (run from tests/tEFF):
  1) stock run() on CH4-bPaired converges below the old ~0.02-0.03 cusp floor
  2) FD-vs-force consistency at the relaxed state (also catches the fsi/fsj swap,
     which produced up to ~400% wrong fsize and wrong signs)
  3) the core electron sitting on the C nucleus carries ~zero force

    python3 test_erfx_refit.py
"""
import os, sys
import numpy as np

sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))
import fire_study_lib as fs
from pyBall import eFF as eff

CH4 = os.path.join(os.path.dirname(__file__), '../../cpp/common_resources/xyz/CH4.xyz')


def fd_consistency_check(h=2e-5):
    """central-FD every DOF; returns (max_abs_diff, worst_dof)."""
    x0 = eff.pDOFs.copy(); eff.eval(); f0 = eff.fDOFs.copy()
    mx, wi = 0.0, -1
    for i in range(len(x0)):
        x1 = x0.copy(); x1[i] += h;   eff.pDOFs[:] = x1; Ep = eff.eval()
        x1[i] = x0[i] - h;            eff.pDOFs[:] = x1; Em = eff.eval()
        d = abs(-(Ep - Em) / (2 * h) - f0[i])
        if d > mx: mx, wi = d, i
    eff.pDOFs[:] = x0; eff.eval()
    return mx, wi


print("=== test_erfx_refit: CH4 bPaired relaxation + FD force check ===", flush=True)
fs.set_atom_params()
eff.setVerbosity(0)

na, ne, fixed = fs.load_a(CH4, bPaired=True)
print("loaded %s: na=%d ne=%d ndof=%d" % (os.path.basename(CH4), na, ne, eff.nDOFs), flush=True)

# 1) stock run() must reach below the old ~0.02-0.03 cusp floor
eff.initOpt(dt=0.005, damping=0.1, f_limit=100.0, bMass=False)
outE, outF = eff.run(30000, dt=0.005, Fconv=1e-4, ialg=2)
Fmin = np.nanmin(outF)
print("run() |F|min = %.3g  (pre-refit floor was ~0.02-0.03)" % Fmin, flush=True)
assert Fmin < 1e-4, "FAIL: cusp floor still present (|F|min=%.3g)" % Fmin

# 2) force/energy consistency at the relaxed state (covers fsi/fsj swap + cusp)
mx, wi = fd_consistency_check()
print("FD force check: max|f - (-dE/dx)| = %.3g (dof %d)" % (mx, wi), flush=True)
# known residual: ~3e-4 esize mismatch (gauss_p8 fs vs erfx_e6 E) - see report sec 4.1
assert mx < 1e-3, "FAIL: force/energy inconsistent (max diff %.3g on dof %d)" % (mx, wi)

# 3) the coincident core electron on the C nucleus must carry ~zero force
epos = eff.epos.reshape(-1, 3); apos = eff.apos.reshape(-1, 3)
d = np.linalg.norm(epos - apos[0], axis=1)
ie = int(np.argmin(d))
fe = np.linalg.norm(eff.fDOFs[na * 3 + ie * 3: na * 3 + ie * 3 + 3])
print("e[%d] nearest nucleus dist=%.2e A  |f|=%.3g (pre-refit ~0.045)" % (ie, d[ie], fe), flush=True)
assert fe < 1e-3, "FAIL: coincident electron still carries force %.3g" % fe

print("ALL CHECKS PASSED", flush=True)
