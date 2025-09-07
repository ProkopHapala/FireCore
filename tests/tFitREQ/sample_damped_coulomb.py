#!/usr/bin/env python3
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf

# C/C++ Python API

sys.path.append("../../")
from pyBall import FitREQ as fit
#from pyBall.FitREQ import sample_funcEF

# -----------------------------
# Python reference implementations
# -----------------------------
SQRT_PI_INV2 = 1.1283791670955126  # 2/sqrt(pi)

def boys_exact_val_deriv(r):
    r = np.asarray(r)
    y = np.empty_like(r)
    dy = np.empty_like(r)
    small = r < 1e-12
    y[small] = SQRT_PI_INV2
    dy[small] = 0.0
    rs = r[~small]
    y[~small] = erf(rs)/rs
    dy[~small] = (SQRT_PI_INV2 * rs * np.exp(-rs*rs) - erf(rs)) / (rs*rs)
    return y, dy

# Boys polynomials (hard-coded for r_min ≈ 1.5)
B_cubic  = dict(c3=0.07607654346400593, c2=-0.3193203709421615, c0=1.12837916709551)
B_quint  = dict(c5=0.0978009599229947,  c4=-0.3186499989199512, c3=0.37187006074364537, c2=-0.3761263890318375, c0=1.12837916709551)
B_qeven  = dict(c4=0.02535884782133531, c2=-0.26226296334415705, c0=1.12837916709551)
B_sext   = dict(c6=0.010677274768021069,c4=-0.022688888634759506, c2=-0.20820925983105038, c0=1.12837916709551)

def boys_poly_val_deriv(r, mode):
    r = np.asarray(r)
    if mode == 1:  # cubic C1
        c3,c2,c0 = B_cubic['c3'], B_cubic['c2'], B_cubic['c0']
        r2 = r*r
        y   = (c3*r + c2)*r2 + c0
        dy  = (3.0*c3*r + 2.0*c2)*r
    elif mode == 2:  # quintic C2
        c5,c4,c3,c2,c0 = B_quint['c5'], B_quint['c4'], B_quint['c3'], B_quint['c2'], B_quint['c0']
        r2 = r*r
        y   = (((c5*r + c4)*r + c3)*r + c2)*r2 + c0
        dy  = (((5.0*c5*r + 4.0*c4)*r + 3.0*c3)*r + 2.0*c2)*r
    elif mode == 3:  # quartic even C1
        c4,c2,c0 = B_qeven['c4'], B_qeven['c2'], B_qeven['c0']
        r2 = r*r
        y   = (c4*r2 + c2)*r2 + c0
        dy  = (4.0*c4*r2 + 2.0*c2)*r
    else:  # 4 sextic even C2
        c6,c4,c2,c0 = B_sext['c6'], B_sext['c4'], B_sext['c2'], B_sext['c0']
        r2 = r*r
        y   = ((c6*r2 + c4)*r2 + c2)*r2 + c0
        dy  = ((6.0*c6*r2 + 4.0*c4)*r2 + 2.0*c2)*r
    return y, dy

def boys_piecewise_ref(r, rmin, mode):
    r = np.asarray(r)
    y = np.empty_like(r)
    dy = np.empty_like(r)
    mask = r < rmin
    if mode == 0:
        y[mask], dy[mask] = boys_exact_val_deriv(r[mask])
    else:
        y[mask], dy[mask] = boys_poly_val_deriv(r[mask], mode)
    # outside: 1/r
    rinv = 1.0/r[~mask]
    y[~mask]  = rinv
    dy[~mask] = -rinv*rinv
    return y, dy

# Soft/smooth clamps (positive and negative)

def soft_clamp_pos_ref(r, y1, y2):
    y0  = 1.0/r
    dy0 = -1.0/(r*r)
    y   = y0.copy()
    dy  = dy0.copy()
    mask = y0 > y1
    if np.any(mask):
        y12  = y2 - y1
        z    = (y0[mask] - y1)/y12
        inv  = 1.0/(1.0+z)
        y[mask]  = y1 + y12*(1.0 - inv)
        dy[mask] = (inv*inv) * dy0[mask]
    return y, dy

def smooth_clamp_pos_ref(r, y1, y2):
    y0  = 1.0/r
    dy0 = -1.0/(r*r)
    y   = y0.copy()
    dy  = dy0.copy()
    mask = y0 > y1
    if np.any(mask):
        y21   = y2 - y1
        z     = (y0[mask] - y1)/y21
        denom = 1.0/(1.0 + z + z*z)
        y[mask]  = y2 - y21*denom
        dy[mask] = ((1.0 + 2.0*z) * denom * denom) * dy0[mask]
    return y, dy

def soft_clamp_neg_ref(r, y1, y2):
    y0  = -1.0/r
    dy0 =  1.0/(r*r)
    y   = y0.copy()
    dy  = dy0.copy()
    mask = y0 < y1
    if np.any(mask):
        y12  = y2 - y1
        z    = (y0[mask] - y1)/y12
        inv  = 1.0/(1.0+z)
        y[mask]  = y1 + y12*(1.0 - inv)
        dy[mask] = (inv*inv) * dy0[mask]
    return y, dy

def smooth_clamp_neg_ref(r, y1, y2):
    y0  = -1.0/r
    dy0 =  1.0/(r*r)
    y   = y0.copy()
    dy  = dy0.copy()
    mask = y0 < y1
    if np.any(mask):
        y21   = y2 - y1
        z     = (y0[mask] - y1)/y21
        denom = 1.0/(1.0 + z + z*z)
        y[mask]  = y2 - y21*denom
        dy[mask] = ((1.0 + 2.0*z) * denom * denom) * dy0[mask]
    return y, dy

# Dispatcher for python reference by kind
# Boys: 10..14 with params=[rmin]
# Clamp: 20..23 with params=[y1,y2]

def python_ref(xs, kind, params):
    xs = np.asarray(xs)
    if kind in (10,11,12,13,14):
        rmin = params[0] if (params is not None and len(params)>0) else 1.5
        mode = {10:0, 11:1, 12:2, 13:3, 14:4}[kind]
        return boys_piecewise_ref(xs, rmin, mode)
    elif kind in (20,21,22,23):
        y1 = params[0] if params is not None else 3.0
        y2 = params[1] if (params is not None and len(params)>1) else 6.0
        if kind==20: return soft_clamp_pos_ref(xs, y1, y2)
        if kind==21: return smooth_clamp_pos_ref(xs, y1, y2)
        if kind==22: return soft_clamp_neg_ref(xs, y1, y2)
        if kind==23: return smooth_clamp_neg_ref(xs, y1, y2)
    elif kind==0:
        y = 1.0/xs
        dy = -1.0/(xs*xs)
        return y, dy
    else:
        raise ValueError(f"Unsupported kind={kind}")

# -----------------------------
# Diagnostics helpers
# -----------------------------

def print_every(xs, Ec, Fc, Ep, Fp, nth=100):
    n = len(xs)
    print("# i\tr\tE_cpp\tF_cpp\tE_py\tF_py\tdE\tdF")
    for i in range(0, n, max(1,nth)):
        print(f"{i}\t{xs[i]:.6f}\t{Ec[i]:.9e}\t{Fc[i]:.9e}\t{Ep[i]:.9e}\t{Fp[i]:.9e}\t{(Ec[i]-Ep[i]):.3e}\t{(Fc[i]-Fp[i]):.3e}")


def summarize(name, y, dy):
    print(f"{name}: E[min,max]=[{np.min(y):.9e},{np.max(y):.9e}]  F[min,max]=[{np.min(dy):.9e},{np.max(dy):.9e}]")

# -----------------------------
# Main
# -----------------------------

def main():
    ap = argparse.ArgumentParser(description="Sample damped Coulomb functions from C++ and compare to Python reference")
    ap.add_argument('--kind', type=int, default=11, help='0 bare; 10..14 Boys (exact,cubic,quintic,quart-even,sext-even); 20..23 clamps (soft+,smooth+,soft-,smooth-)')
    ap.add_argument('--n', type=int, default=1000, help='number of samples along r')
    ap.add_argument('--rmin', type=float, default=1e-3, help='scan start')
    ap.add_argument('--rmax', type=float, default=4.0, help='scan end')
    ap.add_argument('--boys_rmin', type=float, default=1.5, help='matching radius for Boys piecewise modes')
    ap.add_argument('--y1', type=float, default=3.0, help='soft/smooth clamp threshold y1')
    ap.add_argument('--y2', type=float, default=6.0, help='soft/smooth clamp threshold y2')
    ap.add_argument('--print-every', type=int, default=100, help='print every Nth sample')
    ap.add_argument('--save', type=str, default=None, help='save figure path')
    ap.add_argument('--show', type=int, default=1, help='show figure window')
    args = ap.parse_args()

    xs = np.linspace(args.rmin, args.rmax, args.n)

    # params for C++
    if args.kind in (10,11,12,13,14):
        params = np.array([args.boys_rmin], dtype=float)
    elif args.kind in (20,21,22,23):
        params = np.array([args.y1, args.y2], dtype=float)
    else:
        params = None

    # C++ results
    Ec, Fc = fit.sample_funcEF(xs, kind=args.kind, params=params)

    # Python reference
    Ep, Fp = python_ref(xs, args.kind, params);  Fp*=-1.0


    # Diagnostics prints
    print_every(xs, Ec, Fc, Ep, Fp, nth=args.print_every)
    summarize('C++', Ec, Fc)
    summarize('PY ', Ep, Fp)
    print(f"max|dE|={np.max(np.abs(Ec-Ep)):.3e}  max|dF|={np.max(np.abs(Fc-Fp)):.3e}")

    # Plots
    fig, (axE, axF) = plt.subplots(1,2, figsize=(8,3.0), constrained_layout=True)
    axE.plot(xs, Ec, '-r', lw=0.5, label='C++')
    axE.plot(xs, Ep, ':k', lw=1.5, label='Python')
    axE.set_title('Energy E(r)')
    axE.set_xlabel('r')
    axE.set_ylabel('E')
    axE.grid(alpha=0.3, linestyle='--')
    axE.legend(fontsize=8)

    axF.plot(xs, Fc, '-r', lw=0.5, label='C++')
    axF.plot(xs, Fp, ':k', lw=1.5, label='Python')
    axF.set_title('Force F(r) = -dE/dr')
    axF.set_xlabel('r')
    axF.set_ylabel('F')
    axF.grid(alpha=0.3, linestyle='--')
    axF.legend(fontsize=8)

    if args.save:
        print("Saving plot to:", args.save)
        fig.savefig(args.save, dpi=150, bbox_inches='tight')
    if args.show:
        plt.show()

if __name__ == '__main__':
    main()
