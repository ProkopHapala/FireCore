from __future__ import annotations
import sys
import math
import numpy as np
from scipy.special import erf
import matplotlib.pyplot as plt
import sympy as sp

# Local helpers
from func_utils import plot_func, get_polynom_approx, make_poly_approx


# Constants
SQRT_PI    = math.sqrt(math.pi)
BOYS_R0    = 2 / SQRT_PI
BOYS_R0_d1 = 0.0
BOYS_R0_d2 = -4/(3*SQRT_PI)

# ------------------------------------------------------------------------------
# 1) Exact Boys function
# ------------------------------------------------------------------------------
def boys_function(r: np.ndarray):
    """Return (E,F,S) for the Boys function *erf(r)/r* and derivatives."""
    r          = np.asarray(r)
    y, dy, dyy = np.empty_like(r), np.empty_like(r), np.empty_like(r) # THE ONLY FIX IS HERE
    mask       = r > 0
    r_         = r[mask]
    erf_r      = erf(r_)
    exp_r2     = np.exp(-r_**2)
    y  [mask]  = erf_r / r_
    dy [mask]  = (r_ * (2/SQRT_PI * exp_r2) - erf_r) / r_**2
    dyy[mask]  = (2 * erf_r / r_**3) - (4 * exp_r2 / SQRT_PI) * (1/r_**2 + 1)
    mask_zero  = ~mask
    if mask_zero.any():
        y  [mask_zero] = BOYS_R0
        dy [mask_zero] = BOYS_R0_d1
        dyy[mask_zero] = BOYS_R0_d2
    return y, dy, dyy

def coulomb_function(r: np.ndarray):
    """Return (E,F,S) for the Boys function *erf(r)/r* and derivatives."""
    r  = np.asarray(r)
    y  =  1/r
    dy = -1/r**2
    dyy=  2/r**3
    return y, dy, dyy

# ------------------------------------------------------------------------------
# 5) Main – Data-driven configurations
# ------------------------------------------------------------------------------
if __name__ == '__main__':
    r_min = float(sys.argv[1]) if len(sys.argv) > 1 else 1.5

    r_sym   = sp.Symbol('r', positive=True)
    # y_sym   = 1/r_sym 
    # dy_sym  = sp.diff( y_sym , r_sym) 
    # dyy_sym = sp.diff(dy_sym , r_sym)

    y1    =  1/r_min
    dy1   = -1/r_min**2
    ddy1  =  2/r_min**3
    

    C1_HERMITE_CONDITIONS = [
        [(0, BOYS_R0    ) , (r_min,  y1)], # P(0),  P(rmin)
        [(0, BOYS_R0_d1 ) , (r_min, dy1)], # P'(0), P'(rmin)
    ]
    C2_HERMITE_CONDITIONS = [
        [(0, BOYS_R0)   , (r_min,   y1)], # P(0),   P(rmin)
        [(0, BOYS_R0_d1), (r_min,  dy1)], # P'(0),  P'(rmin)
        [(0, BOYS_R0_d2), (r_min, ddy1)], # P''(0), P''(rmin)
    ]
    C2_EVEN_CONDITIONS = [
        [(r_min,   y1),(0, BOYS_R0)], #  P(rmin), P(0),
        [(r_min,  dy1)], # P' (rmin)
        [(r_min, ddy1)], # P''(rmin)
    ]
    C1_EVEN_CONDITIONS = [
        [(r_min,   y1),(0, BOYS_R0)], #  P(rmin), P(0),
        [(r_min,  dy1)], # P' (rmin)
    ]
    cfgs = [
        #   label,              color, powers,        conditions
        ( 'C1 Hermite Cubic',   'g', [3,2,1,0],     C1_HERMITE_CONDITIONS),
        ( 'C2 Hermite Quintic', 'm', [5,4,3,2,1,0], C2_HERMITE_CONDITIONS),
        ( 'C1 Quartic (Even)',  'c', [4,2,0  ],     C1_EVEN_CONDITIONS),
        ( 'C2 Sextic (Even)',   'r', [6,4,2,0],     C2_EVEN_CONDITIONS),
    ]

    xs = np.linspace(0.0, 4.0, 600)
    fig, axs = plot_func(boys_function,     xs, labels=('Boys E','Boys F','Boys S'),          lw=1.0,         colors=('k','k','k'), figsize=(10,15))
    plot_func            (coulomb_function, xs, labels=('Coulomb E','Coulomb F','Coulomb S'), lw=1.5, ls=':', colors=('gray','gray','gray'), axs=axs)

    ref_F     = coulomb_function(xs)
    poly_mask = xs < r_min
    for cfg in cfgs:
        label, color, powers, conditions = cfg
        coeffs, _  = get_polynom_approx(r_sym, powers, conditions)
        print(f"get_polynom_approx): ", coeffs )
        approx_fun = make_poly_approx(ref_F, poly_mask, coeffs, powers, label)
        plot_func( approx_fun, xs, axs=axs, labels=(label,None,None), colors=(color,color,color), lw=1.0, linestyle='-')
        
    for ax in axs:
        ax.axvline(r_min, ls=':', c='gray', label=f'r_min={r_min}' if ax is axs[0] else "")
    axs[0].legend()  
    axs[0].set_ylim(-0.1, 1.5)
    axs[1].set_ylim(-1, 0.2)
    axs[2].set_ylim(-1, 1)
    axs[0].axhline(0, ls='--', c='k')
    axs[1].axhline(0, ls='--', c='k')
    axs[2].axhline(0, ls='--', c='k')
    fig.suptitle(f'Boys Function Approximations (r_min={r_min})')
    plt.tight_layout()
    plt.show()