#!/usr/bin/env python3
"""Boys-function approximations using generic utilities from *func_utils*.

This module demonstrates how *func_utils.match_poly_at_point* can be employed
in order to construct low-order polynomial patches that smoothly connect the
Boys function ``B(r)=1/r`` (for *r ≥ r_min*) to a regularised polynomial form
(for *r < r_min*).

Compared to *boys_approx2.py* this version is MUCH more modular:

* Symbolic matching is handled by :func:`func_utils.match_poly_at_point`.
* Plotting relies on :func:`func_utils.plot_func`.
* Adding yet another approximation now only requires specifying *powers* and
  *continuity* – no boiler-plate sympy code needed.

Run ``python boys_approx3.py`` to get a comparative plot of two example
approximations (C¹ quartic and C² sextic).
"""

from __future__ import annotations
import sys
import math
from typing import Sequence
import numpy as np
from scipy.special import erf
import matplotlib.pyplot as plt
import sympy as sp

# Local helpers ----------------------------------------------------------------
from func_utils import plot_func, match_poly_at_point  # relative import to doc/py folder package context

SQRT_PI    = math.sqrt(math.pi)
BOYS_R0    = 2 / SQRT_PI
BOYS_R0_d2 = -4/(3*SQRT_PI)  # limit r→0 of erf(r)/r


# ------------------------------------------------------------------------------
# 1) Exact Boys function (+ derivatives) on the half-line r>0
# ------------------------------------------------------------------------------

def boys_function(r: np.ndarray):
    """Return (E,F,S) for the Boys function *erf(r)/r* and derivatives.

    Parameters
    ----------
    r : (N,) ndarray
        Radii where the function shall be evaluated.  r may include the value
        0 – appropriate limits are taken.
    """
    r   = np.asarray(r)
    y   = np.empty_like(r)
    dy  = np.empty_like(r)
    dyy = np.empty_like(r)
    mask = r > 0
    # r>0 branch
    r_      = r[mask]
    erf_r   = erf(r_)
    exp_r2  = np.exp(-r_**2)
    y  [mask] = erf_r / r_
    dy [mask] = (r_ * (2/SQRT_PI * exp_r2) - erf_r) / r_**2
    #dyy[mask] = (2*erf_r / r_**3) - (4/SQRT_PI * exp_r2 / r_) - (4/SQRT_PI * exp_r2 * r_)
    dyy[mask] = (2 * erf(r_) / r_**3 ) - (4 * exp_r2 / SQRT_PI ) * (1/r_**2 + 1)
    # r=0 limit
    mask_zero = ~mask
    if mask_zero.any():
        y  [mask_zero] = BOYS_R0
        dy [mask_zero] = 0.0
        dyy[mask_zero] = BOYS_R0_d2
    return y, dy, dyy

# ------------------------------------------------------------------------------
# 2) Symbolic helper – obtain polynomial coefficients for a given configuration
# ------------------------------------------------------------------------------
def poly_coeffs(r_min: float, powers: Sequence[int], cont_order: int):
    """Return *numeric* polynomial coefficients matching 1/r at *r_min*."""
    # Build extra conditions: fix constant term c_i for power 0 to BOYS_R0
    _r      = sp.Symbol('r', positive=True)  # global sympy symbol
    _f_expr = 1/_r                           # Boys tail expression (same for all fits)
    extra = None
    if 0 in powers:
        cs   = sp.symbols(f'c0:{len(powers)}')
        idx0 = powers.index(0)
        extra = [sp.Eq(cs[idx0], BOYS_R0)]
    sol = match_poly_at_point(_f_expr, _r, r_min, powers, cont_order, extra)
    print(f"[_poly_coeffs] symbolic solution: {sol}")
    coeffs = []
    for c, v in sol.items():
        v_num = v.evalf()
        print(f"[_poly_coeffs] coeff {c} = {v_num}")
        coeffs.append(float(v_num))
    print(f"[_poly_coeffs] numeric coeffs: {coeffs}")
    return coeffs

# ------------------------------------------------------------------------------
# 3) Factory for a polynomial Boys approximation
# ------------------------------------------------------------------------------
def make_boys_poly_approx(r_min: float, powers: Sequence[int], cont_order: int):
    """Create a callable *approx(r_values) → (E,F,S)*.

    Inside *r < r_min* a polynomial with given *powers* is used.  Outside the
    exact Boys tail (1/r) is returned.
    """
    coeffs = poly_coeffs(r_min, powers, cont_order)

    # Pre-compute derivative coefficients
    p_arr = np.array(powers)
    c_arr = np.array(coeffs)
    d1_c = c_arr * p_arr
    d2_c = d1_c  * (p_arr - 1)

    def _eval(r_vals: np.ndarray):
        r_vals = np.asarray(r_vals)

        # 1. Initialize output arrays to the correct shape and type
        y   = np.zeros_like(r_vals, dtype=float)
        dy  = np.zeros_like(r_vals, dtype=float)
        ddy = np.zeros_like(r_vals, dtype=float)

        # 2. Create boolean masks for the two regions
        mask_outer = r_vals >= r_min
        mask_inner = ~mask_outer

        # 3. Handle outer region (r >= r_min) using the asymptotic 1/r form
        if np.any(mask_outer):
            r_out = r_vals[mask_outer]
            # Since r_min > 0, r_out will be positive, so no division by zero occurs.
            y[mask_outer]   =  1.0 / r_out
            dy[mask_outer]  = -1.0 / r_out**2
            ddy[mask_outer] =  2.0 / r_out**3

        # 4. Handle inner region (r < r_min) using the polynomial approximation
        if np.any(mask_inner):
            r_in = r_vals[mask_inner]
            
            rpows_y   = np.power(r_in[:, None], p_arr)
            rpows_dy  = np.power(r_in[:, None], p_arr - 1)
            rpows_ddy = np.power(r_in[:, None], p_arr - 2)
        
            y  [mask_inner] = np.dot( rpows_y   , c_arr )
            dy [mask_inner] = np.dot( rpows_dy  , d1_c  )
            ddy[mask_inner] = np.dot( rpows_ddy , d2_c  )

        return y, dy, ddy

    _eval.__name__ = f'boys_poly_deg{max(powers)}_C{cont_order}'
    return _eval

# ------------------------------------------------------------------------------
# 4) Main – quick visual comparison
# ------------------------------------------------------------------------------

if __name__ == '__main__':

    r_min = float(sys.argv[1]) if len(sys.argv) > 1 else 1.5

    cfgs = [
        # powers    C-order   color    label
        ([3,2,1,0],    2,      'g', 'C2 cubic'  ),  # this does not work
        ([4,2,0],      1,      'b', 'C1 quartic'),
        ([6,4,2,0],    2,      'r', 'C2 sextic' ),
    ]

    xs = np.linspace(0.0, 4.0, 600)

    # Plot exact Boys function first -----------------------------------------
    fig, axs = plot_func(boys_function, xs, labels=('Boys E','Boys F','Boys S'), colors=('k','k','k'))
    
    # Add approximations ------------------------------------------------------
    for cfg in cfgs:
        powers, C, color, label = cfg
        approx_fun = make_boys_poly_approx(r_min, powers, C)
        fig, (axY, axDY, axDYY) = plot_func(approx_fun, xs, axs=axs, labels=(label,None,None), colors=(color,color,color), lw=1.5, linestyle='--')
    plt.legend()
        

    #axDYY.set_ylim(-10, 10)

    # Decorate ---------------------------------------------------------------
    for ax in axs:  ax.axvline(r_min, ls=':', c='red')
    fig.suptitle(f'Boys approximations (r_min={r_min})')
    plt.tight_layout()
    plt.show()
