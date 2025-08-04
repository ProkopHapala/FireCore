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
import sympy as _sp

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

_r      = _sp.Symbol('r', positive=True)  # global sympy symbol
_f_expr = 1/_r                        # Boys tail expression (same for all fits)

def _poly_coeffs(r_min: float, powers: Sequence[int], cont_order: int):
    """Return *numeric* polynomial coefficients matching 1/r at *r_min*."""
    # Build extra conditions: fix constant term c_i for power 0 to BOYS_R0
    extra = None
    if 0 in powers:
        _cs = _sp.symbols(f'c0:{len(powers)}')
        idx0 = powers.index(0)
        extra = [_sp.Eq(_cs[idx0], BOYS_R0)]
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
    coeffs = _poly_coeffs(r_min, powers, cont_order)

    # Pre-compute derivative coefficients -------------------------------------
    # P(r) = Σ c_i r^p_i
    # P'(r) = Σ c_i p_i r^(p_i-1)
    # P''(r)= Σ c_i p_i (p_i-1) r^(p_i-2)
    p_arr = np.array(powers)
    c_arr = np.array(coeffs)

    d1_c = c_arr * p_arr
    d2_c = d1_c * (p_arr - 1)

    def _eval(r_vals: np.ndarray):
        r_vals = np.asarray(r_vals)
        E = np.where(r_vals >= r_min,  1.0 / np.where(r_vals==0, np.inf, r_vals),    0.0)
        F = np.where(r_vals >= r_min, -1.0 / np.where(r_vals==0, np.inf, r_vals)**2, 0.0)
        S = np.where(r_vals >= r_min,  2.0 / np.where(r_vals==0, np.inf, r_vals)**3, 0.0)
        mask_inner = (r_vals < r_min)
        r_in = r_vals[mask_inner]
        if r_in.size:
            E[mask_inner] = (c_arr * r_in[:,None]**p_arr).sum(axis=1)
            F[mask_inner] = (d1_c * r_in[:,None]**(p_arr-1)).sum(axis=1)
            S[mask_inner] = (d2_c * r_in[:,None]**(p_arr-2)).sum(axis=1)
        return E, F, S
    _eval.__name__ = f'boys_poly_deg{max(powers)}_C{cont_order}'
    return _eval

# ------------------------------------------------------------------------------
# 4) Main – quick visual comparison
# ------------------------------------------------------------------------------

if __name__ == '__main__':

    r_min = float(sys.argv[1]) if len(sys.argv) > 1 else 1.5

    # Two sample approximations ------------------------------------------------
    # cfgs = [
    #     dict(powers=[4,2,0],   C=1, color='tab:blue',  label='C1 quartic'),
    #     dict(powers=[6,4,2,0], C=2, color='tab:orange',label='C2 sextic'),
    # ]

    cfgs = [
        # powers    C-order   color    label
        ([4,2,0],      1,      'b', 'C1 quartic'),
        ([6,4,2,0],    2,      'r', 'C2 sextic' ),
    ]

    xs = np.linspace(0.0, 4.0, 600)

    # Plot exact Boys function first -----------------------------------------
    fig, axs = plot_func(boys_function, xs, labels=('Boys E','Boys F','Boys S'), colors=('k','k','k'))

    # for cfg in cfgs:
    #     approx_fun = make_boys_poly_approx(r_min, cfg['powers'], cfg['C'])
    #     fig, (axY, axDY, axDYY) = plot_func(approx_fun, xs, axs=axs, labels=(cfg['label'],None,None), colors=(cfg['color'],cfg['color'],cfg['color']), lw=1.5, linestyle='--')

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
