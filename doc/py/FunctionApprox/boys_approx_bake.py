from __future__ import annotations
import sys
import math
import numpy as np
from scipy.special import erf
import matplotlib.pyplot as plt
import sympy as sp

# ==============================================================================
# NOTE: The func_utils.py file was not provided, so I have recreated the
# essential functions here based on their usage in the main script.
# You can replace this section with your actual func_utils.py import.
# ==============================================================================

def get_polynom_approx(r_sym, powers, conditions):
    """
    Solves for polynomial coefficients that satisfy the given Hermite interpolation conditions.
    """
    print(f"--- Solving for powers {powers} ---")
    coeffs_sym = [sp.Symbol(f'c{i}') for i in range(len(powers))]
    
    # Build the polynomial P(r)
    poly = sum(c * r_sym**p for c, p in zip(coeffs_sym, powers))
    
    # Build system of linear equations from conditions
    equations = []
    for deriv_order, points in enumerate(conditions):
        # Differentiate the polynomial if needed
        p_deriv = sp.diff(poly, r_sym, deriv_order)
        for r_val, target_val in points:
            eq = p_deriv.subs(r_sym, r_val) - target_val
            equations.append(eq)
            
    # Solve the system
    solution = sp.solve(equations, coeffs_sym)
    print(f"Symbolic solution: {solution}")
    
    # Return numeric coefficients in the order of the provided powers
    numeric_coeffs = [float(solution[c]) for c in coeffs_sym]
    print(f"Numeric coeffs:    {numeric_coeffs}")
    
    return numeric_coeffs, poly.subs(solution)

def make_poly_approx(ref_func_tuple, mask, coeffs, powers, label):
    """
    Creates a callable function that returns the polynomial approximation
    where the mask is True, and the reference function otherwise.
    coeffs: list of coefficients c_i corresponding to powers[i]
    powers: list of integer powers p_i (may be non-consecutive)
    """
    ref_y, ref_dy, ref_dyy = ref_func_tuple

    coeffs = np.asarray(coeffs, dtype=float)
    powers = np.asarray(powers, dtype=int)

    def approx_func(r_in):
        r_in = np.asarray(r_in)
        # Start with the reference function values
        y_out, dy_out, dyy_out = ref_y.copy(), ref_dy.copy(), ref_dyy.copy()

        # Evaluate the polynomial and its derivatives on the masked region
        r_masked = r_in[mask]
        if r_masked.size:
            # y = sum c_i * r^p
            y_val = np.zeros_like(r_masked, dtype=float)
            dy_val = np.zeros_like(r_masked, dtype=float)
            dyy_val = np.zeros_like(r_masked, dtype=float)
            for c, p in zip(coeffs, powers):
                if c == 0.0:
                    continue
                rp = r_masked**p if p != 0 else np.ones_like(r_masked)
                y_val += c * rp
                if p >= 1:
                    dy_val += c * p * (r_masked**(p-1))
                if p >= 2:
                    dyy_val += c * p * (p-1) * (r_masked**(p-2))

            y_out[mask] = y_val
            dy_out[mask] = dy_val
            dyy_out[mask] = dyy_val

        return y_out, dy_out, dyy_out

    return approx_func

def plot_func(func, xs, axs=None, labels=(None, None, None), **kwargs):
    """
    Plots a function and its derivatives. Always plots all three series.
    If a label element is None, the line will be plotted without a legend entry.
    """
    if axs is None:
        fig, axs = plt.subplots(3, 1, sharex=True, figsize=(10, 15))
    else:
        fig = axs[0].get_figure()

    ys, dys, dyys = func(xs)

    # Plot on all subplots; include label only if provided
    kw0 = dict(kwargs)
    if labels[0] is not None:
        axs[0].plot(xs, ys, label=labels[0], **kw0)
    else:
        axs[0].plot(xs, ys, **kw0)

    kw1 = dict(kwargs)
    if labels[1] is not None:
        axs[1].plot(xs, dys, label=labels[1], **kw1)
    else:
        axs[1].plot(xs, dys, **kw1)

    kw2 = dict(kwargs)
    if labels[2] is not None:
        axs[2].plot(xs, dyys, label=labels[2], **kw2)
    else:
        axs[2].plot(xs, dyys, **kw2)

    axs[0].set_ylabel('Energy (E)')
    axs[1].set_ylabel('Force (F) = -E\'')
    axs[2].set_ylabel('Stiffness (S) = -F\' = E\'\'')
    axs[2].set_xlabel('r')

    return fig, axs

# ==============================================================================
# End of recreated func_utils
# ==============================================================================


# Constants
SQRT_PI    = math.sqrt(math.pi)
BOYS_R0    = 2 / SQRT_PI
BOYS_R0_d1 = 0.0
BOYS_R0_d2 = -4/(3*SQRT_PI)

# ------------------------------------------------------------------------------
# 1) Exact and Reference functions
# ------------------------------------------------------------------------------
def boys_function(r: np.ndarray):
    """Return (E,F,S) for the Boys function *erf(r)/r* and derivatives."""
    r          = np.asarray(r)
    y, dy, dyy = np.empty_like(r), np.empty_like(r), np.empty_like(r)
    mask       = r > 1e-9 # Use a small epsilon to avoid division by zero
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
    """Return (E,F,S) for the Coulomb function *1/r* and derivatives."""
    r_safe = np.asarray(r).copy()
    mask_zero = r_safe < 1e-9
    r_safe[mask_zero] = 1e-9 # Avoid division by zero for plotting
    y   = 1 / r_safe
    dy  = -1 / r_safe**2
    dyy = 2 / r_safe**3
    return y, dy, dyy

# ------------------------------------------------------------------------------
# 2) Explicit hardcoded polynomial approximations for r_min = 1.5
# ------------------------------------------------------------------------------

def approx_c1_hermite_cubic(r: np.ndarray):
    c3, c2, c0 = 0.07607654346400593, -0.3193203709421615, 1.12837916709551
    r2 = r * r
    y   = (c3 * r + c2) * r2 + c0
    dy  = (3.0 * c3 * r + 2.0 * c2) * r
    dyy = 6.0 * c3 * r + 2.0 * c2
    return y, dy, dyy

def approx_c2_hermite_quintic(r: np.ndarray):
    c5, c4, c3, c2, c0 = 0.0978009599229947, -0.3186499989199512, 0.37187006074364537, -0.3761263890318375, 1.12837916709551
    r2 = r * r
    y   = (((c5 * r + c4) * r + c3) * r + c2) * r2 + c0
    dy  = (((5.0 * c5 * r + 4.0 * c4) * r + 3.0 * c3) * r + 2.0 * c2) * r
    dyy = ((20.0 * c5 * r + 12.0 * c4) * r + 6.0 * c3) * r + 2.0 * c2
    return y, dy, dyy

def approx_c1_quartic_even(r: np.ndarray):
    c4, c2, c0 = 0.02535884782133531, -0.26226296334415705, 1.12837916709551
    r2 = r * r
    y   = (c4 * r2 + c2) * r2 + c0
    dy  = (4.0 * c4 * r2 + 2.0 * c2) * r
    dyy = 12.0 * c4 * r2 + 2.0 * c2
    return y, dy, dyy

def approx_c2_sextic_even(r: np.ndarray):
    c6, c4, c2, c0 = 0.010677274768021069, -0.022688888634759506, -0.20820925983105038, 1.12837916709551
    r2 = r * r
    y   = ((c6 * r2 + c4) * r2 + c2) * r2 + c0
    dy  = ((6.0 * c6 * r2 + 4.0 * c4) * r2 + 2.0 * c2) * r
    dyy = (30.0 * c6 * r2 + 12.0 * c4) * r2 + 2.0 * c2
    return y, dy, dyy

def combine_approx_with_coulomb(approx_func, r_min):
    """Wrapper to combine a polynomial approx (for r<r_min) with Coulomb (for r>=r_min)."""
    def combined_func(r):
        r = np.asarray(r)
        y, dy, dyy = np.empty_like(r), np.empty_like(r), np.empty_like(r)
        
        mask_approx = r < r_min
        mask_coulomb = ~mask_approx
        
        # Approximation part
        y_a, dy_a, dyy_a = approx_func(r[mask_approx])
        y[mask_approx], dy[mask_approx], dyy[mask_approx] = y_a, dy_a, dyy_a
        
        # Coulomb part
        y_c, dy_c, dyy_c = coulomb_function(r[mask_coulomb])
        y[mask_coulomb], dy[mask_coulomb], dyy[mask_coulomb] = y_c, dy_c, dyy_c
        
        return y, dy, dyy
    return combined_func

# ------------------------------------------------------------------------------
# 3) Main – Data-driven configurations and plotting
# ------------------------------------------------------------------------------
if __name__ == '__main__':
    r_min = float(sys.argv[1]) if len(sys.argv) > 1 else 1.5

    r_sym = sp.Symbol('r', positive=True)
    
    y1   = 1/r_min
    dy1  = -1/r_min**2
    ddy1 = 2/r_min**3
    
    C2_HERMITE_CONDITIONS = [
        [(0, BOYS_R0)   , (r_min,   y1)], # P(0),   P(rmin)
        [(0, BOYS_R0_d1), (r_min,  dy1)], # P'(0),  P'(rmin)
        [(0, BOYS_R0_d2), (r_min, ddy1)], # P''(0), P''(rmin)
    ]
    C2_EVEN_CONDITIONS = [
        [(r_min,   y1),(0, BOYS_R0)], #  P(rmin), P(0),
        [(r_min,  dy1)],             # P'(rmin)
        [(r_min, ddy1)],             # P''(rmin)
    ]
    # C1 conditions (for cubic Hermite and quartic even)
    C1_HERMITE_CONDITIONS = [
        [(0, BOYS_R0)   , (r_min,   y1)], # P(0),  P(rmin)
        [(0, BOYS_R0_d1), (r_min,  dy1)], # P'(0), P'(rmin)
    ]
    C1_EVEN_CONDITIONS = [
        [(r_min,   y1),(0, BOYS_R0)], #  P(rmin), P(0),
        [(r_min,  dy1)],              # P'(rmin)
    ]
    
    # Sympy-based configurations
    cfgs = [
        ('C1 Hermite Cubic',   'g', [3,2,1,0],     C1_HERMITE_CONDITIONS),
        ('C2 Hermite Quintic', 'm', [5,4,3,2,1,0], C2_HERMITE_CONDITIONS),
        ('C1 Quartic (Even)',  'c', [4,2,0],       C1_EVEN_CONDITIONS),
        ('C2 Sextic (Even)',   'r', [6,4,2,0],     C2_EVEN_CONDITIONS),
    ]

    # Explicit function configurations for verification
    explicit_cfgs = [
        ('C1 Hermite Cubic (Explicit)',   'g', approx_c1_hermite_cubic),
        ('C2 Hermite Quintic (Explicit)', 'm', approx_c2_hermite_quintic),
        ('C1 Quartic (Even) (Explicit)',  'c', approx_c1_quartic_even),
        ('C2 Sextic (Even) (Explicit)',   'r', approx_c2_sextic_even),
    ]

    xs = np.linspace(0.0, 4.0, 600)
    fig, axs = plot_func(boys_function, xs, labels=('Boys', 'Boys F', 'Boys S'), lw=1.5, color='k')
    plot_func(coulomb_function, xs, labels=('Coulomb', 'Coulomb F', 'Coulomb S'), lw=1.5, ls=':', color='gray', axs=axs)

    ref_F = coulomb_function(xs)
    poly_mask = xs < r_min
    
    # Plot Sympy-generated functions
    print("--- Generating approximations with Sympy ---")
    for cfg in cfgs:
        label, color, powers, conditions = cfg
        coeffs, _ = get_polynom_approx(r_sym, powers, conditions)
        approx_fun = make_poly_approx(ref_F, poly_mask, coeffs, powers, label)
        plot_func(approx_fun, xs, axs=axs, labels=(label, None, None), color=color, lw=0.5, linestyle='-')
        
    # Plot explicit hardcoded functions for verification
    print("\n--- Plotting explicit functions for verification ---")
    for label, color, func in explicit_cfgs:
        combined_func = combine_approx_with_coulomb(func, r_min)
        plot_func(combined_func, xs, axs=axs, labels=(label, None, None), color=color, lw=1.5, linestyle=':')
        
    for ax in axs:
        ax.axvline(r_min, ls=':', c='gray', label=f'r_min={r_min}' if ax is axs[0] else "")
    
    axs[0].legend()
    axs[0].set_ylim(-0.1, 1.5)
    axs[1].set_ylim(-1.0, 0.2)
    axs[2].set_ylim(-1.0, 1.0)
    for ax in axs: ax.grid(True, linestyle='--', alpha=0.6)
    fig.suptitle(f'Boys Function Approximations (r_min={r_min})')
    plt.tight_layout()
    plt.show()