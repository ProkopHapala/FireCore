import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from typing import Sequence, List, Any, Callable, Tuple, Dict

def numDeriv(x, y):
    """Numerical derivative using central difference"""
    dx = x[2:]-x[:-2]
    dy = y[2:]-y[:-2]
    x_ = x[1:-1]
    return -dy/dx, x_

# --- Modular Plotting Functions ---
def plot_with_deriv(ax1, ax2, x, y, y_deriv, label, color, linestyle='-'):
    ax1.plot(x, y,       label=label, color=color, linestyle=linestyle)
    ax2.plot(x, y_deriv, label=label, color=color, linestyle=linestyle)

def plot1d(x, ys, derivs=None, labels=None, colors=None, bNumDeriv=True, linestyle='-', linewidth=2, ax1=None, ax2=None, bGrid=True, bLegend=True ):
    """
    Plots a function and its derivative on given axes
    
    Args:
        ax_func: axis for function plot
        ax_deriv: axis for derivative plot
        x: x values
        y: y values (function)
        deriv: analytical derivative values (optional)
        label: legend label
        color: plot color
        show_num_deriv: whether to show numerical derivative
        linestyle: line style
        linewidth: line width
    """
    if ax1 is None:
        fig,(ax1,ax2) = plt.subplots(2,1)
    label=None
    color=None
    for i,y in enumerate(ys):
        if labels is not None: label = labels[i]
        if colors is not None: color = colors[i]
        ax1.plot(x, y, label=label, color=color, linestyle=linestyle, linewidth=linewidth)
        if bGrid: ax1.grid(True)
        if bLegend: ax1.legend()
    if derivs is not None:
        for i,dy in enumerate(derivs):
            if labels is not None: label = labels[i]
            if colors is not None: color = colors[i]
            ax2.plot(x, dy, label=f'{label} (analytical)',  color=color, linestyle=linestyle, linewidth=linewidth)
            if bNumDeriv:
                num_deriv, num_x = numDeriv(x, y)
                ax2.plot(num_x, num_deriv, label=f'{label} (numerical)', color=color, linestyle=':', linewidth=1.5, alpha=0.7)
            if bGrid: ax2.grid(True)
            if bLegend: ax2.legend()
    return fig, (ax1, ax2)

def plot1d_zip(funcs):
    ys,dys,labels = [],[],[]
    for func in funcs:
        ys.append    (func[1][0])
        dys.append   (func[1][1])
        labels.append(func[0])
    return plot1d(x, ys, derivs=dys, labels=labels)

# Example usage:
"""
Basic usage example for plot_with_deriv:

# Create figure
fig, (ax1, ax2) = plt.subplots(2, 1)

# Plot function and derivative
plot_with_deriv(ax1, ax2, x, y, dy, label='Function 1', color='blue')
plot_with_deriv(ax1, ax2, x, y2, dy2, label='Function 2', color='red', linestyle='--')

# Customize axes
ax1.set_title('Function Comparison')
ax2.set_title('Derivative Comparison')
ax1.legend()
ax2.legend()
plt.show()

Advanced usage with automatic figure creation:

# This will create a new figure automatically
fig, (ax1, ax2) = plot_with_deriv(x, y, dy, label='Function', color='green')

# Add reference lines
ax1.axhline(0, ls='--', c='k')
ax2.axhline(0, ls='--', c='k')
plt.show()
"""

# =====================================================================
# Additional generic utilities (added by Cascade on 2025-07-20)
# =====================================================================



__all__ = ['numDeriv', 'plot1d', 'plot_func', 'match_poly_at_point']

# ---------------------------------------------------------------------
# 1) plot_func : evaluate *func* on a grid and visualise E, F, S
# ---------------------------------------------------------------------

def plot_func(func, xs, params=None, axs=None, labels=('E','F','S'), colors=None, figsize=(6,9), **plot_kwargs):
    """Visualise a scalar function together with its derivatives.

    The *func* callable must support the signature::

        E, F, *rest = func(xs, **params)

    where *xs* is a 1-D `numpy.ndarray` and *E*, *F* (and optionally the second
    derivative *S*) are arrays of the same length.

    Parameters
    ----------
    func : callable
        Function returning at least *E* and *F*.
    xs : array-like
        Points at which to evaluate *func*.
    params : dict, optional
        Extra keyword arguments forwarded to *func*.
    axs : tuple(matplotlib.axes.Axes), optional
        Axes ``(axE, axF, axS)`` to plot on.  If *None* a new figure is
        created.
    labels : tuple(str), optional
        y-labels used for the three sub-plots.
    colors : tuple(str), optional
        Matplotlib colours.
    **plot_kwargs
        Further keyword arguments forwarded to ``Axes.plot``.

    Returns
    -------
    fig, axs : (matplotlib.figure.Figure, list(matplotlib.axes.Axes))
    """
    if params is None:
        params = {}
    xs = np.asarray(xs)
    out = func(xs, **params)
    if len(out) < 2:  raise ValueError('func must return at least (E, F).')
    y,dy = out[:2]
    dyy = out[2] if len(out) > 2 else None
    if axs is None:
        fig, (axY, axDY, axDYY) = plt.subplots(3, 1, sharex=True, figsize=figsize)
    else:
        axY, axDY, axDYY = axs
        fig = axY.figure
    if colors is None: colors = (None, None, None)
    axY                      .plot(xs, y,   color=colors[0], label=labels[0], **plot_kwargs); axY.legend()
    axDY                     .plot(xs, dy,  color=colors[1], label=labels[1], **plot_kwargs); axDY.legend()
    if dyy is not None: axDYY.plot(xs, dyy, color=colors[2], label=labels[2], **plot_kwargs); axDYY.legend()
    for ax, ylabel in zip((axY, axDY, axDYY), labels):
        ax.set_ylabel(ylabel)
        ax.grid(True)
    axDYY.set_xlabel('r')
    return fig, (axY, axDY, axDYY)

# ---------------------------------------------------------------------
# 2) match_poly_at_point : symbolic polynomial matching utility
# ---------------------------------------------------------------------

def match_poly_at_point(f_expr, r_sym, r0, powers, order, conds=None):
    """Construct a polynomial that matches *f_expr* up to the given derivative order.

    A polynomial ``P(r) = Σ c_i r**powers[i]`` is built and the unknown
    coefficients *c_i* are determined such that the equality

    ``d^k P/d r^k (r0) = d^k f/d r^k (r0)`` holds for ``k = 0 .. continuity_order``.

    Parameters
    ----------
    f_expr : sympy.Expr  Symbolic expression of the target function.
    r_sym  : sympy.Symbol Differentiation variable.
    r0     : float | sympy.Symbol Matching point.
    powers : Sequence[int] Exponents used in the polynomial.
    order  : int Highest derivative order to be matched (C^order continuity).
    conds  : list[sympy.Eq], optional Additional constraints (e.g. symmetry conditions).

    Returns
    -------
    dict
        Mapping *coeff_symbol → value* for the solved coefficients.
    """
    coeffs = sp.symbols(f'c0:{len(powers)}')
    poly = sum(c * r_sym**p for c, p in zip(coeffs, powers))
    eqs = [ sp.Eq(
                sp.diff( poly,   r_sym, k).subs(r_sym, r0),  
                sp.diff( f_expr, r_sym, k).subs(r_sym, r0)
            ) for k in range(order+1)  
    ]
    if conds: eqs.extend(conds)
    sol = sp.solve(eqs, coeffs, dict=True)
    if not sol: raise ValueError('Polynomial matching system has no solution.')
    return sol[0]

def solve_poly_coeffs(powers: Sequence[int], r_sym: sp.Symbol, equations: List[sp.Eq]) -> Tuple[List[float], List[int]]:
    """General symbolic solver for polynomial coefficients."""
    coeffs_sym = sp.symbols(f'c0:{len(powers)}')
    if len(equations) != len(coeffs_sym): raise ValueError(f"Number of equations ({len(equations)}) must match number of coefficients ({len(coeffs_sym)})")
    sol = sp.solve(equations, coeffs_sym)
    if not sol:raise RuntimeError(f"Symbolic solver could not find a solution for powers {powers}")
    print(f"\n--- Solving for powers {powers} ---")
    coeffs_ordered = sol
    if isinstance(sol, dict): coeffs_ordered = [sol[c] for c in coeffs_sym]        
    coeffs_num = [float(v.evalf()) for v in coeffs_ordered]
    print(f"Symbolic solution: {sol}")
    print(f"Numeric coeffs:    {coeffs_num}")
    return coeffs_num, powers

# def get_coeffs_from_recipe(r_sym: sp.Symbol, r_min_val: float, powers: Sequence[int], conditions_template: List) -> Tuple[List[float], List[int]]:
#     """Builds and solves a system of equations from a user-specified recipe."""
#     cs = sp.symbols(f'c0:{len(powers)}')
#     poly = sum(c * r_sym**p for c, p in zip(cs, powers))
#     eqs = []
#     for deriv_order, deriv_conds in enumerate(conditions_template):
#         poly_d_n = sp.diff(poly, r_sym, deriv_order)
#         for x_sym, rhs_sym in deriv_conds:
#             x_val     = x_sym   .subs(rmin, r_min_val) if hasattr(x_sym, 'subs')     else x_sym
#             rhs_temp  = rhs_sym .subs(rmin, r_min_val) if hasattr(rhs_sym, 'subs')   else rhs_sym
#             rhs_final = rhs_temp.subs(r_sym, x_val)    if hasattr(rhs_temp, 'subs')  else rhs_temp
#             eqs.append(sp.Eq(poly_d_n.subs(r_sym, x_val), rhs_final))
#     return solve_poly_coeffs(powers, r_sym, eqs)

def get_poly(x: sp.Symbol, powers: Sequence[int]):
    cs = sp.symbols(f'c0:{len(powers)}')
    return sum(c * x**p for c, p in zip(cs, powers))

def subs_conds( expr, x_sym, conditions: List ):
    """Builds and solves a system of equations from a recipe and a general substitution map."""
    eqs = []
    for x_val, rhs in conditions:
        eq = sp.Eq( expr.subs(x_sym,x_val), rhs )
        eqs.append(eq)
    return eqs

def get_polynom_approx(x_sym: sp.Symbol, powers: Sequence[int], conditions: List ):
    poly = get_poly(x_sym, powers)
    eqs = []
    for deriv_order, deriv_conds in enumerate(conditions):
        poly_d_n = sp.diff(poly, x_sym, deriv_order)
        eqs.extend( subs_conds(poly_d_n, x_sym, deriv_conds) )
    return solve_poly_coeffs(powers, x_sym, eqs)

def make_poly_approx(ref_F: Tuple[np.ndarray, ...], poly_mask: np.ndarray, coeffs: List[float], powers: List[int], label: str):
    """Create a callable approx(r) -> (E,F,S) using pre-computed data."""
    y_ref, dy_ref, ddy_ref = ref_F
    p_arr = np.array(powers)
    c_arr = np.array(coeffs)
    d1_c  = c_arr * p_arr
    d2_c  = c_arr * p_arr * (p_arr - 1)

    def _eval(x: np.ndarray):
        y, dy, ddy = np.zeros_like(x), np.zeros_like(x), np.zeros_like(x)
        outer_mask = ~poly_mask
        if np.any(outer_mask):
            y  [outer_mask] = y_ref  [outer_mask] 
            dy [outer_mask] = dy_ref [outer_mask] 
            ddy[outer_mask] = ddy_ref[outer_mask] 
        if np.any(poly_mask):
            x_in = x[poly_mask]
            with np.errstate(divide='ignore'):
                rpows_y   = np.power(x_in[:, None], p_arr)
                rpows_dy  = np.power(x_in[:, None], p_arr - 1)
                rpows_ddy = np.power(x_in[:, None], p_arr - 2)
            y  [poly_mask] = np.dot(rpows_y,   c_arr)
            dy [poly_mask] = np.nan_to_num(np.dot(rpows_dy,  d1_c))
            ddy[poly_mask] = np.nan_to_num(np.dot(rpows_ddy, d2_c))
        return y, dy, ddy
    _eval.__name__ = label
    return _eval    