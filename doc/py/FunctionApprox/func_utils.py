import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

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

def plot_func(func, xs, params=None, axs=None, labels=('E','F','S'), colors=None, **plot_kwargs):
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
        fig, (axY, axDY, axDYY) = plt.subplots(3, 1, sharex=True, figsize=(6, 9))
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