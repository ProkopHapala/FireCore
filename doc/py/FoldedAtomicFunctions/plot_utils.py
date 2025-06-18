# --- 5. Plotting & Analysis Utilities ---

import numpy as np
import matplotlib.pyplot as plt

def plot1D(
    xs: np.ndarray,
    ys: np.ndarray, 
    title: str,
    ylims: tuple = None,
    scMin: float = 1.5,
    bLogY: bool = False,
    ls='-',lw=0.5,
    ax=None,
    labels=None,
):
    """Plots 1D profiles."""
    n = len(ys)
    if ax is None: 
        fig, ax = plt.subplots(figsize=(12, 7))
    ymin=1e+300
    for i in range(n):
        label = labels[i] if labels is not None else f'{i+1}'
        ax.plot(xs, ys[i], ls, label=label, lw=lw, alpha=1.0)
        ymini = np.min(ys[i])
        if ymini < ymin: ymin = ymini

    if bLogY:
        ax.set_yscale('log')
    elif ylims is not None:
        ax.set_ylim(ylims)
    elif scMin is not None:
        vmin=ymin * scMin
        ax.set_ylim(vmin, -vmin)
    ax.set_title(title)
    ax.set_xlabel('z (Å)')
    ax.set_ylabel('Potential / Value')
    ax.legend() # Legend for main plot
    ax.grid(True)
    #if filename: plt.savefig(filename)
    #plt.show()
    return ax

# def makeTwin(ax, xs, ys, label=None, ls=':', lw=0.5, c='k'): # Keep existing function
#     ax2 = ax.twinx()
#     ax2.plot(xs, ys, ls, c=c, label=label, lw=lw)
#     ax2.set_ylabel(label)
#     ax2.tick_params(axis='y')
#     return ax2

def plot_SV(s_vals, K_opt):
    """Plots singular values."""
    if s_vals.size == 0:
        print("No singular values to plot.")
        return
    plt.figure(figsize=(8, 5))
    plt.plot(range(1, len(s_vals) + 1), s_vals, 'o-')
    plt.title('Singular Values from SVD of Sample Coefficients Matrix (S)')
    plt.xlabel('Component Number')
    plt.ylabel('Singular Value')
    if K_opt > 0 and K_opt <= len(s_vals):
      plt.axvline(K_opt, color='r', linestyle='--', label=f'Selected K={K_opt}')
      plt.legend()
    plt.grid(True)
    plt.yscale('log')
    # if filename: plt.savefig(filename)
    # plt.show()

def plotMultiFunctionApprox(xs, data_pairs, bError=False, colors=None, errMax=None, scMin=None, title='Function Approximation'):
    """
    Plots multiple reference functions and their approximations.

    Parameters
    ----------
    xs : np.ndarray
        The x-coordinates (z-values).
    data_pairs : list of tuples (y_ref, y_approx)
        List where each tuple contains a reference function and its approximation.
    bError : bool, default: False
        Whether to plot the error (difference) on a secondary y-axis.
    colors : list of colors, optional
        List of colors to cycle through for each pair of ref/approx.
    errMax : float, optional
        Maximum absolute value for the error y-axis limits.
    scMin : float, optional
        Scaling factor for the main y-axis limits based on the minimum reference value.
    title : str, default: 'Function Approximation'
        Title for the main plot.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object.
    axes : tuple of matplotlib.axes.Axes
        A tuple containing the main axis and the error axis (if bError is True).
    """
    fig, ax1 = plt.subplots(figsize=(12, 7))
    ax2 = None
    if bError:
        ax2 = ax1.twinx()
        ax2.set_ylabel('Error', color='r')
        ax2.tick_params(axis='y', labelcolor='r')
        ax2.axhline(0, color='gray', linestyle=':', linewidth=0.5)
    n_pairs = len(data_pairs)
    if colors is None: colors = plt.cm.get_cmap('tab10', n_pairs)
    for i, (y_ref, y_app) in enumerate(data_pairs):
        color = colors(i) if callable(colors) else colors[i % len(colors)]
        ax1.plot(xs, y_ref, label=f'Sample {i+1}', ls=':', lw=1.5, c=color)
        ax1.plot(xs, y_app, label=f'Approx {i+1}', ls='-', lw=1.0, c=color)
        if bError: ax2.plot(xs, y_ref - y_app, color=color, linestyle='--', lw=0.8, label=f'Err {i+1}')
    ax1.set_xlabel('z (Å)'); ax1.set_ylabel('Value'); ax1.set_title(title); ax1.grid(True, linestyle=':', alpha=0.7)
    if scMin is not None: vmin = np.min(np.concatenate([pair[0] for pair in data_pairs])); ax1.set_ylim(vmin * scMin, -vmin * scMin)
    if bError and errMax is not None: ax2.set_ylim(-errMax, errMax)
    lines1, labels1 = ax1.get_legend_handles_labels(); lines2, labels2 = ax2.get_legend_handles_labels() if bError else ([], [])
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='best'); fig.tight_layout()
    return fig, (ax1, ax2) if bError else (ax1,)

def plotFunctionApprox( xs, y_ref, ys_approx, bError=False, colors=None, errMax=1.0e-3, scMin=1.5 ):
    fig, ax1 = plt.subplots(figsize=(12, 7))
    n_approx = len(ys_approx)
    if colors is None:
        colors = plt.cm.jet(np.linspace(0, 1, n_approx))
    # Setup secondary axis for errors
    if bError:
        ax2 = ax1.twinx()
        #ax2.set_ylabel(f'Error * {scErr:.0f}', color='gray')
        ax2.set_ylabel(f'Error *', color='r')
        ax2.tick_params(axis='y', labelcolor='gray')
        ax2.axhline(0, color='gray', linestyle=':', linewidth=0.5)
    for i in range(n_approx):
        #poly_approx, z_cut_n = eval_poly_at_z0(xs, k0, n_val, z0_approx)
        y_app, z_cut_n, label = ys_approx[i]
        ax1.plot(xs, y_app, label=label+f' z_cut={z_cut_n:.3f}', ls='-', lw=1.0, c=colors[i])
        # Mark z_cut if it's within the plot range
        ax1.axvline(z_cut_n, c=colors[i], ls=':', lw=0.8, ymax=0.7)
        # Plot error for this approximation
        if bError:
            error = (y_ref - y_app) #* scErr
            ax2.plot(xs, error, color=colors[i], linestyle='-.', lw=0.8, label=f'Err {label}')
            if bError: ax2.set_ylim(-errMax, errMax)
    ax1.axhline(0, color='k', linestyle='--', linewidth=0.5)
    ax1.plot(xs, y_ref,':k', label=f'Target', lw=2.0)
    ax1.set_xlabel('z')
    ax1.set_ylabel('Value')
    if scMin is not None:
        vmin = np.min(y_ref)
        ax1.set_ylim(vmin * scMin, -vmin * scMin)
    ax1.set_xlim(xs[0], xs[-1])
    ax1.grid(True, linestyle=':', alpha=0.7)
    lines, labels = ax1.get_legend_handles_labels()
    if bError:
        lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, loc='center right')
    fig.tight_layout()
    #plt.savefig("exp_kz_approximation_sequence.png")
    #plt.show()
    return fig,(ax1,ax2) 


# ===================================================================== # Keep existing function
# 2-D plotting – compact helpers (imshow based)
# =====================================================================

def imshow_grid(grid, extent, title="", atoms=None, cmap="RdBu_r"):
    """Quick imshow of a 2-D grid with optional atom markers.

    Parameters
    ----------
    grid    : 2-D ndarray (shape (Nx, Nz)).
    extent  : [xmin, xmax, zmin, zmax] for imshow.
    atoms   : iterable of dicts with keys `x`, `z`, optional `r0` (size) & `color`.
    """
    fig, ax = plt.subplots(figsize=(10, 4))
    im = ax.imshow(grid.T, origin="lower", aspect="auto", extent=extent, cmap=cmap)
    plt.colorbar(im, ax=ax, shrink=0.8)

    if atoms is not None:
        for a in atoms:
            sz = (a.get("r0", 1.0) * 7) ** 2
            ax.scatter(a["x"], a["z"], s=sz, c=a.get("color", "k"), edgecolors="w", linewidths=0.5)
    ax.set_xlabel("x (Å)"); ax.set_ylabel("z (Å)")
    ax.set_title(title)
    #if fname: plt.savefig(fname); print("saved", fname)
    #plt.tight_layout(); plt.show(); 
    return ax

def plot2Dapprox(ref, fit, extent, title="Fit vs Ref", cmap="RdBu_r"):
    """Show reference, fit and error in one row of imshows."""
    err = fit - ref
    data = [ref, fit, err]
    lbls = ["Reference", "Fit", "Error"]
    vmin = min(ref.min(), fit.min())
    vmax = max(ref.max(), fit.max())
    fig, axs = plt.subplots(1, 3, figsize=(14, 4))
    for ax, d, l in zip(axs, data, lbls):
        im = ax.imshow(d.T, origin="lower", aspect="auto", extent=extent, cmap=cmap, vmin=vmin, vmax=vmax)
        ax.set_title(l)
        ax.set_axis_off()
    plt.colorbar(im, ax=axs.ravel().tolist(), shrink=0.7)
    fig.suptitle(title)
    #if fname:   plt.savefig(fname); print("saved", fname)
    #plt.tight_layout(); plt.show(); 
    return axs

def plot2Dbasis(rows, shape, extent, coeffs=None, labels=None, ncol=4, cmap="RdBu_r"):
    """Grid of imshows for basis rows (each row flattened)."""
    K = rows.shape[0]
    ncol = min(ncol, K)
    nrow = (K + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(3*ncol, 3*nrow))
    axes = np.atleast_2d(axes)
    for k in range(K):
        ax = axes[k//ncol, k % ncol]
        im = ax.imshow(rows[k].reshape(shape).T, origin="lower", aspect="auto", extent=extent, cmap=cmap)
        t = labels[k] if labels else f"ϕ{k}"
        if coeffs is not None and k < len(coeffs):
            t += f"\nC={coeffs[k]:.1e}"
        ax.set_title(t, fontsize=7); ax.set_axis_off()
    for ax in axes.ravel()[K:]: ax.set_visible(False)
    fig.suptitle("Basis functions")
    #if fname: plt.savefig(fname); print("saved", fname)
    #plt.tight_layout(); plt.show(); 
    return axes
