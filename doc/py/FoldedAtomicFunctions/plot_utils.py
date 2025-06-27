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

def plotMultiFunctionApprox(xs, data_pairs, bError=False, colors=None, errMax=None, scMin=None, title='Function Approximation', label_pairs=None, ax1=None):
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

    label_pairs : list[tuple[str, str]] | None
        Optional list of (ref_label, approx_label) tuples matching the order in *data_pairs*.

    ax1 : matplotlib.axes.Axes, optional
        Axis to plot on. If None, a new figure & axis are created.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object.
    axes : tuple of matplotlib.axes.Axes
        A tuple containing the main axis and the error axis (if bError is True).
    """
    if ax1 is None:
        fig, ax1 = plt.subplots(figsize=(10,6))
    else:
        fig = ax1.figure
    # Secondary axis for error if requested
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
        # Custom labels if provided
        if label_pairs and i < len(label_pairs):
            ref_lab, app_lab = label_pairs[i]
        else:
            ref_lab, app_lab = f'Sample {i+1}', f'Approx {i+1}'
        ax1.plot(xs, y_ref, label=ref_lab, ls=':', lw=1.5, c=color)
        ax1.plot(xs, y_app, label=app_lab, ls='-', lw=1.0, c=color)
        if bError:
            err_label = f'Err {i+1}' if not (label_pairs and i < len(label_pairs)) else f'Err {i+1}'
            ax2.plot(xs, y_ref - y_app, color=color, linestyle='--', lw=0.8, label=err_label)
    ax1.set_xlabel('z (Å)'); ax1.set_ylabel('Value'); ax1.set_title(title); ax1.grid(True, linestyle=':', alpha=0.7)
    if scMin is not None: vmin = np.min(np.concatenate([pair[0] for pair in data_pairs])); ax1.set_ylim(vmin * scMin, -vmin * scMin)
    if bError and errMax is not None: ax2.set_ylim(-errMax, errMax)
    lines1, labels1 = ax1.get_legend_handles_labels(); lines2, labels2 = ax2.get_legend_handles_labels() if (bError and ax2 is not None) else ([], [])
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='best', fontsize=8); fig.tight_layout()
    return fig, (ax1, ax2) if bError else (ax1,)

def plotFunctionApprox( xs, y_ref, ys_approx, bError=False, colors=None, errMax=1.0e-3, scMin=1.5 ):
    fig, ax1 = plt.subplots(figsize=(12,6))
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
    ax2.legend(lines + lines2, labels + labels2, loc='center right', fontsize=8)
    fig.tight_layout()
    #plt.savefig("exp_kz_approximation_sequence.png")
    #plt.show()
    return fig,(ax1,ax2) 

diverging_cmaps = {'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu', 'RdYlBu','RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic','berlin', 'managua', 'vanimo' }


def is_diverging_cmap(cmap_name):
    cmap = cmap_name.split('_')[0]
    return cmap in diverging_cmaps


# ===================================================================== # Keep existing function
# 2-D plotting – compact helpers (imshow based)
# =====================================================================

def imshow_grid(grid, extent=None, cmap="seismic", ax=None, figsize=(10, 4),title=None):
    if ax is None: fig, ax = plt.subplots(figsize=figsize)
    if is_diverging_cmap(cmap):
        vmin = np.nanmin(grid); vmax = -vmin; print("imshow_grid() vmin, vmax:", vmin, vmax, title)
    im = ax.imshow(grid.T, origin="lower", aspect="auto", extent=extent, cmap=cmap, vmin=vmin, vmax=vmax)
    plt.colorbar(im, ax=ax, shrink=0.8)
    ax.set_xlabel("x (Å)"); ax.set_ylabel("z (Å)")
    if title is not None: ax.set_title(title)
    fig.tight_layout()
    #if fname: plt.savefig(fname); print("saved", fname)
    #plt.tight_layout(); plt.show(); 
    return ax

def plot_atoms( apos, colors, sz=10, ax=None, bEqual=True, axes=(0,2), figsize=(10,4) ):
    if ax is None: fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(apos[:,axes[0]], apos[:,axes[1]], s=sz, c=colors, edgecolors="w", linewidths=0.5)
    ax.set_xlabel("x (Å)"); ax.set_ylabel("z (Å)") 
    if bEqual: ax.set_aspect("equal")
    #fig.tight_layout()
    #if fname: plt.savefig(fname); print("saved", fname)
    return ax

def plot2Dapprox(ref, fit, err=None, extent=None, title="Fit vs Ref", cmap="seismic", scErr=None):
    """Show reference, fit and error in one row of imshows."""
    if err is None: err = fit - ref
    fig, axs = plt.subplots(1, 3, figsize=(14, 4))
    if is_diverging_cmap(cmap):
        vmin = min(np.nanmin(ref),np.nanmin(fit)); vmax = -vmin; # print("imshow_grid() vmin, vmax:", vmin, vmax, title)
    im_ref = axs[0].imshow(ref.T, origin="lower", aspect="auto", extent=extent, cmap=cmap, vmin=vmin, vmax=vmax)
    im_fit = axs[1].imshow(fit.T, origin="lower", aspect="auto", extent=extent, cmap=cmap, vmin=vmin, vmax=vmax)
    if scErr is None: 
        scErr = max(-np.nanmin(err),np.nanmax(err));
    im_err = axs[2].imshow(err.T, origin="lower", aspect="auto", extent=extent, cmap="seismic", vmin=scErr, vmax=-scErr)
    axs[0].set_title("Reference")
    axs[1].set_title("Fit")
    axs[2].set_title("Error")
    # This colorbar is for error plot
    plt.colorbar(im_ref, ax=axs[0], shrink=0.7)
    plt.colorbar(im_fit, ax=axs[1], shrink=0.7)
    plt.colorbar(im_err, ax=axs[2], shrink=0.7)  
    fig.suptitle(title)
    #if fname:   plt.savefig(fname); print("saved", fname)
    #plt.tight_layout(); plt.show(); 
    return axs

def plot2Dbasis(cols, shape, extent, coeffs=None, labels=None, nrow=4, cmap="RdBu_r"):
    """Grid of imshows for basis columns (each column flattened)."""
    K = cols.shape[0]
    nrow = min(nrow, K)
    ncol = (K + nrow - 1) // nrow
    fig, axes = plt.subplots(nrow, ncol, figsize=(3*ncol, 3*nrow))
    axes = np.atleast_2d(axes)
    for k in range(K):
        ax = axes[k % nrow, k // nrow]
        im = ax.imshow(cols[k].reshape(shape).T, origin="lower", aspect="auto", extent=extent, cmap=cmap)
        t = labels[k] if labels else f"ϕ{k}"
        if coeffs is not None and k < len(coeffs):
            t += f"\nC={coeffs[k]:.1e}"
        ax.set_title(t, fontsize=7); ax.set_axis_off()
    for ax in axes.ravel()[K:]: ax.set_visible(False)
    fig.suptitle("Basis functions")
    #if fname: plt.savefig(fname); print("saved", fname)
    #plt.tight_layout(); plt.show(); 
    return axes


'''
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
'''