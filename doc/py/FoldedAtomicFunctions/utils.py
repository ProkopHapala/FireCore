# --- 5. Plotting & Analysis Utilities ---

import numpy as np
import matplotlib.pyplot as plt

def plot_1d_profiles(
    zs: np.ndarray,
    profiles_T: np.ndarray, # (num_profiles, Nz)
    title: str,
    max_plot: int = 10,
    ws: np.ndarray = None,
    potential_plot_yscale_factor: float = 1.5,
    filename: str = None
):
    """Plots 1D profiles."""
    if profiles_T.size == 0:
        print(f"No profiles to plot for: {title}")
        return
        
    num_profiles = profiles_T.shape[0]
    plt.figure(figsize=(10, 6))
    
    actual_profiles_to_plot = profiles_T[:min(num_profiles, max_plot), :]
    for i in range(actual_profiles_to_plot.shape[0]):
        plt.plot(zs, actual_profiles_to_plot[i, :], label=f'Profile {i+1}', alpha=0.7)
    
    if ws is not None:
        # Scale weights to fit nicely on the plot
        ax2 = plt.gca().twinx()
        min_prof_disp = np.min(actual_profiles_to_plot) if actual_profiles_to_plot.size > 0 else 0
        max_prof_disp = np.max(actual_profiles_to_plot) if actual_profiles_to_plot.size > 0 else 1
        # Plot weights such that they are visible, e.g., in top 20% of y-axis range
        # Ensure plot_ws is calculated robustly even if max_prof_disp == min_prof_disp
        range_prof_disp = max_prof_disp - min_prof_disp
        if range_prof_disp < 1e-9 : range_prof_disp = 1.0 # Avoid division by zero or tiny range

        plot_ws = min_prof_disp + 0.8 * range_prof_disp + 0.2 * range_prof_disp * (ws / np.max(ws) if np.max(ws) > 1e-9 else ws)
        ax2.plot(zs, ws, 'k--', label='Weights (scaled)', alpha=0.5, linewidth=1)
        ax2.set_ylabel('Weights (scaled)')
        ax2.tick_params(axis='y')

    # Y-axis scaling for potentials
    if actual_profiles_to_plot.size > 0:
        overall_min_val = np.min(actual_profiles_to_plot)
        if overall_min_val < -1e-9: # If there's a negative minimum (potential well)
            plot_vmin = overall_min_val * potential_plot_yscale_factor
            plot_vmax = -plot_vmin
            plt.ylim(plot_vmin, plot_vmax)

    plt.title(title)
    plt.xlabel('z (Å)')
    plt.ylabel('Potential / Value')
    if actual_profiles_to_plot.shape[0] > 0 : plt.legend() # Legend for main plot
    plt.grid(True)
    if filename:
        plt.savefig(filename)
        print(f"Plot saved to {filename}")
    plt.show()

def plot_singular_values(s_vals, K_opt, filename: str = None):
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
    if filename:
        plt.savefig(filename)
        print(f"Plot saved to {filename}")
    plt.show()


def print_analytical_form_polynomial(
    U_k_coeffs: np.ndarray,  # (P, K)
    z_scale_info: tuple | None,  # (z_min, z_range)
    basis_labels: list | None = None,
    K_to_print: int = -1,
):
    """Pretty-print analytical form of optimal basis if *Phi* was polynomial.

    The implementation is intentionally concise and skips negligible terms.
    """
    if U_k_coeffs.size == 0:
        return

    P, K_actual = U_k_coeffs.shape
    if K_to_print < 0 or K_to_print > K_actual:
        K_to_print = K_actual

    if z_scale_info:
        z_min, z_range = z_scale_info
        print(
            f"Note: z_scaled = (z − {z_min:.3f}) / {z_range:.3f} (used in original basis)"
        )
    else:
        print("Note: original z-coordinates used (no scaling).")

    for k in range(K_to_print):
        parts: list[str] = []
        for p in range(P):
            c = U_k_coeffs[p, k]
            if abs(c) < 1e-6:
                continue
            label = basis_labels[p] if basis_labels and p < len(basis_labels) else f"phi_{p}"
            parts.append(f"({c:+.3e} * {label})")
        print(f"B_opt_{k+1}(z) = " + " ".join(parts))

def plotFunctionApprox( xs, y_ref, ys_approx, bError=False, colors=None, errMax=1.0e-3 ):
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
        #poly_approx, z_cut_n = eval_poly_at_z0(zs, k0, n_val, z0_approx)
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
    ax1.set_ylim(-0.1, 1.2)
    ax1.set_xlim(xs[0], xs[-1])
    ax1.grid(True, linestyle=':', alpha=0.7)
    lines, labels = ax1.get_legend_handles_labels()
    if bError:
        lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, loc='center right')
    fig.tight_layout()
    #plt.savefig("exp_kz_approximation_sequence.png")
    #plt.show()
    return ax1, ax2