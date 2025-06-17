import numpy as np
import matplotlib.pyplot as plt
# Assuming potentials.py is in the same directory or accessible via PYTHONPATH
from potentials import morse_potential
from plot_utils import plotFunctionApprox
from basis_utils import eval_poly_exp_native_form, eval_poly_at_z0, get_taylor_approx_data_for_plot

# infinite numpy line length
np.set_printoptions(linewidth=np.inf)


if __name__ == "__main__":
    # --- General Parameters ---
    z_min, z_max = 0.0, 10.0
    nPts = 400
    zs = np.linspace(z_min, z_max, nPts)

    # --- Morse Potential Parameters ---
    D_m     = 0.1  # Depth (eV)
    alpha_m = 1.5 # Width (1/Angstrom)
    r0_eq   = 2.5 # Equilibrium distance (Angstrom)

    # =========================================================
    # --- Step 1: Approximate exp(-kz) by Taylor expansion ---
    # =========================================================
    print("--- Step 1: Approximating exp(-kz) by Taylor expansion ---")
    k0_taylor = 1.0
    z0_taylor_approx = r0_eq # Expand around r0_eq
    y_exp_taylor = np.exp(-k0_taylor * zs)
    #order_seq_taylor = [0, 2, 4, 6, 8] # Orders of Taylor polynomial
    order_seq_taylor = [ 1, 3, 5, 7] # Orders of Taylor polynomial

    print(f"\nApproximating exp(-{k0_taylor:.1f}z) around z0={z0_taylor_approx:.1f} using Taylor expansion:")
    
    ys_approx_taylor_data = get_taylor_approx_data_for_plot(zs, k0_taylor, z0_taylor_approx, order_seq_taylor)
    
    ax1_taylor, ax2_taylor = plotFunctionApprox(zs, y_exp_taylor, ys_approx_taylor_data, True, errMax=0.1) # errMax can be adjusted
    ax1_taylor.axvline(z0_taylor_approx, color='magenta', linestyle='--', linewidth=1.5, label=f'z0_expansion={z0_taylor_approx:.1f}')
    ax1_taylor.set_title(f'Taylor Approx of exp(-{k0_taylor:.1f}z) around z0={z0_taylor_approx:.1f}')
    ax1_taylor.set_ylim(-0.2, 1.5) # Adjust ylim for Taylor plot visibility
    #plt.savefig("taylor_exp_kz_approximation.png")
    plt.show()

    # =========================================================
    # --- Step 2: Approximate exp(-kz) by (1-kz/(2n))^(2n) ---
    # =========================================================
    print("--- Step 2: Approximating exp(-kz) ---")
    k0 = 1.0  # Example decay constant
    
    y_exp = np.exp(-k0 * zs)
    
    n_seq = [1,2,4,8,16] # Sequence of n for (1-kz/(2n))^(2n)
    #z0_approx = 3.0 # Point where we want the approximation to be good
    z0_approx = r0_eq

    print(f"\nApproximating exp(-{k0:.1f}z) around z0={z0_approx:.1f} using C*(z_cut-z)^(2n):")
    print("  where z_cut = z0 + 2n/k")

    ys_approx = [ (eval_poly_at_z0(zs, k0, n_val, z0_approx)) + (f'Approx n={n_val}',) for n_val in n_seq ]
    ax1,ax2 = plotFunctionApprox(zs, y_exp, ys_approx, True)
    ax1.axvline(z0_approx, color='magenta', linestyle='--', linewidth=1.5, label=f'z0_approx={z0_approx:.1f}')
    ax1.set_title(f'Approximation of exp(-{k0:.1f}z) around z0={z0_approx:.1f}')
    #plt.savefig("exp_kz_approximation_sequence.png")
    #plt.show()

    # ===========================================
    # --- Step 3: Approximate Morse Potential ---
    # ===========================================
    print("\n--- Step 3: Approximating Morse Potential ---")
    # Morse Potential Parameters
    # Calculate original Morse potential
    V_m = morse_potential(zs, D_m, alpha_m, r0_eq)

    # Use a sequence of n values for Morse approximation
    n_seq_morse = [1, 2, 4, 8, 16] # Different n values for Morse
    morse_approx_data = []

    print(f"\nApproximating Morse potential using V_app = D*(X_app^2 - 2*X_app):")
    print(f"  where X_app = (1 - alpha_m*(z-r0_eq)/(2n))^(2n)")

    for n_val_morse in n_seq_morse:
        z_arg_morse        = zs - r0_eq
        X_app, z_cut_arg_X = eval_poly_exp_native_form(z_arg_morse, alpha_m, n_val_morse)
        z_cut_for_X_abs    = r0_eq + z_cut_arg_X
        V_m_app_n          = D_m * (X_app**2 - 2 * X_app)
        morse_approx_data.append( (V_m_app_n, z_cut_for_X_abs, f'Approx n={n_val_morse}'))
        print(f"  n={n_val_morse}: z_cut for X_app at z={z_cut_for_X_abs:.3f} Å")

    ax1_morse, ax2_morse = plotFunctionApprox(zs, V_m, morse_approx_data, True, errMax=0.1*D_m if D_m > 0 else 0.01)
    ax1_morse.set_title(f'Morse Approx (D={D_m}, α={alpha_m}, r0={r0_eq})')
    ax1_morse.set_xlabel('z (Å)')
    ax1_morse.set_ylabel('Potential (eV)')
    ax1_morse.axvline(r0_eq, color='green', linestyle=':', linewidth=1.5, label=f'r0_eq={r0_eq:.1f}')
    
    # Y-axis limits: vmin = 1.5 * (-D_morse), vmax = -vmin
    plot_vmin = -1.5 * D_m
    plot_vmax = 1.5 * D_m
    ax1_morse.set_ylim(plot_vmin, plot_vmax)
        
    plt.savefig("morse_potential_polynomial_approx.png")
    plt.show()

    print("\nScript completed. Check generated PNG files.")
