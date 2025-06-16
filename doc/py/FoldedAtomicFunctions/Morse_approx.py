import numpy as np
import matplotlib.pyplot as plt
# Assuming potentials.py is in the same directory or accessible via PYTHONPATH
from potentials import morse_potential
from utils import plotFunctionApprox

# infinite numpy line length
np.set_printoptions(linewidth=np.inf)


def get_base_zero_arg(k: float, n: int) -> float:
    """
    Calculates the value of z_arg where the base (1 - k*z_arg/(2n)) becomes zero.
    This is the "natural" cutoff for the (1 - k*z_arg/(2n))**(2n) form.
    """
    if k <= 0: raise ValueError("k must be positive.")
    if n <= 0: raise ValueError("n must be a positive integer.")
    return (2.0 * n) / k

def get_zcut(k: float, n: int, z0_approx_point: float) -> float:
    """Calculates z_cut for polynomial (z_cut - z)^(2n) to match exp(-kz) decay at z0_approx_point."""
    if k <= 0: raise ValueError("k must be positive.")
    if n <= 0: raise ValueError("n must be a positive integer.")
    return z0_approx_point + get_base_zero_arg(k, n)

def eval_poly_exp_native_form(z_arg: np.ndarray | float, k: float, n: int) -> tuple[np.ndarray | float, float]:
    """
    Calculates the polynomial approximation term (max(0, 1 - k_decay*z_arg/(2*n_factor)))**(2*n_factor).
    This approximates exp(-k_decay*z_arg).

    # Parameters:
    z_arg : np.ndarray    Argument of the exponential decay (e.g., 'z' or 'z-r0').
    k : float             The decay constant 'k'.
    n : int               The 'n' in the polynomial (1 - kz/(2n))**(2n).

    # Returns:
    approx_values : The polynomial approximation.
    z_cut_arg :  The value of z_arg where the base (1 - k*z_arg/(2n)) becomes zero.         z_cut_arg = (2 * n_factor) / k_decay.
    """
    if k <= 0: raise ValueError("k must be positive.")
    if n <= 0: raise ValueError("n must be a positive integer.")
    term_base = 1.0 - (k * z_arg) / (2.0 * n)
    approx_values = np.maximum(0, term_base)**(2 * n)
    z_cut_arg = get_base_zero_arg(k, n) # Natural cutoff for this form's argument
    return approx_values, z_cut_arg

def eval_poly_at_z0(zs: np.ndarray, k: float, n: int, z0_approx_point: float) -> tuple[np.ndarray | float, float]:
    """
    Evaluates polynomial C*(z_cut - z)^(2n) to approximate exp(-kz) optimally at z0_approx_point.

    Parameters:
    -----------
    zs : np.ndarray
        Array of z-coordinates for evaluation.
    k : float
        Decay constant in exp(-kz).
    n : int
        Power factor 'n' in (z_cut - z)^(2n).
    z0_approx_point : float
        The z-coordinate around which the approximation is optimized.

    Returns:
    --------
    poly_values : np.ndarray
        Values of the polynomial approximation.
    z_cut : float
        The calculated z_cut for the polynomial.
    """
    z_cut = get_zcut(k, n, z0_approx_point)
    
    # C = exp(-k*z0) / (z_cut - z0)^(2n) = exp(-k*z0) / ( (2n)/k )**(2n)
    C_coeff = np.exp(-k * z0_approx_point) / (get_base_zero_arg(k, n)**(2 * n))
    
    poly_values = C_coeff * (np.maximum(0, z_cut - zs)**(2 * n))
    return poly_values, z_cut

def get_taylor_coeffs(k: float, z0: float, order: int) -> np.ndarray:
    """
    Calculates Taylor expansion coefficients for exp(-kz) around z0 up to a given order.
    The polynomial is in terms of (z-z0).
    Coefficient c_i corresponds to the term (z-z0)^i / i!
    So, c_i = f^(i)(z0) = (-k)^i * exp(-kz0)
    The full term is (c_i / i!) * (z-z0)^i.
    This function returns coeffs_for_poly1d = [c_order/order!, c_{order-1}/(order-1)!, ..., c_1/1!, c_0/0!]
    for np.poly1d, which expects coefficients from highest power to lowest.
    """
    if order < 0: raise ValueError("Order must be non-negative.")
    exp_kz0 = np.exp(-k * z0)
    # Coefficients for P(x) = a_n x^n + ... + a_1 x + a_0
    # where x = (z-z0)
    coeffs_for_poly1d = [((-k)**i * exp_kz0) / np.math.factorial(i) for i in range(order)]
    return np.array(coeffs_for_poly1d)

def eval_taylor_exp_approx(zs: np.ndarray, k: float, z0: float, order: int) -> np.ndarray:
    """
    Evaluates Taylor expansion of exp(-kz) around z0 up to a given order.
    f(z) approx sum_{i=0 to order} [ f^(i)(z0) * (z-z0)^i / i! ]
    f^(i)(z0) = (-k)^i * exp(-kz0)
    """
    if order < 0: raise ValueError("Order must be non-negative.")
    
    # Get coefficients for polynomial in (z-z0)
    # coeffs are [c_order, c_{order-1}, ..., c_0] where c_i is for (z-z0)^i
    taylor_coeffs = get_taylor_coeffs(k, z0, order)
    print(f"Taylor O({order}) at z0={z0:.1f}: Coeffs (for (z-z0)^i terms): {taylor_coeffs}")
    poly = np.poly1d(taylor_coeffs[::-1])
    return poly(zs - z0)

def get_taylor_approx_data_for_plot(zs, k, z0, order_seq):
    """
    Helper to generate data for plotFunctionApprox for Taylor series.
    The 'z_cut' equivalent for Taylor isn't a sharp cutoff, so we pass z0
    as a reference point if the plotting function needs a coordinate.
    """
    ys_approx_taylor = []
    for order_val in order_seq:
        taylor_approx_values = eval_taylor_exp_approx(zs, k, z0, order_val)
        # Get and print coefficients for inspection
        coeffs = get_taylor_coeffs(k, z0, order_val)
        # print(f"Taylor O({order_val}) at z0={z0:.1f}:")
        # print(f"  Coeffs (for (z-z0)^i terms, highest power first): {coeffs}")
        # print(f"  Evaluated: {taylor_approx_values[:3]}...") # Print a few values
        # For plotFunctionApprox, the second element in the tuple is z_cut.
        # For Taylor, there isn't a direct equivalent. We can pass z0 or None.
        # Let's pass z0 as it's the expansion point.
        ys_approx_taylor.append((taylor_approx_values, z0, f'Taylor O({order_val})'))
    return ys_approx_taylor

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
    
    n_seq = [1, 8, 10, 20 ] # Sequence of n for (1-kz/(2n))^(2n)
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
    n_seq_morse = [1, 5, 10, 20] # Different n values for Morse
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
