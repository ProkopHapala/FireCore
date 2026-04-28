import numpy             as np
import matplotlib.pyplot as plt
import func_utils as fu

# ==============================================================================
#  Reference Potential
# ==============================================================================

def Coulomb(r, A=1.0, b=1.0):
    """True Coulomb potential for reference."""
    r_safe = np.maximum(r, 1e-9)
    y = 1/r_safe
    E = A * y
    F = A * y * (b + 1/r_safe)
    return E, F, y

def Yukawa(r, A=1.0, b=1.0):
    """True Yukawa potential for reference."""
    r_safe = np.maximum(r, 1e-9)
    y      = np.exp(-b * r_safe) / r_safe
    E      = A * y
    F      = A * y * (b + 1/r_safe)
    return E, F, y

# ==============================================================================
#  Yukawa Approximation with Power-Transform and Automatic Cutoff
# ==============================================================================

def Yukawa_pow_cubic_autoRc(r, A, b, n, r1, r2):
    """
    Approximates the Yukawa potential using the power-transform method.
    (CORRECTED VERSION)
    """
    # 1. Define the target function g(r) and its derivative g'(r)
    def g(r_val):
        r_safe = np.maximum(r_val, 1e-9)
        return (np.exp(-b * r_safe) / r_safe)**(1/n)
        
    def g_prime(r_val):
        r_safe = np.maximum(r_val, 1e-9)
        return -g(r_safe) * (b + 1/r_safe) / n

    # 2. Get target values and derivatives at the fit points
    g1, d1 = g(r1), g_prime(r1)
    g2, d2 = g(r2), g_prime(r2)
    
    # 3. Set up and solve the 4x4 linear system for the coefficients of p(r)
    #    ***** THE FIX IS HERE: Renamed the matrix from 'A' to 'A_matrix' *****
    A_matrix = np.array([
        [r1**3,   r1**2,  r1, 1], [3*r1**2, 2*r1,   1,  0],
        [r2**3,   r2**2,  r2, 1], [3*r2**2, 2*r2,   1,  0]
    ])
    B = np.array([g1, d1, g2, d2])
    coeffs = np.linalg.solve(A_matrix, B)
    
    # 4. Find the natural cutoff Rc by finding the root of p(r)
    all_roots = np.roots(coeffs)
    real_roots = all_roots[np.isreal(all_roots)].real
    physical_roots = real_roots[real_roots > 0]
    if len(physical_roots) == 0:
        raise ValueError("Could not find a physical cutoff root for Yukawa.")
    Rc = np.min(physical_roots)

    # 5. Evaluate the base polynomial p(r) and its derivative
    deriv_coeffs = coeffs[:-1] * np.array([3, 2, 1])
    p_r = np.polyval(coeffs, r)
    p_prime_r = np.polyval(deriv_coeffs, r)
    
    p_r_clamped = np.maximum(0, p_r)
    
    # 6. Reconstruct the Yukawa term and the final Energy
    y_approx = p_r_clamped**n
    # Now this 'A' correctly refers to the scalar function argument.
    E = A * y_approx
    
    # 7. Calculate the Force using the chain rule
    de_dr = n * (p_r_clamped**(n - 1)) * p_prime_r
    F = -A * de_dr
    
    return E, F, p_r, Rc


if __name__ == "__main__":
    # Define the physical parameters
    xs = np.linspace(0.01, 20.0, 1000)
    A_val = 5.0  # Amplitude
    b_val = 0.05  # Screening parameter
    n_val = 8    # Integer power for the transform
    
    # Define fit points in absolute coordinates
    r1 = 4.0
    r2 = 7.0
    
    print("--- Model Parameters ---")
    print(f"Yukawa: A = {A_val}, b = {b_val}")
    print(f"Power-Transform: n = {n_val}")
    print(f"Fit Points (absolute r): r1={r1:.2f}, r2={r2:.2f}")

    # Call the function once to get the calculated Rc for the legend
    _, _, _, Rc_calc = Yukawa_pow_cubic_autoRc(xs, A_val, b_val, n_val, r1, r2)
    print(f"\n[INFO] Calculated Cutoff Rc = {Rc_calc:.4f}")

    functions_to_plot = [
        ('Coulomb (True)',        Coulomb,                   'gray', {'A': A_val, 'b': b_val}, None),
        ('Yukawa (True)',         Yukawa,                    'k',    {'A': A_val, 'b': b_val}, None),
        (f'PowCubic r2={r1+2.0:.2f}', Yukawa_pow_cubic_autoRc,   'r',    {'A': A_val, 'b': b_val, 'n':n_val, 'r1':r1, 'r2':r1+2.0}, None),
        (f'PowCubic r2={r1+3.0:.2f}', Yukawa_pow_cubic_autoRc,   'g',    {'A': A_val, 'b': b_val, 'n':n_val, 'r1':r1, 'r2':r1+3.0}, None),
        (f'PowCubic r2={r1+4.0:.2f}', Yukawa_pow_cubic_autoRc,   'b',    {'A': A_val, 'b': b_val, 'n':n_val, 'r1':r1, 'r2':r1+4.0}, None),
    ]

    # Plot functions and their derivatives
    fig, axs = fu.plot_funcs(
        functions_to_plot,
        xs,
        figsize=(10, 12),
        nderivs=1, 
        titles=['Potential Energy (E)', 'Force (F = -dE/dr)', "Approximation of the n-th root: p(r)"]

    )

    # --- Customize plots ---
    axs[0].set_ylim(-0.1,2.0)
    axs[0].axvline(Rc_calc, ls='--', c='r', label=f'Calculated Rc = {Rc_calc:.2f}')
    axs[0].legend()

    axs[0].axvline(r1, ls='--', c='g', label=f'Fit Point r1 = {r1:.2f}')
    axs[0].axvline(r2, ls='--', c='g', label=f'Fit Point r2 = {r1+2.0:.2f}')
    axs[0].axvline(r2, ls='--', c='g', label=f'Fit Point r2 = {r1+3.0:.2f}')
    axs[0].axvline(r2, ls='--', c='g', label=f'Fit Point r2 = {r1+4.0:.2f}')
    
    axs[1].axhline(0, ls=':', c='k', lw=0.5)
    axs[1].legend()
    axs[1].set_ylim(-0.1,1.0)

    # # On the 3rd plot, we validate the fit of p(r) to the true n-th root
    # target_g = (np.exp(-b_val * xs) / xs)**(1/n_val)
    # axs[2].plot(xs, target_g, 'k:', label=f'True (exp(-br)/r)^(1/n)')
    # axs[2].axvline(r1_val, ls=':', c='r', lw=1)
    # axs[2].axvline(r2_val, ls=':', c='r', lw=1)
    # axs[2].axhline(0, ls=':', c='k', lw=0.5)
    # axs[2].set_xlabel("Internuclear Distance (r)")
    # axs[2].set_ylim(-0.2, 1.5)
    # axs[2].legend()

    

    
    plt.tight_layout()
    plt.show()