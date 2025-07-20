import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf
import math

# --- Constants ---
# Arbitrary r_min for demonstration. You can change this value.
# R_MIN = 1.5

# --- Coefficients for C2-Smooth Symmetric Approximation (r^6 + r^4 + r^2 + const) ---

SQRT_PI = math.sqrt(math.pi)
boys_r0 = 2 / SQRT_PI

def boys_approx_sympy():
    import sympy
    import math

    r_min = sympy.Symbol('r_min')
    c6, c4, c2 = sympy.symbols('c6 c4 c2')

    c0_val = 2 / sympy.sqrt(sympy.pi)

    f   = 1/r_min
    df  = sympy.diff(f, r_min)
    ddf = sympy.diff(df, r_min)

    print(f"f   = {f}")
    print(f"df  = {df}")
    print(f"ddf = {ddf}")

    poly = c6 * r_min**6 + c4 * r_min**4 + c2 * r_min**2 + c0_val
    d_poly  = sympy.diff(  poly, r_min)
    dd_poly = sympy.diff(d_poly, r_min)

    print(f"poly = {poly}")
    print(f"d_poly = {d_poly}")
    print(f"dd_poly = {dd_poly}")

    # Equations
    eq1 = sympy.Eq(     poly, f   )
    eq2 = sympy.Eq(   d_poly, df  )
    eq3 = sympy.Eq(  dd_poly, ddf )

    # Solve the system
    solution = sympy.solve([eq1, eq2, eq3], (c6, c4, c2))

    print(f"c6 = {solution[c6]}")
    print(f"c4 = {solution[c4]}")
    print(f"c2 = {solution[c2]}")


# --- Boys Function and its Analytical Derivatives (for r > 0) ---
def boys_function_positive_r(r_values):
    """
    Evaluates the Boys function B(r) = erf(r)/r and its derivatives for r > 0.
    Handles r=0 case if present in r_values by using the limit.
    Returns: (E, F, S) where E is function values, F is first derivatives, S is second derivatives.
    """
    r_values = np.asarray(r_values)
    E = np.zeros_like(r_values)
    F = np.zeros_like(r_values)
    S = np.zeros_like(r_values)

    mask_nonzero = r_values > 0
    r_masked = r_values[mask_nonzero]

    # Function value: erf(r)/r
    erf_r = erf(r_masked)
    E[mask_nonzero] = erf_r / r_masked

    # First derivative: (r * (2/sqrt(pi) * exp(-r^2)) - erf(r)) / r^2
    exp_r_sq = np.exp(-r_masked**2)
    F[mask_nonzero] = (r_masked * (2/SQRT_PI * exp_r_sq) - erf_r) / r_masked**2

    # Second derivative: d/dr(F)
    # This is more complex. Let's use sympy for verification if needed, or stick to numerical if not exact.
    # d/dr [ (2/sqrt(pi) * exp(-r^2))/r - erf(r)/r^2 ]
    # d/dr [ C * exp(-r^2)/r ] = C * [ -2r*exp(-r^2)*r - exp(-r^2) ] / r^2 = C * exp(-r^2) * (-2r^2 - 1) / r^2
    # d/dr [ -erf(r)/r^2 ] = -[ (2/sqrt(pi)*exp(-r^2))*r^2 - erf(r)*2r ] / r^4 = -[ (2/sqrt(pi)*exp(-r^2)) - 2*erf(r)/r ] / r^2
    # S = (2/SQRT_PI * exp_r_sq * (-2*r_masked**2 - 1) / r_masked**2) - ((2/SQRT_PI * exp_r_sq) - 2 * erf_r / r_masked) / r_masked**2
    # Let's simplify and use the common form for B''(r)
    S[mask_nonzero] = (2 * erf_r / r_masked**3) - (4/SQRT_PI * exp_r_sq / r_masked) - (4/SQRT_PI * exp_r_sq * r_masked)

    # Handle r=0 case if it's explicitly passed (though we only plot r>0)
    mask_zero = r_values == 0
    E[mask_zero] = boys_r0
    F[mask_zero] = 0.0
    S[mask_zero] = 0.0 # Limit of B''(r) as r->0 is 0

    return E, F, S

# --- C2-Smooth Symmetric Approximation and its Analytical Derivatives (for r > 0) ---
def c2_symmetric_approximation(r_values, R_MIN=1.5):
    """
    Evaluates the C2-smooth symmetric approximation and its derivatives for r > 0.
    Returns: (E, F, S) where E is function values, F is first derivatives, S is second derivatives.
    """

    C0 = boys_r0
    C6 = (15 * SQRT_PI - 16 * R_MIN) / (8 * SQRT_PI * R_MIN**7)
    C4 = (24 * R_MIN - 21 * SQRT_PI) / (4 * SQRT_PI * R_MIN**5)
    C2 = (35 * SQRT_PI - 48 * R_MIN) / (8 * SQRT_PI * R_MIN**3)

    r_values = np.asarray(r_values)
    E = np.zeros_like(r_values)
    F = np.zeros_like(r_values)
    S = np.zeros_like(r_values)

    mask_inner = (r_values > 0) & (r_values < R_MIN)
    mask_outer = r_values >= R_MIN
    
    r_inner = r_values[mask_inner]
    r_outer = r_values[mask_outer]

    # For 0 < r < R_MIN (polynomial region)
    E[mask_inner] =      C6 * r_inner**6 +      C4 * r_inner**4 +     C2 * r_inner**2 + C0
    F[mask_inner] = 6  * C6 * r_inner**5 +  4 * C4 * r_inner**3 + 2 * C2 * r_inner
    S[mask_inner] = 30 * C6 * r_inner**4 + 12 * C4 * r_inner**2 + 2 * C2

    # For r >= R_MIN (1/r region)
    E[mask_outer] =  1 / r_outer
    F[mask_outer] = -1 / r_outer**2
    S[mask_outer] =  2 / r_outer**3

    # Handle r=0 if explicitly in r_values (though plotting for r>0)
    mask_zero = r_values == 0
    E[mask_zero] = C0
    F[mask_zero] = 0.0
    S[mask_zero] = 0.0

    return E, F, S

if __name__ == "__main__":

    import sys

    R_MIN = 1.5
    if len(sys.argv) > 1: R_MIN = float(sys.argv[1])

    boys_approx_sympy()

    # Generate x values for r > 0
    r_values = np.linspace(0, 4.0, 500) # Start from 0 to show the match at r=0

    # Evaluate Boys function and approximation
    boys_E, boys_F, boys_S = boys_function_positive_r(r_values)
    approx_E, approx_F, approx_S = c2_symmetric_approximation(r_values, R_MIN)

    

    # Create figure
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 12), sharex=True)
    fig.suptitle(f'C2-Smooth Symmetric Potential Approximation (r_min = {R_MIN:.3f})', fontsize=16)

    # Plot Potential (Energy)
    ax1.plot(r_values, boys_E, label='Boys Function (erf(r)/r)', color='black', linewidth=2)
    ax1.plot(r_values, approx_E, label='C2 Symmetric Approx', color='blue', linestyle='--')
    ax1.axvline(R_MIN, color='red', linestyle=':', label=f'r_min = {R_MIN:.3f}')
    ax1.set_ylabel('V(r)', fontsize=12)
    ax1.legend()
    ax1.grid(True)
    ax1.set_ylim(bottom=0) # Potential should be positive

    # Plot Force (-dV/dr)
    ax2.plot(r_values, -boys_F, label='Force [Boys]', color='black', linewidth=2)
    ax2.plot(r_values, -approx_F, label='Force [C2 Symmetric Approx]', color='blue', linestyle='--')
    ax2.axvline(R_MIN, color='red', linestyle=':')
    ax2.set_ylabel('Force (-dV/dr)', fontsize=12)
    ax2.legend()
    ax2.grid(True)
    # Add a horizontal line at 0 for force
    ax2.axhline(0, color='gray', linestyle='--', linewidth=0.5)

    # Plot Force Derivative (-d^2V/dr^2)
    ax3.plot(r_values, -boys_S, label='Force Derivative [Boys]', color='black', linewidth=2)
    ax3.plot(r_values, -approx_S, label='Force Derivative [C2 Symmetric Approx]', color='blue', linestyle='--')
    ax3.axvline(R_MIN, color='red', linestyle=':')
    ax3.set_ylabel('Force Derivative (-d²V/dr²)', fontsize=12)
    ax3.set_xlabel('r', fontsize=12)
    ax3.legend()
    ax3.grid(True)
    # Add a horizontal line at 0 for force derivative
    ax3.axhline(0, color='gray', linestyle='--', linewidth=0.5)

    plt.tight_layout(rect=[0, 0.03, 1, 0.96])
    plt.show()