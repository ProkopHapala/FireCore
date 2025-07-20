import numpy as np
import matplotlib.pyplot as plt

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf
import math

# --- Constants for the Optimized Piecewise Approximation ---
# These values are derived by:
# 1. Matching function value at x=0: f1(0) = B(0) => c = 2/sqrt(pi)
# 2. Matching function value at r_min: f1(r_min) = f2(r_min) => a*r_min^2 + c = 1/r_min
# 3. Matching derivative at r_min: f1'(r_min) = f2'(r_min) => 2*a*r_min = -1/r_min^2
# Solving these three equations simultaneously gives the unique r_min, a, c.

R_MIN = (3 * math.sqrt(math.pi)) / 4
C_APPROX = 2 / math.sqrt(math.pi)
A_APPROX = -1 / (2 * R_MIN**3)

# --- Boys Function and its Analytical Derivative ---

def boys_function(r):
    """
    Evaluates the Boys function B(r) = erf(r)/r and its derivative.
    Handles the singularity at r=0 using masks.
    Returns: (E, F) where E is function values, F is derivatives
    """
    r = np.asarray(r)
    E = np.zeros_like(r)
    F = np.zeros_like(r)
    
    mask = r != 0
    r_masked = r[mask]
    
    # Compute for non-zero values
    erf_r = erf(r_masked)
    exp_r = np.exp(-r_masked**2)
    E[mask] = erf_r / r_masked
    F[mask] = (r_masked * (2/np.sqrt(np.pi) * exp_r) - erf_r) / r_masked**2
    
    # Handle r=0 case
    E[~mask] = 2/np.sqrt(np.pi)
    F[~mask] = 0.0
    
    return E, F

def piecewise_approximation(r):
    """
    Evaluates the piecewise approximation and its derivative.
    Returns: (E, F) where E is function values, F is derivatives
    """
    r = np.asarray(r)
    E = np.zeros_like(r)
    F = np.zeros_like(r)
    
    mask = np.abs(r) < R_MIN
    
    # Compute for |r| < R_MIN
    E[mask] = A_APPROX * r[mask]**2 + C_APPROX
    F[mask] = 2 * A_APPROX * r[mask]
    
    # Compute for |r| >= R_MIN
    E[~mask] = 1 / np.abs(r[~mask])
    F[~mask] = -1 / r[~mask]**2
    
    return E, F


# --- Numerical Derivative Function ---

def numDeriv(x, y):
    """
    Calculates the numerical derivative using central difference.
    x: array of x values
    y: array of corresponding y values
    Returns: (numerical_derivative_values, x_values_for_derivative)
    """
    dx = x[2:] - x[:-2]
    dy = y[2:] - y[:-2]
    x_ = x[1:-1]
    return dy / dx, x_

# --- Modular Plotting Functions ---

def plot_with_deriv(ax1, ax2, x, y, y_deriv, label, color, linestyle='-'):
    ax1.plot(x, y, label=label, color=color, linestyle=linestyle)
    ax2.plot(x, y_deriv, label=label, color=color, linestyle=linestyle)

# --- Plotting ---

def plot_comparison():
    # Generate x values
    r_values = np.linspace(0.0, 5.0, 1000)

    boys_vals, boys_deriv_vals     = boys_function(r_values)
    approx_vals, approx_deriv_vals = piecewise_approximation(r_values)

    # Create figure
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), sharex=True)
    fig.suptitle('Boys Function and its Piecewise Parabolic Approximation', fontsize=16)

    # Plot using modular functions
    plot_with_deriv(ax1, ax2, r_values, boys_vals, boys_deriv_vals, label='Boys Function (erf(r)/r)', color='blue')
    
    plot_with_deriv(ax1, ax2, r_values, approx_vals, approx_deriv_vals, label='Piecewise Approximation', color='red', linestyle='--')

    # Add vertical lines at R_MIN
    ax1.axvline( R_MIN,  color='gray',  linestyle=':', label=f'r_min = {R_MIN:.3f}')
    #ax1.axvline(-R_MIN, color='gray', linestyle=':')

    # Customize axes
    ax1.set_ylabel('Function Value (Energy)',  fontsize=12)
    ax2.set_ylabel('Derivative Value (Force)', fontsize=12)
    ax2.set_xlabel('r', fontsize=12)

    # Add legends and grid
    ax1.legend()
    ax2.legend()
    ax1.grid(True)
    ax2.grid(True)

    plt.tight_layout(rect=[0, 0.03, 1, 0.96])
    plt.show()

if __name__ == "__main__":

    plot_comparison()