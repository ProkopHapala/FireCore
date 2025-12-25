import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import convolve

class WaveletBasis:
    def __init__(self, name='coif2'):
        """
        Initialize with standard filter coefficients.
        Symlets and Coiflets are orthonormal.
        """
        self.name = name
        
        # Coefficients for Coiflet 2 (approximate, often called Coif4 in some libs)
        # These are the Low Pass Decomposition filters
        if name == 'coif2':
            # Length 12, fairly smooth, symmetric-ish
            self.h = np.array([
                -0.000720549445, -0.001823208871,  0.005611434819,  0.023680171947,
                -0.059434418646, -0.076488599079,  0.417005184424,  0.812723635450,
                 0.386110066823, -0.067372554722, -0.041464936782,  0.016387336465
            ])
        
        # Coefficients for Symlet 4 (Least Asymmetric Daubechies 4)
        elif name == 'sym4':
            self.h = np.array([
                -0.075765714789, -0.029637522770,  0.497618667632,  0.803738751805,
                 0.297857795605, -0.099219543576, -0.012609967269,  0.032223100604
            ])
            
        else:
            raise ValueError("Unknown wavelet. Try 'coif2' or 'sym4'")

        # In standard orthonormal wavelet theory, the reconstruction filter 
        # is the time-reverse of decomposition.
        self.g = self.h[::-1] # Reconstruction Low Pass (Scaling Function)
        
        # Compute the Laplacian Filter (Kinetic Energy Operator) numerically
        self.laplacian_kernel = self._compute_laplacian_stencil()

    def get_scaling_function(self, resolution=8):
        """
        Uses the Cascade Algorithm to compute the shape of the scaling function phi(x).
        resolution: number of iterations (higher = smoother curve)
        """
        # Start with a "delta" in coefficient space
        coeffs = np.array([1.0])
        
        for _ in range(resolution):
            # 1. Upsample (insert zeros between samples)
            upsampled = np.zeros(len(coeffs) * 2)
            upsampled[::2] = coeffs
            # 2. Convolve with reconstruction filter
            coeffs = np.convolve(upsampled, self.g) * np.sqrt(2)
            
        x = np.linspace(0, len(self.h) - 1, len(coeffs))
        return x, coeffs

    def _compute_laplacian_stencil(self):
        """
        Computes the matrix elements T_k = <phi'' | phi(x-k)>
        This creates a small filter we can convolve with to get Kinetic Energy.
        """
        # 1. Generate very high res scaling function
        _, phi = self.get_scaling_function(resolution=10)
        
        # 2. Compute 2nd derivative via Finite Differences on this fine grid
        dx = 1.0 / (2**10) # Grid spacing at this resolution
        d2phi = np.gradient(np.gradient(phi, dx), dx)
        
        # 3. Compute overlap integral <phi''(x) | phi(x-k)>
        # We do this by cross-correlating d2phi with phi
        # This results in a filter of size roughly 2*len(h)
        # We normalize by the fine grid spacing
        T_kernel = np.correlate(d2phi, phi, mode='full') * dx
        
        # Downsample back to integer grid shifts
        # The center of the correlation is the self-interaction (k=0)
        center_idx = len(T_kernel) // 2
        
        # We extract every (2^resolution)th point to get integer shifts
        step = 2**10
        
        # Extract the stencil (e.g., -2, -1, 0, 1, 2)
        # We grab enough points to cover the support
        num_coeffs = len(self.h)
        indices = np.arange(center_idx - (num_coeffs-1)*step, center_idx + (num_coeffs)*step, step)
        indices = indices[ (indices >= 0) & (indices < len(T_kernel)) ]
        
        stencil = T_kernel[indices]
        
        # The stencil must be symmetric for a Hermitian operator. 
        # Numerical noise might break it slightly, so enforce symmetry:
        stencil = 0.5 * (stencil + stencil[::-1])
        
        return stencil

# --- SIMULATION HELPERS ---

def grid_to_coeffs(func_on_grid):
    """
    Trivial Identity transform for this demo.
    In a 'Scaling Function Basis' code, the coefficients c_i 
    are essentially the values on the grid points for smooth functions.
    (Approximation holds well for Coiflets due to vanishing moments).
    """
    return func_on_grid

def compute_kinetic_energy(coeffs, wavelet):
    """
    Calculates <psi | -0.5 * d^2/dx^2 | psi>
    using the precomputed Laplacian convolution kernel.
    """
    # Convolution with periodic boundary conditions
    # T_psi = Laplacian * psi
    T_psi = convolve(coeffs, wavelet.laplacian_kernel, mode='same', method='direct')
    
    # Energy = sum( psi * (-0.5 * T_psi) )
    # Note: The kernel already contains the derivatives.
    # We add -0.5 factor for physics.
    kinetic_energy_density = coeffs * (-0.5 * T_psi)
    return np.sum(kinetic_energy_density)

def compute_potential_energy(coeffs, V_on_grid):
    """
    Calculates <psi | V | psi> using Dual Space approach.
    Since we are in the Scaling Function basis (uniform grid),
    and Coiflets are 'interpolating', this is simply a point-wise product.
    """
    # If using Coiflets, coeffs approx value at grid point.
    # E_pot = sum( |c_i|^2 * V_i )
    rho = coeffs**2
    return np.sum(rho * V_on_grid)

# ==========================================
# MAIN EXECUTION
# ==========================================

# 1. Setup
N = 100  # Grid points
x_grid = np.linspace(-10, 10, N)
dx = x_grid[1] - x_grid[0]

# Create Wavelet Solvers
sym4 = WaveletBasis('sym4')
coif2 = WaveletBasis('coif2')

# 2. Plotting the Scaling Functions (Visual Check)
plt.figure(figsize=(10, 5))

x_s, y_s = sym4.get_scaling_function()
# Center the plot
plt.plot(x_s - np.mean(x_s), y_s, label='Symlet 4 (Least Asymmetric)', lw=2)

x_c, y_c = coif2.get_scaling_function()
plt.plot(x_c - np.mean(x_c), y_c, label='Coiflet 2 (Smoother)', lw=2, linestyle='--')

plt.title("Scaling Functions $\phi(x)$ generated via Cascade Algorithm")
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()

# 3. Physics Demo: Harmonic Oscillator
# Let's define a Gaussian Wavepacket (Ground state approximation)
sigma = 1.0
psi_real = np.exp(-x_grid**2 / (2*sigma**2))
psi_real /= np.sqrt(np.sum(psi_real**2)) # Normalize on grid (sum c^2 = 1)

# In scaling function basis, coefficients ~= function values (normalized)
coeffs = psi_real 

# Define Potential V = 0.5 * k * x^2
k_spring = 1.0
V_grid = 0.5 * k_spring * x_grid**2

# Calculate Energies using Coiflets
# Note: We must scale the Laplacian Kernel by 1/dx^2 because our grid 
# is not spaced by 1.0, but by dx.
# The stencil was computed for unit spacing.
T_val = compute_kinetic_energy(coeffs, coif2) / (dx**2)
V_val = compute_potential_energy(coeffs, V_grid)
Total_E = T_val + V_val

print(f"--- Simulation Results (Harmonic Oscillator) ---")
print(f"Grid Spacing (dx): {dx:.4f}")
print(f"Kinetic Energy:    {T_val:.6f} Ha (Exact: 0.25)")
print(f"Potential Energy:  {V_val:.6f} Ha (Exact: 0.25)")
print(f"Total Energy:      {Total_E:.6f} Ha (Exact: 0.50)")
print("\nNote: Accuracy improves with Coiflet order and Grid density.")

# 4. Show the Laplacian Filter (The "Stencil")
print(f"\nComputed Laplacian Stencil for Coiflet (size {len(coif2.laplacian_kernel)}):")
print(np.round(coif2.laplacian_kernel, 4))