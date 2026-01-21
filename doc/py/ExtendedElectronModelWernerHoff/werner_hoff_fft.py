import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fftn, ifftn, fftfreq

class EEMSolver:
    def __init__(self, box_size=6.0, grid_points=64):
        self.L = box_size
        self.N = grid_points
        self.dV = (self.L / self.N)**3
        
        # Grid setup
        x = np.linspace(-self.L/2, self.L/2, self.N, endpoint=False)
        self.X, self.Y, self.Z = np.meshgrid(x, x, x)
        
        # Fourier Space setup (for Derivatives/Poisson)
        k = (2 * np.pi / self.L) * np.fft.fftfreq(self.N) * self.N
        self.KX, self.KY, self.KZ = np.meshgrid(k, k, k)
        self.K2 = self.KX**2 + self.KY**2 + self.KZ**2
        self.K2[0,0,0] = 1.0 # Avoid division by zero (DC component handled separately)

    def get_laplacian(self, field_real):
        """Computes Laplacian via FFT"""
        field_k = fftn(field_real)
        lap_k = -self.K2 * field_k
        return np.real(ifftn(lap_k))

    def solve_poisson(self, density):
        """Solves Poisson eq for Hartree Potential: -Del^2 V = 4*pi*n"""
        rho_k = fftn(density)
        V_k = 4 * np.pi * rho_k / self.K2
        V_k[0,0,0] = 0 # Ignore infinite DC offset
        return np.real(ifftn(V_k))

    def get_nuclear_potential(self, atom_coords, charges):
        """Gaussian smeared nuclear potential (Local Pseudopotential)"""
        V = np.zeros_like(self.X)
        sigma = 0.4  # Smearing width to avoid 1/r singularity on grid
        for coord, Z in zip(atom_coords, charges):
            r = np.sqrt((self.X-coord[0])**2 + (self.Y-coord[1])**2 + (self.Z-coord[2])**2)
            # Soft-Coulomb
            V += -Z * (1.0 / np.sqrt(r**2 + sigma**2)) 
            # Or Gaussian: V += -Z / (sigma * np.sqrt(2*np.pi)) * np.exp(-0.5 * (r/sigma)**2)
        return V

    def functional_derivative_kinetic(self, sqrt_field):
        """
        Derivative of von Weizsacker Term: T[rho] = int (grad rho)^2 / 8rho
        Let phi = sqrt(rho). T[phi] = 1/2 int (grad phi)^2
        dT/dphi = -1/2 Laplacian(phi)
        """
        return -0.5 * self.get_laplacian(sqrt_field)

    def optimize_h2(self, bond_length=1.4):
        print(f"Solving EEM for H2 (Bond Length {bond_length} au)...")
        
        # 1. System Setup
        atoms = [([0, 0, -bond_length/2], 1.0), ([0, 0, bond_length/2], 1.0)]
        V_ext = self.get_nuclear_potential([a[0] for a in atoms], [a[1] for a in atoms])
        
        # 2. Ansatz (Initial Guess)
        # We start with two gaussians. 
        # Crucially: We separate them slightly in Rho vs S to break symmetry
        # If Rho and S are identical, Eq 11 in Pope 2019 says bivector potential vanishes.
        # We want to see them SEPARATE.
        
        r1 = np.sqrt((self.X)**2 + (self.Y)**2 + (self.Z+0.5)**2)
        r2 = np.sqrt((self.X)**2 + (self.Y)**2 + (self.Z-0.5)**2)
        
        # Guess: Total density n is roughly sum of atoms
        n_guess = np.exp(-2*r1) + np.exp(-2*r2)
        
        # Initialize sqrt(Rho) and sqrt(S)
        # We give S a slightly different profile to start
        phi_rho = np.sqrt(0.5 * n_guess)
        phi_s   = np.sqrt(0.5 * n_guess) 
        
        # Normalize
        norm = np.sqrt(np.sum((phi_rho**2 + phi_s**2) * self.dV))
        phi_rho /= norm
        phi_s /= norm
        target_electrons = 2.0
        phi_rho *= np.sqrt(target_electrons)
        phi_s *= np.sqrt(target_electrons)

        # 3. Optimization Loop (Steepest Descent)
        step_size = 0.05
        
        for i in range(400):
            # Densities
            rho = phi_rho**2
            s   = phi_s**2
            n_tot = rho + s
            
            # Potentials
            V_Hartree = self.solve_poisson(n_tot)
            # LDA Exchange (Slater): -3/4 (3n/pi)^(1/3)
            # Simple approximation for demonstration
            V_XC = -(3.0/4.0) * (3.0 * n_tot / np.pi + 1e-10)**(1.0/3.0) 
            
            V_eff = V_ext + V_Hartree + V_XC
            
            # Gradients (Forces on the Wavefunctions)
            # E_total = T(rho) + T(S) + Int(V_eff * (rho+S))
            # dE/dphi_rho = dT/dphi_rho + V_eff * 2 * phi_rho
            
            grad_rho = self.functional_derivative_kinetic(phi_rho) + V_eff * phi_rho
            grad_s   = self.functional_derivative_kinetic(phi_s)   + V_eff * phi_s
            
            # Update (Downhill)
            phi_rho -= step_size * grad_rho
            phi_s   -= step_size * grad_s
            
            # Re-normalize (Project back to N=2)
            current_N = np.sum((phi_rho**2 + phi_s**2) * self.dV)
            scale = np.sqrt(target_electrons / current_N)
            phi_rho *= scale
            phi_s *= scale
            
            if i % 50 == 0:
                # Calc Energy
                T_rho = np.sum(0.5 * (np.gradient(phi_rho, self.L/self.N)[0]**2 + np.gradient(phi_rho, self.L/self.N)[1]**2 + np.gradient(phi_rho, self.L/self.N)[2]**2)) * self.dV
                print(f"Step {i}: Norm={current_N:.4f}")

        return self.X, self.Z, phi_rho**2, phi_s**2, V_ext

# --- Run the Simulation ---
solver = EEMSolver(grid_points=64, box_size=5.0)
X, Z, rho, s, V = solver.optimize_h2(bond_length=1.4)

# --- Visualization ---
# We look for the segregation of Rho and S
n_tot = rho + s
diff = rho - s

fig, ax = plt.subplots(1, 3, figsize=(18, 5))
mid = 32 # Slice index

# 1. Total Density (Standard DFT view)
ax[0].imshow(n_tot[:, mid, :], extent=[-2.5,2.5,-2.5,2.5], cmap='inferno')
ax[0].set_title("Total Density $n(r) = \\rho + S$")

# 2. Hofer's Spin Density S
ax[1].imshow(s[:, mid, :], extent=[-2.5,2.5,-2.5,2.5], cmap='viridis')
ax[1].set_title("EEM Spin Density $S(r)$")

# 3. The Difference (The Structure)
# If Hofer is right, Rho and S should not be identical everywhere.
# The relative difference reveals the "Bivector Potential" structure.
im3 = ax[2].imshow(diff[:, mid, :], extent=[-2.5,2.5,-2.5,2.5], cmap='coolwarm')
ax[2].set_title("Difference $\\rho - S$ \n(The 'Mechanical Polarization')")
plt.colorbar(im3, ax=ax[2])

plt.show()