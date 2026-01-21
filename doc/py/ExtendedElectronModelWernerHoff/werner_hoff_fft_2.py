import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fftn, ifftn
from mpl_toolkits.mplot3d import Axes3D

def get_water_geometry():
    # O at center
    # H-O-H angle approx 104.5 deg
    return [
        {'pos': [0,0,0], 'Z': 8.0},
        {'pos': [1.4, 0, 1.1], 'Z': 1.0}, # Approx coords in Bohr
        {'pos': [-1.4, 0, 1.1], 'Z': 1.0}
    ]

def get_acetylene_geometry():
    return [
        {'pos': [0,0, 1.1], 'Z': 6.0},
        {'pos': [0,0,-1.1], 'Z': 6.0},
        {'pos': [0,0, 3.1], 'Z': 1.0},
        {'pos': [0,0,-3.1], 'Z': 1.0}
    ]

class EEMGeneralSolver:
    def __init__(self, box_size=8.0, grid_points=40): # Grid reduced for speed
        self.L = box_size
        self.N = grid_points
        self.dV = (self.L / self.N)**3
        
        # Real Space Grid
        x = np.linspace(-self.L/2, self.L/2, self.N, endpoint=False)
        self.X, self.Y, self.Z = np.meshgrid(x, x, x, indexing='ij')
        
        # Reciprocal (Fourier) Space Grid
        k = (2 * np.pi / self.L) * np.fft.fftfreq(self.N) * self.N
        self.KX, self.KY, self.KZ = np.meshgrid(k, k, k, indexing='ij')
        self.K2 = self.KX**2 + self.KY**2 + self.KZ**2
        self.K2[0,0,0] = 1.0 # Avoid singularity at DC component

    def fft_laplacian(self, field):
        """Computes Laplacian using FFT."""
        return np.real(ifftn(-self.K2 * fftn(field)))

    def fft_gradient(self, field):
        """Computes Gradient vector field using FFT."""
        F = fftn(field)
        gx = np.real(ifftn(1j * self.KX * F))
        gy = np.real(ifftn(1j * self.KY * F))
        gz = np.real(ifftn(1j * self.KZ * F))
        return gx, gy, gz

    def solve_poisson(self, density):
        """Solves Poisson Equation: -Del^2 V = 4*pi*n"""
        rho_k = fftn(density)
        # V(k) = 4 * pi * rho(k) / k^2
        V_k = 4 * np.pi * rho_k / self.K2
        V_k[0,0,0] = 0.0 # DC component
        return np.real(ifftn(V_k))

    def get_nuclear_potential(self, molecule_def):
        """Generates smoothed nuclear potential for arbitrary geometry."""
        V = np.zeros_like(self.X)
        sigma = 0.5 # Gaussian smearing width (Pseudopotential-like)
        
        print("System Geometry:")
        for atom in molecule_def:
            pos = atom['pos']
            Z = atom['Z']
            print(f"  Atom Z={Z} at {pos}")
            # Distance squared
            r2 = (self.X - pos[0])**2 + (self.Y - pos[1])**2 + (self.Z - pos[2])**2
            # Gaussian nuclear charge distribution
            V += -Z * (1.0 / np.sqrt(r2 + sigma**2))
            
        return V

    def functional_derivative_kinetic(self, psi):
        """
        Derivative of Weizsacker Kinetic Energy w.r.t wavefunction psi = sqrt(rho)
        dT/dpsi = -1/2 * Laplacian(psi)
        """
        return -0.5 * self.fft_laplacian(psi)

    def solve(self, molecule_def, n_electrons, steps=300, alpha=0.02):
        print(f"\nStarting EEM Minimization for {n_electrons} electrons...")
        
        # 1. Initialize Potentials
        V_ext = self.get_nuclear_potential(molecule_def)
        
        # 2. Initialize Wavefunctions (sqrt(rho) and sqrt(S))
        # Ansatz: Sum of Gaussians on atoms + Random noise to break symmetry
        psi_guess = np.zeros_like(self.X)
        for atom in molecule_def:
            pos = atom['pos']
            r2 = (self.X - pos[0])**2 + (self.Y - pos[1])**2 + (self.Z - pos[2])**2
            psi_guess += atom['Z'] * np.exp(-1.5 * r2)
            
        # IMPORTANT: Initialize Rho and S slightly differently.
        # If they are identical, the Hofer coupling term vanishes.
        # We add "Topological Noise" to S to seed the filaments.
        psi_rho = np.sqrt(psi_guess)
        psi_s   = np.sqrt(psi_guess) * (1.0 + 0.1 * np.cos(3*self.X)*np.sin(3*self.Y))
        
        # Normalize
        psi_rho = self.normalize(psi_rho, n_electrons)
        psi_s   = self.normalize(psi_s, n_electrons)

        # 3. Minimization Loop (Steepest Descent)
        for i in range(steps):
            rho = psi_rho**2
            S   = psi_s**2
            n_tot = rho + S # Total particle density
            
            # --- Potentials ---
            V_Hartree = self.solve_poisson(n_tot)
            
            # EEM Bivector Coupling (Approximated)
            # In the 2019 paper, the potential v_b depends on the mismatch of Kinetic Densities.
            # This pushes S away from Rho.
            # We model this phenomenologically as a repulsive force between the two fluids
            # proportional to the local density.
            V_coupling_rho = 0.5 * S / (rho + S + 1e-6)
            V_coupling_s   = 0.5 * rho / (rho + S + 1e-6)
            
            V_eff_rho = V_ext + V_Hartree + V_coupling_rho
            V_eff_s   = V_ext + V_Hartree + V_coupling_s
            
            # --- Gradients ---
            # Gradient = dT/dpsi + V_eff * psi
            grad_rho = self.functional_derivative_kinetic(psi_rho) + V_eff_rho * psi_rho
            grad_s   = self.functional_derivative_kinetic(psi_s)   + V_eff_s * psi_s
            
            # --- Update ---
            psi_rho -= alpha * grad_rho
            psi_s   -= alpha * grad_s
            
            # --- Renormalize ---
            psi_rho = self.normalize(psi_rho, n_electrons)
            psi_s   = self.normalize(psi_s, n_electrons)
            
            if i % 50 == 0:
                print(f"  Step {i}")

        return psi_rho**2, psi_s**2

    def normalize(self, psi, target_N):
        current_N = np.sum(psi**2) * self.dV
        return psi * np.sqrt(target_N / current_N)

# --- Define Molecules ---

def get_ethylene_geometry():
    # C=C bond length approx 1.34 Angstrom (approx 2.5 Bohr)
    # H-C-H angle 120 deg
    cc_dist = 2.5
    ch_dist = 2.0
    
    # C atoms on Z axis
    c1 = [0, 0, cc_dist/2]
    c2 = [0, 0, -cc_dist/2]
    
    # H atoms in XZ plane (planar molecule)
    # cos(120) = -0.5, sin(120) = 0.866
    h1 = [ch_dist * 0.866, 0,  cc_dist/2 + ch_dist * 0.5]
    h2 = [-ch_dist * 0.866, 0, cc_dist/2 + ch_dist * 0.5]
    h3 = [ch_dist * 0.866, 0, -cc_dist/2 - ch_dist * 0.5]
    h4 = [-ch_dist * 0.866, 0, -cc_dist/2 - ch_dist * 0.5]
    
    return [
        {'pos': c1, 'Z': 6.0}, {'pos': c2, 'Z': 6.0},
        {'pos': h1, 'Z': 1.0}, {'pos': h2, 'Z': 1.0},
        {'pos': h3, 'Z': 1.0}, {'pos': h4, 'Z': 1.0}
    ]

# --- Main Execution ---

# 1. Setup Solver
# Note: box_size needs to be big enough to hold the molecule
solver = EEMGeneralSolver(box_size=12.0, grid_points=64) 

# 2. Select System: Ethylene (C2H4)
# Total electrons = 6*2 + 4*1 = 16
molecule = get_ethylene_geometry()
rho, S = solver.solve(molecule, n_electrons=16.0, steps=250)

# 3. Calculate "Hofer Vorticity" (Filaments)
# Filament = | Grad(Rho) x Grad(S) |
grx, gry, grz = solver.fft_gradient(rho)
gsx, gsy, gsz = solver.fft_gradient(S)

# Cross Product
vx = gry*gsz - grz*gsy
vy = grz*gsx - grx*gsz
vz = grx*gsy - gry*gsx
vorticity_mag = np.sqrt(vx**2 + vy**2 + vz**2)

# --- Visualization ---

# Slice index (middle of box)
mid = solver.N // 2

fig = plt.figure(figsize=(15, 5))

# Plot 1: Standard Density (C=C bond vertical)
ax1 = fig.add_subplot(131)
ax1.set_title("Total Density (C=C vertical)")
n_tot = rho + S
# Slice through Y=0 (Plane of the molecule)
im1 = ax1.imshow(n_tot[:, mid, :], extent=[-6,6,-6,6], cmap='inferno', origin='lower')
plt.colorbar(im1, ax=ax1)

# Plot 2: Mechanical Polarization (Rho - S)
# This shows where Mass and Spin try to separate
ax2 = fig.add_subplot(132)
ax2.set_title("Mechanical Polarization\n$(\\rho - S)$")
diff = rho - S
im2 = ax2.imshow(diff[:, mid, :], extent=[-6,6,-6,6], cmap='coolwarm', origin='lower')
plt.colorbar(im2, ax=ax2)

# Plot 3: The FILAMENTS (Vorticity)
# Looking for topological defects
ax3 = fig.add_subplot(133)
ax3.set_title("EEM Vortex Filaments\n$|\\nabla \\rho \\times \\nabla S|$")
im3 = ax3.imshow(vorticity_mag[:, mid, :], extent=[-6,6,-6,6], cmap='cubehelix', origin='lower')
plt.colorbar(im3, ax=ax3)

plt.tight_layout()
plt.show()

# Optional: 3D Isosurface of the Filaments
# If you run this locally, this generates a pop-up 3D plot
try:
    from skimage import measure
    verts, faces, _, _ = measure.marching_cubes(vorticity_mag, level=np.max(vorticity_mag)*0.3)
    
    fig_3d = plt.figure(figsize=(8,8))
    ax_3d = fig_3d.add_subplot(111, projection='3d')
    ax_3d.plot_trisurf(verts[:, 0], verts[:, 1], faces, verts[:, 2], color='cyan', alpha=0.5)
    ax_3d.set_title("3D Structure of EEM Spin-Vortices (Ethylene)")
    plt.show()
except ImportError:
    print("skimage not installed, skipping 3D isosurface plot.")