import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fftn, ifftn, fftfreq
from mpl_toolkits.mplot3d import art3d

class HoferSmokeSolver:
    def __init__(self, box_size=10.0, grid_points=64):
        """
        Initializes the Spectral Solver Grid.
        """
        self.L = box_size
        self.N = grid_points
        self.dx = self.L / self.N
        self.dV = self.dx**3
        
        # 1. Real Space Grid
        x = np.linspace(-self.L/2, self.L/2, self.N, endpoint=False)
        self.X, self.Y, self.Z = np.meshgrid(x, x, x, indexing='ij')
        
        # 2. Fourier Space Grid (for Kinetic Energy operator)
        k = (2 * np.pi / self.L) * fftfreq(self.N, d=1/self.N)
        self.KX, self.KY, self.KZ = np.meshgrid(k, k, k, indexing='ij')
        self.K2 = self.KX**2 + self.KY**2 + self.KZ**2
        
        # Avoid division by zero for Poisson solve
        self.K2_poisson = self.K2.copy()
        self.K2_poisson[0,0,0] = 1.0 

    def get_nuclear_potential(self, atom_list):
        """
        Creates the static external potential V_ext from nuclei.
        Uses Soft-Coulomb to avoid singularities on the grid.
        """
        V = np.zeros_like(self.X)
        softening = 0.5 # Width of the nucleus (Pseudo-potential)
        
        print("Building Molecular Potential...")
        for atom in atom_list:
            pos = atom['pos']
            charge = atom['Z']
            r = np.sqrt((self.X - pos[0])**2 + (self.Y - pos[1])**2 + (self.Z - pos[2])**2)
            # V = -Z / r
            V += -charge / np.sqrt(r**2 + softening**2)
            
        return V

    def solve_hartree(self, density):
        """
        Solves Poisson Equation for Electron-Electron Repulsion.
        Del^2 V_H = -4 * pi * density
        """
        rho_k = fftn(density)
        # V_k = 4 * pi * rho_k / k^2
        V_k = 4 * np.pi * rho_k / self.K2_poisson
        V_k[0,0,0] = 0.0 # DC component (neutral background)
        return np.real(ifftn(V_k))

    def imaginary_time_propagation(self, psi, V_ext, n_electrons, steps=200, dt=0.05):
        """
        Relaxes the wavefunction to the Ground State using Split-Step FFT.
        Equation: dPsi/dt = -(T + V)Psi  (Imaginary Time Schrodinger)
        """
        print(f"Relaxing electronic structure ({steps} steps)...")
        
        # Pre-compute Kinetic Propagator (in k-space)
        # T = k^2 / 2.  Propagator = exp(-T * dt)
        kinetic_op = np.exp(-0.5 * self.K2 * dt)
        
        for i in range(steps):
            # 1. Density
            density = np.abs(psi)**2
            
            # 2. Potential Step (Half)
            # Calculate Hartree (Repulsion)
            V_hartree = self.solve_hartree(density)
            V_total = V_ext + V_hartree
            
            # Apply Potential in Real Space: exp(-V * dt/2)
            psi = psi * np.exp(-0.5 * V_total * dt)
            
            # 3. Kinetic Step (Full)
            # FFT -> Apply Kinetic -> IFFT
            psi_k = fftn(psi)
            psi_k = psi_k * kinetic_op
            psi = ifftn(psi_k)
            
            # 4. Potential Step (Half)
            # Re-evaluate potential (optional for higher accuracy, usually skip in simple codes)
            # reusing previous V_total for speed
            psi = psi * np.exp(-0.5 * V_total * dt)
            
            # 5. Renormalize (Chemistry Constraint)
            # Ensure total electron number is conserved
            current_N = np.sum(np.abs(psi)**2) * self.dV
            psi *= np.sqrt(n_electrons / current_N)
            
            if i % 50 == 0:
                E_pot = np.sum(density * V_total) * self.dV
                print(f"  Step {i}: Norm={current_N:.4f}, PotEnergy={E_pot:.2f}")
                
        return psi

# --- Visualization Helpers ---
def plot_isosurfaces(X, Y, Z, rho, psi_phase, title):
    """
    Visualizes the filament structure.
    """
    from skimage import measure
    
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # 1. Plot Density Isosurface (The "Cloud")
    # Low opacity to see inside
    try:
        verts, faces, _, _ = measure.marching_cubes(rho, level=np.max(rho)*0.2)
        # Transform grid indices to physical coordinates
        verts_phys = verts * (10.0/64.0) - 5.0 
        mesh = art3d.Poly3DCollection(verts_phys[faces], alpha=0.1, color='blue')
        ax.add_collection3d(mesh)
    except:
        pass

    # 2. Plot THE FILAMENT (The Node)
    # The filament is where density is very low, but surrounded by high density.
    # Or structurally: where the phase flips.
    # We plot a very low density isosurface to find the "holes" in the orbital.
    try:
        # Level set close to zero
        verts_node, faces_node, _, _ = measure.marching_cubes(rho, level=np.max(rho)*0.01)
        verts_node_phys = verts_node * (10.0/64.0) - 5.0
        
        # Only plot the node if it's "inside" the molecule bounds
        mesh_node = art3d.Poly3DCollection(verts_node_phys[faces_node], alpha=0.4, color='red')
        ax.add_collection3d(mesh_node)
    except:
        print("Could not extract filament isosurface.")

    ax.set_xlim(-5, 5)
    ax.set_ylim(-5, 5)
    ax.set_zlim(-5, 5)
    ax.set_title(title)
    ax.set_xlabel("X (Bohr)")
    ax.set_ylabel("Y (Bohr)")
    ax.set_zlabel("Z (Bohr)")
    
    # Legend trick
    import matplotlib.patches as mpatches
    blue_patch = mpatches.Patch(color='blue', alpha=0.2, label='Electron Density')
    red_patch = mpatches.Patch(color='red', label='Vortex Filament (Node)')
    plt.legend(handles=[blue_patch, red_patch])
    
    plt.show()

# --- Main Simulation: Ethylene ---

if __name__ == "__main__":
    solver = HoferSmokeSolver(box_size=10.0, grid_points=64)
    
    # 1. Define Ethylene (C2H4) Geometry
    # C=C bond along Z axis. Molecule in XZ plane.
    # This geometry forces the Pi-bond to have a node in the YZ plane (X=0).
    cc_dist = 2.5
    h_dist = 2.0
    theta = np.deg2rad(120)
    
    atoms = [
        {'pos': [0, 0,  cc_dist/2], 'Z': 6.0}, # C1
        {'pos': [0, 0, -cc_dist/2], 'Z': 6.0}, # C2
        # Hydrogens
        {'pos': [h_dist*np.sin(theta), 0, cc_dist/2 + h_dist*np.cos(theta)], 'Z': 1.0},
        {'pos': [-h_dist*np.sin(theta), 0, cc_dist/2 + h_dist*np.cos(theta)], 'Z': 1.0},
        {'pos': [h_dist*np.sin(theta), 0, -cc_dist/2 - h_dist*np.cos(theta)], 'Z': 1.0},
        {'pos': [-h_dist*np.sin(theta), 0, -cc_dist/2 - h_dist*np.cos(theta)], 'Z': 1.0},
    ]
    
    # 2. Get Potential
    V_ext = solver.get_nuclear_potential(atoms)
    
    # 3. Initialize Wavefunction (The "Pi" Ansatz)
    # To see the filament, we must look at the Pi-bond electrons.
    # We initialize with a "p-orbital" shape (antisymmetric along X).
    # If we started symmetric, we would just get the boring Sigma bond.
    
    # Pi-orbital ansatz: x * exp(-r)
    # We add a small complex component (i*y) to give it "vorticity" or "spin"
    # This turns the planar node into a vortex core.
    psi_guess = (solver.X + 0.1j*solver.Y) * np.exp(-1.0 * (solver.X**2 + solver.Y**2 + solver.Z**2))
    
    # Normalize initial guess
    psi_guess /= np.sqrt(np.sum(np.abs(psi_guess)**2) * solver.dV)
    psi_guess *= np.sqrt(2.0) # 2 electrons in the Pi bond
    
    # 4. Run "Hofer Smoke" Relaxation
    # This corresponds to finding the ground state of the EEM Hamiltonian.
    psi_final = solver.imaginary_time_propagation(psi_guess, V_ext, n_electrons=2.0, steps=150)
    
    # 5. Analyze Results
    rho = np.abs(psi_final)**2
    phase = np.angle(psi_final)
    
    # Plot Slices
    mid = solver.N // 2
    
    plt.figure(figsize=(12, 5))
    
    # Plot A: The Density (Showing the lobes)
    plt.subplot(1, 2, 1)
    # Slice through XZ plane (Side view of molecule)
    plt.imshow(rho[:, mid, :], extent=[-5,5,-5,5], cmap='inferno', origin='lower')
    plt.title("Electron Density $\\rho$ (Side View)\nNote the gap in the middle")
    plt.colorbar()
    
    # Plot B: The Phase (Showing the Vortex/Node)
    plt.subplot(1, 2, 2)
    # Slice through XY plane (Cross section of the bond)
    # This should look like a vortex (phase winding)
    plt.imshow(phase[:, :, mid], extent=[-5,5,-5,5], cmap='twilight', origin='lower')
    plt.title("Hofer Phase Structure $\\phi$\n(Cross Section of C-C bond)")
    plt.colorbar()
    
    plt.show()
    
    print("Generating 3D Filament Plot...")
    plot_isosurfaces(solver.X, solver.Y, solver.Z, rho, phase, "Ethylene Pi-Bond: Density & Vortex Filament")