import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# --- CONFIGURATION: The "Chemistry" Parameters ---
# Adjust these energies (eigenvalues) to change the beat frequencies
E_SIGMA_CC     = 0.0   # Ground state (C-C bonding)
E_SIGMA_CC_STAR= 2.0   # Antibonding (C-C*)
E_SIGMA_CH     = 0.5   # C-H bonding energy
E_SIGMA_CH_STAR= 2.5   # C-H Antibonding

# Simulation Speed
TIME_STEP = 0.05
FRAMES = 300

class ChemicalCanvas:
    def __init__(self, res=200, bounds=5.0):
        # 1. Setup Grid
        x = np.linspace(-bounds, bounds, res)
        y = np.linspace(-bounds, bounds, res)
        self.X, self.Y = np.meshgrid(x, y)
        self.shape = self.X.shape
        
        # 2. Define Atoms (C2H4 geometry approx)
        # Bond length ~ 2.5 Bohr. Atoms at -1.25 and +1.25
        self.c1_pos = np.array([-1.25, 0.0]) 
        self.c2_pos = np.array([ 1.25, 0.0])
        
        # Precompute Basis Functions (STO-3G style) to save FPS
        # We need s, px, py on both atoms
        self.basis = {}
        self._precompute_basis(self.c1_pos, 'C1')
        self._precompute_basis(self.c2_pos, 'C2')
        
    def _precompute_basis(self, center, label):
        # Radial distance
        R = np.sqrt((self.X - center[0])**2 + (self.Y - center[1])**2)
        # Radial decay (alpha ~ 1.5 for Carbon valence)
        radial = np.exp(-1.5 * R)
        
        # Angular parts (unnormalized but qualitative)
        self.basis[f'{label}_s']  = radial
        self.basis[f'{label}_px'] = (self.X - center[0]) * radial
        self.basis[f'{label}_py'] = (self.Y - center[1]) * radial

    def get_sp2_hybrid(self, atom_label, angle_deg):
        """
        Constructs a directed sp2 orbital at a specific angle.
        angle 0 points +x, 180 points -x.
        """
        rad = np.radians(angle_deg)
        c_s  = 1.0 / np.sqrt(3)
        c_p  = np.sqrt(2.0/3.0)
        
        # sp2 = c_s*s + c_p*(cos(th)*px + sin(th)*py)
        s  = self.basis[f'{atom_label}_s']
        px = self.basis[f'{atom_label}_px']
        py = self.basis[f'{atom_label}_py']
        
        return c_s * s + c_p * (np.cos(rad) * px + np.sin(rad) * py)

# --- Initialize System ---
sim = ChemicalCanvas(res=300, bounds=6.0)

# --- 1. Construct Molecular Orbitals (Spatial Part) ---

# C1 Hybrids (Left Atom)
# Pointing Right (towards C2)
phi_L_in  = sim.get_sp2_hybrid('C1', 0)    
# Pointing Top-Left (towards H)
phi_L_up  = sim.get_sp2_hybrid('C1', 120)  
# Pointing Bottom-Left
phi_L_dn  = sim.get_sp2_hybrid('C1', 240)  

# C2 Hybrids (Right Atom)
# Pointing Left (towards C1)
phi_R_in  = sim.get_sp2_hybrid('C2', 180)  
# Pointing Top-Right
phi_R_up  = sim.get_sp2_hybrid('C2', 60)   
# Pointing Bottom-Right
phi_R_dn  = sim.get_sp2_hybrid('C2', 300)  

# --- 2. Build Eigenstates (Linear Combinations) ---

# State A: C-C Sigma Bonding (Symmetric)
# Energy: Low
psi_CC_bond = (phi_L_in + phi_R_in)
psi_CC_bond /= np.max(np.abs(psi_CC_bond)) # Normalize max to 1

# State B: C-C Antibonding (Asymmetric / Node in center)
# Energy: High
psi_CC_anti = (phi_L_in - phi_R_in)
psi_CC_anti /= np.max(np.abs(psi_CC_anti))

# State C: C-H Skeleton (Symmetric sum of all H-pointing hybrids)
psi_CH_skel = (phi_L_up + phi_L_dn + phi_R_up + phi_R_dn)
psi_CH_skel /= np.max(np.abs(psi_CH_skel))

# --- 3. Define the Superposition (The Animation) ---

# Choose what to mix! 
# Here we mix the Bonding and Antibonding C-C states plus the C-H skeleton.
# This creates a "Sloshing" of charge between the carbons.

active_states = [
    # (Spatial_Wavefunction, Energy_Eigenvalue, Amplitude_Coefficient)
    
    # 1. The Main Bond (Stationary Base)
    {'psi': psi_CC_bond, 'E': E_SIGMA_CC,      'amp': 1.0},
    
    # 2. The Antibonding (Adds the "Beating/Sloshing")
    # Because E is different, the relative phase rotates, causing interference.
    {'psi': psi_CC_anti, 'E': E_SIGMA_CC_STAR, 'amp': 0.8},
    
    # 3. The C-H frame (Just to make it look like Ethylene)
    {'psi': psi_CH_skel, 'E': E_SIGMA_CH,      'amp': 0.5}
]

# --- Visualization Setup ---
fig, ax = plt.subplots(figsize=(8, 6), facecolor='black')
ax.set_facecolor('black')

# Initial Calculation
rho = np.zeros(sim.shape)
im = ax.imshow(rho, cmap='magma', origin='lower', animated=True, 
               extent=[-6,6,-6,6], vmin=0, vmax=2.5)

# Plot Atom Positions markers
ax.scatter([sim.c1_pos[0], sim.c2_pos[0]], [0,0], color='cyan', marker='o', s=50, alpha=0.5)
ax.text(sim.c1_pos[0], -0.5, "C1", color='cyan', ha='center')
ax.text(sim.c2_pos[0], -0.5, "C2", color='cyan', ha='center')

# Labels
time_text = ax.text(0.05, 0.95, '', transform=ax.transAxes, color='white', fontsize=12)
ax.set_title(f"Dynamic Density Evolution (Energy-Driven)", color='white')

def update(frame):
    t = frame * TIME_STEP
    
    # Construct total complex wavefunction
    # Psi(t) = Sum [ c_n * psi_n(r) * exp(-i * E_n * t) ]
    psi_total = np.zeros(sim.shape, dtype=complex)
    
    for state in active_states:
        phase = np.exp(-1j * state['E'] * t)
        psi_total += state['amp'] * state['psi'] * phase
        
    # Density = |Psi|^2
    density = np.abs(psi_total)**2
    
    # Update Image
    im.set_data(density)
    time_text.set_text(f"Time: {t:.2f} au")
    
    return im, time_text

ani = FuncAnimation(fig, update, frames=FRAMES, interval=30, blit=True)
plt.tight_layout()
plt.show()