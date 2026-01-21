import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# --- CONFIGURATION ---
# We split the orbitals to preserve the Y-axis geometry (Triangular shape)

# State 1: The Central C-C Bond
E_CC    = 1.0  

# State 2: The UPPER C-H Bonds (Top Left + Top Right)
E_UPPER = 1.5  

# State 3: The LOWER C-H Bonds (Bot Left + Bot Right)
E_LOWER = 2.0  

# Mixing Coefficients (Try changing these!)
c_cc    = 1.0        # Central Bond
c_upper = 0.6        # Top 'V' shape
c_lower = 0.6 * 1.0j # Bot 'V' shape (Imaginary -> 90 deg phase shift)

# Simulation Speed
TIME_STEP = 0.1
FRAMES = 200

# ---------------------

class Sp2Canvas:
    def __init__(self, res=300, bounds=5.0):
        # 1. Grid (XY Plane - Molecule is flat on screen)
        x = np.linspace(-bounds, bounds, res)
        y = np.linspace(-bounds, bounds, res)
        self.X, self.Y = np.meshgrid(x, y)
        self.dV = (x[1]-x[0])**2
        
        # 2. Atoms positions (C1 Left, C2 Right)
        self.c1 = np.array([-1.3, 0.0])
        self.c2 = np.array([ 1.3, 0.0])
        
        # 3. Precompute Basis (Slater Type Orbitals)
        self.basis = {}
        self._make_basis(self.c1, 'C1')
        self._make_basis(self.c2, 'C2')
        
    def _make_basis(self, pos, label):
        # Distance from atom
        R = np.sqrt((self.X - pos[0])**2 + (self.Y - pos[1])**2)
        # Decay (slower decay to make lobes visible)
        radial = np.exp(-1.2 * R)
        
        # Atomic Orbitals
        self.basis[f'{label}_s']  = radial
        self.basis[f'{label}_px'] = (self.X - pos[0]) * radial
        self.basis[f'{label}_py'] = (self.Y - pos[1]) * radial

    def get_lobe(self, atom, angle_deg):
        """ Creates one sp2 lobe pointing at angle_deg """
        rad = np.radians(angle_deg)
        # sp2 mixing: s + sqrt(2)*p
        s = self.basis[f'{atom}_s']
        px = self.basis[f'{atom}_px']
        py = self.basis[f'{atom}_py']
        
        # Coefficients for sp2 hybridization
        # s_coeff = 1/sqrt(3), p_coeff = sqrt(2/3)
        # We boost p_coeff slightly to make the "bunny ears" sharper for vis
        psi = 0.5 * s + 0.8 * (np.cos(rad)*px + np.sin(rad)*py)
        return psi

# Setup Simulation
sim = Sp2Canvas(res=256, bounds=5.0)

# --- CONSTRUCT THE GEOMETRY ---

# 1. Central C-C Bond (Sigma)
# C1 points right (0), C2 points left (180)
phi_cc = sim.get_lobe('C1', 0) + sim.get_lobe('C2', 180)

# 2. Upper "V" Shape (120 deg geometry)
# C1 Top-Left (120), C2 Top-Right (60)
phi_upper = sim.get_lobe('C1', 120) + sim.get_lobe('C2', 60)

# 3. Lower "V" Shape (120 deg geometry)
# C1 Bot-Left (240), C2 Bot-Right (300)
phi_lower = sim.get_lobe('C1', 240) + sim.get_lobe('C2', 300)

# Normalize for brightness consistency
def norm(psi): return psi / np.sqrt(np.sum(np.abs(psi)**2)*sim.dV)

phi_cc    = norm(phi_cc)
phi_upper = norm(phi_upper)
phi_lower = norm(phi_lower)

# --- ANIMATION ---

fig, ax = plt.subplots(figsize=(8, 6), facecolor='black')
ax.set_facecolor('black')

# Initial frame
rho = np.abs(phi_cc)**2
im = ax.imshow(rho, cmap='gist_heat', origin='lower', animated=True,
               extent=[-5,5,-5,5], vmin=0, vmax=0.4)

# Draw Skeleton (Cyan lines) to show ideal geometry
# C1
ax.plot([-1.3, -1.3 + np.cos(np.deg2rad(120))*2], [0, np.sin(np.deg2rad(120))*2], 'c--', alpha=0.3)
ax.plot([-1.3, -1.3 + np.cos(np.deg2rad(240))*2], [0, np.sin(np.deg2rad(240))*2], 'c--', alpha=0.3)
# C2
ax.plot([1.3, 1.3 + np.cos(np.deg2rad(60))*2], [0, np.sin(np.deg2rad(60))*2], 'c--', alpha=0.3)
ax.plot([1.3, 1.3 + np.cos(np.deg2rad(300))*2], [0, np.sin(np.deg2rad(300))*2], 'c--', alpha=0.3)
# Atoms
ax.scatter([-1.3, 1.3], [0, 0], color='cyan', s=50)

title = ax.set_title("Ethylene sp2 Density Dynamics", color='white')

def update(frame):
    t = frame * TIME_STEP
    
    # Superposition of states with different energies
    # This causes the density to "beat" or flow between the lobes
    psi = (c_cc    * phi_cc    * np.exp(-1j * E_CC * t) +
           c_upper * phi_upper * np.exp(-1j * E_UPPER * t) +
           c_lower * phi_lower * np.exp(-1j * E_LOWER * t))
    
    rho = np.abs(psi)**2
    im.set_data(rho)
    title.set_text(f"t={t:.2f} | 120-degree Geometry Visible")
    return im, title

ani = FuncAnimation(fig, update, frames=FRAMES, interval=40, blit=True)
plt.show()